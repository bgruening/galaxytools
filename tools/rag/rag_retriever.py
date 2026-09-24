import hashlib
import json
import os
import sys
import time
from functools import partial
from pathlib import Path

import httpx
import torch
import yaml
from llama_index.core import Document, SimpleDirectoryReader, VectorStoreIndex
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.embeddings.openai import OpenAIEmbedding
from llama_index.readers.file import PDFReader
from llama_index.readers.json import JSONReader
from sentence_transformers import CrossEncoder


# How many chunks a reranker re-scores when the form does not say. The
# galaxy-rag-project evaluation found gains up to about 20 candidates and none
# beyond.
DEFAULT_RERANK_CANDIDATES = 20
# Cross-encoders judge the question and one chunk together; 512 tokens is the
# window the evaluated models (ms-marco-MiniLM, MedCPT, bge-reranker) were
# trained on.
RERANK_MAX_LENGTH = 512


# --- LiteLLM proxy config resolution -------------------------------------
# The LiteLLM proxy exposes an OpenAI-compatible /v1/embeddings endpoint. The
# YAML config may be flat (global LITELLM_API_KEY / LITELLM_BASE_URL) or expose
# a ``servers`` mapping that keys provider names to per-server credentials; the
# ``provider`` argument selects which server to use.

def load_litellm_config() -> dict:
    """Read the LiteLLM YAML config referenced by ``LITELLM_CONFIG_FILE``.

    Exits with a clear message if the env var is unset or the file is missing,
    since neither is recoverable for a tool job.
    """
    config_file = os.environ.get("LITELLM_CONFIG_FILE")
    if not config_file:
        sys.exit("LITELLM_CONFIG_FILE environment variable is not set.")
    if not os.path.isfile(config_file):
        sys.exit(f"LiteLLM config file does not exist: {config_file}")
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
    if not config:
        sys.exit(
            f"LiteLLM config file is empty or contains no entries: {config_file}"
        )
    return config


def resolve_server(config: dict, provider: str) -> dict:
    """Resolve the server config for ``provider`` and validate its credentials.

    Returns the per-provider ``servers[provider]`` dict when a ``servers`` block
    exists, otherwise the global config (backward compatibility). Exits with a
    message if the provider is unknown or the API key / base URL is missing.
    """
    servers = config.get("servers", {})
    if servers:
        if provider not in servers:
            sys.exit(f"Provider '{provider}' not found in LiteLLM configuration.")
        source = servers[provider]
    else:
        source = config
    if not source.get("LITELLM_API_KEY"):
        sys.exit(
            "LiteLLM API key is not configured! Please set LITELLM_API_KEY "
            "in the configuration."
        )
    if not source.get("LITELLM_BASE_URL"):
        sys.exit(
            "LiteLLM base URL is not configured! Please set LITELLM_BASE_URL "
            "in the configuration."
        )
    return source


def attribution_user(galaxy_user_id: str, galaxy_url: str) -> str:
    """Proxy-side identity for the Galaxy user, sent as the OpenAI ``user`` field.

    Proxies (e.g. LiteLLM) use it to meter usage and apply per-user budgets and
    rate limits. The id is namespaced by the Galaxy instance URL and hashed: the
    URL keeps ids unique when several Galaxy instances share one proxy (e.g.
    usegalaxy.eu), and hashing means no instance-identifying or personal data
    leaves for the provider. Anonymous users (no id) fall back to a
    per-instance shared "anonymous" bucket. Mirrors the attribution in
    llm_hub.py so both tools map one Galaxy user to one proxy identity.
    """
    raw_user = galaxy_user_id or "anonymous"
    return hashlib.sha256(f"{galaxy_url}|{raw_user}".encode()).hexdigest()


def build_embed_model(embed_cfg: dict, user: str):
    """Embedding model for the LiteLLM proxy or a local HuggingFace model."""
    if embed_cfg["source"] == "litellm":
        model = embed_cfg.get("model")
        provider = embed_cfg.get("provider")
        if not model:
            sys.exit("No LiteLLM embedding model selected.")
        if not provider:
            sys.exit("No LiteLLM provider selected.")
        server = resolve_server(load_litellm_config(), provider)
        # OpenAIEmbedding validates ``model`` against a hardcoded enum of OpenAI
        # model names; passing the LiteLLM model id via ``model_name`` instead is
        # the class's intended escape hatch (it overrides the enum-derived
        # engine), so arbitrary proxy-hosted models (BGE-M3, nomic, Qwen3, ...)
        # can be used with the framework's batching, retry and client handling.
        # ``additional_kwargs`` is forwarded as ``**kwargs`` to the OpenAI SDK's
        # ``embeddings.create`` call, carrying the ``user`` field through.
        return OpenAIEmbedding(
            model_name=model,
            api_key=server["LITELLM_API_KEY"],
            api_base=server["LITELLM_BASE_URL"],
            timeout=float(os.environ.get("LITELLM_REQUEST_TIMEOUT", "600")),
            max_retries=int(os.environ.get("LITELLM_REQUEST_MAX_RETRIES", "3")),
            embed_batch_size=100,
            additional_kwargs={"user": user},
        )
    # Local HuggingFace model (preinstalled path or uploaded archive).
    model_path = embed_cfg.get("path")
    if not model_path:
        sys.exit("No embedding model path given.")
    if not os.path.exists(model_path):
        sys.exit(f"Embedding model path does not exist: {model_path}")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    return HuggingFaceEmbedding(model_name=model_path, normalize=True, device=device)


# --- Reranking --------------------------------------------------------------
# Retrieval embeds the question and each chunk separately and compares the two
# vectors. A cross-encoder reads the question and one chunk together, so it can
# judge whether the chunk answers *this* question: better, but one model pass
# per pair. So retrieval first narrows the corpus to a few candidates, the
# reranker re-sorts them, and the best ``top_k`` are kept.

def rerank_local(model_path: str, question: str, chunks: list, top_k: int) -> list:
    """Re-sort ``chunks`` with a local cross-encoder; return the best ``top_k``."""
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = CrossEncoder(model_path, max_length=RERANK_MAX_LENGTH, device=device)
    scores = model.predict(
        [(question, chunk) for chunk in chunks], batch_size=16, show_progress_bar=False
    )
    order = sorted(range(len(chunks)), key=lambda i: -float(scores[i]))
    return [chunks[i] for i in order[:top_k]]


def rerank_litellm(server: dict, model: str, user: str, question: str, chunks: list, top_k: int) -> list:
    """Re-sort ``chunks`` with a proxy-hosted reranker via LiteLLM's /v1/rerank.

    Transport errors, rate limits (429) and 5xx answers are retried with
    backoff, like the embedding requests; any other HTTP error ends the job
    with the proxy's message, since retrying a bad request cannot succeed.
    """
    url = server["LITELLM_BASE_URL"].rstrip("/") + "/v1/rerank"
    payload = {"model": model, "query": question, "documents": chunks, "top_n": top_k, "user": user}
    headers = {"Authorization": f"Bearer {server['LITELLM_API_KEY']}"}
    timeout = float(os.environ.get("LITELLM_REQUEST_TIMEOUT", "600"))
    max_retries = max(0, int(os.environ.get("LITELLM_REQUEST_MAX_RETRIES", "3")))
    for attempt in range(max_retries + 1):
        try:
            response = httpx.post(url, json=payload, headers=headers, timeout=timeout)
        except httpx.TransportError as e:
            error = f"{type(e).__name__}: {e}"
        else:
            if response.status_code < 400:
                break
            error = f"HTTP {response.status_code}: {response.text[:500]}"
            if response.status_code < 500 and response.status_code != 429:
                sys.exit(f"Reranker request to model '{model}' failed with {error}")
        if attempt == max_retries:
            sys.exit(f"Reranker request to model '{model}' failed after {attempt + 1} attempts: {error}")
        time.sleep(min(2 ** attempt, 30))
    try:
        results = response.json()["results"]
        order = [r["index"] for r in sorted(results, key=lambda r: -float(r["relevance_score"]))]
        return [chunks[i] for i in order[:top_k]]
    except (ValueError, KeyError, TypeError, IndexError):
        sys.exit(f"Unexpected response from reranker model '{model}': {response.text[:500]}")


def prepare_reranker(rerank_cfg: dict, user: str):
    """Validate the reranker settings before any document is embedded.

    Returns ``(rerank, candidates)``: a function ``rerank(question, chunks,
    top_k)`` and how many chunks to retrieve for it, or ``(None, 0)`` when
    reranking is off. Bad settings exit here, so a job fails in seconds rather
    than after embedding the whole corpus.
    """
    source = rerank_cfg["source"]
    if source == "none":
        return None, 0
    if source not in ("litellm", "local"):
        sys.exit(f"Unknown reranker source: {source}")
    try:
        candidates = int(rerank_cfg.get("candidates", DEFAULT_RERANK_CANDIDATES))
    except (TypeError, ValueError):
        candidates = 0
    if candidates < 1:
        sys.exit("Reranker candidates must be a positive integer.")
    if source == "litellm":
        model = rerank_cfg.get("model")
        provider = rerank_cfg.get("provider")
        if not model:
            sys.exit("No LiteLLM reranker model selected.")
        if not provider:
            sys.exit("No LiteLLM reranker provider selected.")
        server = resolve_server(load_litellm_config(), provider)
        return partial(rerank_litellm, server, model, user), candidates
    model_path = rerank_cfg.get("path")
    if not model_path:
        sys.exit("No reranker model path given.")
    if not os.path.exists(model_path):
        sys.exit(f"Reranker model path does not exist: {model_path}")
    return partial(rerank_local, model_path), candidates


def main():
    context_files = json.loads(sys.argv[1])
    question = (sys.argv[2] or "").strip()
    embed_cfg = json.loads(sys.argv[3])
    top_k = int(sys.argv[4])
    # Galaxy user id + instance URL, for request attribution on the proxy.
    # The literal "Anonymous" is rendered for anonymous sessions, in which
    # case no per-user id is sent (the request falls back to a per-instance
    # shared "anonymous" bucket).
    galaxy_user_id = sys.argv[5] if len(sys.argv) > 5 else ""
    if galaxy_user_id == "Anonymous":
        galaxy_user_id = ""
    galaxy_url = sys.argv[6] if len(sys.argv) > 6 else ""
    # Reranker settings; absent in command lines from before reranking existed.
    rerank_cfg = json.loads(sys.argv[7]) if len(sys.argv) > 7 and sys.argv[7].strip() else {"source": "none"}

    if not question:
        sys.exit("Question is empty.")
    if not context_files:
        sys.exit("No input files given.")
    if top_k <= 0:
        sys.exit("Top K must be a positive integer.")

    if not isinstance(embed_cfg, dict) or "source" not in embed_cfg:
        sys.exit("Invalid embedding configuration: expected a JSON object with a 'source' key.")
    if not isinstance(rerank_cfg, dict) or "source" not in rerank_cfg:
        sys.exit("Invalid reranker configuration: expected a JSON object with a 'source' key.")

    user = attribution_user(galaxy_user_id, galaxy_url)
    rerank, candidates = prepare_reranker(rerank_cfg, user)
    embed_model = build_embed_model(embed_cfg, user)

    docs: list[Document] = []

    valid_file_types = ["pdf", "json", "txt", "csv", "markdown"]

    for file_path, file_type in context_files:
        if file_type not in valid_file_types:
            sys.exit(f"Unsupported file type: {file_type} for file {file_path}")
        if file_type == "pdf":
            docs.extend(
                SimpleDirectoryReader(
                    input_files=[file_path],
                    file_extractor={".pdf": PDFReader()},
                ).load_data()
            )
        elif file_type == "json":
            docs.extend(JSONReader(levels_back=1).load_data(file_path))
        else:
            docs.extend(SimpleDirectoryReader(input_files=[file_path]).load_data())

    if not docs:
        sys.exit("No documents loaded.")

    index = VectorStoreIndex.from_documents(docs, embed_model=embed_model)
    # A reranker needs a wider pool to choose from; never fewer than top_k.
    fetch_k = max(candidates, top_k) if rerank else top_k
    retriever = index.as_retriever(similarity_top_k=fetch_k)
    nodes = retriever.retrieve(question)

    chunks = []
    for n in nodes:
        node = getattr(n, "node", n)
        text = node.get_content()
        if text:
            chunks.append(text.strip())

    if rerank and chunks:
        chunks = rerank(question, chunks, top_k)

    context_text = "\n\n---\n\n".join(chunks).strip()

    out = [
        "## Retrieved context\n",
        context_text if context_text else "(No context retrieved)",
    ]
    Path("rag_context.txt").write_text("\n".join(out), encoding="utf-8")


if __name__ == "__main__":
    main()
