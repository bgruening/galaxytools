import hashlib
import json
import logging
import os
import sys
from functools import cache, partial
from pathlib import Path

import torch
import yaml
from llama_index.core import Document, SimpleDirectoryReader, VectorStoreIndex
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.embeddings.openai import OpenAIEmbedding
from llama_index.readers.file import PDFReader
from llama_index.readers.json import JSONReader
from openai import APIError, OpenAI
from sentence_transformers import CrossEncoder


# Limits for every request to the LiteLLM proxy (embeddings and reranking).
REQUEST_TIMEOUT = float(os.environ.get("LITELLM_REQUEST_TIMEOUT", "600"))
REQUEST_MAX_RETRIES = int(os.environ.get("LITELLM_REQUEST_MAX_RETRIES", "3"))
# Largest pool a reranker may re-score (the form's limit). Gains level off at
# about 20 candidates, and every candidate costs one more model pass.
MAX_RERANK_CANDIDATES = 100
# Cross-encoders judge the question and one chunk together; 512 tokens is the
# window the evaluated models (ms-marco-MiniLM, MedCPT, bge-reranker) were
# trained on.
RERANK_MAX_LENGTH = 512


# --- LiteLLM proxy config resolution -------------------------------------
# The LiteLLM proxy exposes OpenAI-compatible embeddings and rerank endpoints.
# The YAML config may be flat (global LITELLM_API_KEY / LITELLM_BASE_URL) or
# expose a ``servers`` mapping that keys provider names to per-server
# credentials; the ``provider`` argument selects which server to use.

@cache
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
    leaves for the provider. Anonymous sessions (Galaxy renders the literal
    "Anonymous", or no id) share one per-instance "anonymous" bucket. Mirrors
    the attribution in llm_hub.py so both tools map one Galaxy user to one
    proxy identity.
    """
    if not galaxy_user_id or galaxy_user_id == "Anonymous":
        galaxy_user_id = "anonymous"
    return hashlib.sha256(f"{galaxy_url}|{galaxy_user_id}".encode()).hexdigest()


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
            timeout=REQUEST_TIMEOUT,
            max_retries=REQUEST_MAX_RETRIES,
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

def is_cross_encoder(model_path: str) -> bool:
    """Whether ``model_path`` holds a single-score cross-encoder.

    That is a ``*ForSequenceClassification`` model with one output label.
    CrossEncoder loads other models too, e.g. an embedding model gets a new,
    randomly initialised scoring head and then ranks at random, so the check
    reads the model's config.json instead of trusting the load.
    """
    try:
        config = json.loads((Path(model_path) / "config.json").read_text(encoding="utf-8"))
        labels = config.get("id2label")
        num_labels = len(labels) if labels else config.get("num_labels", 2)
        architectures = config.get("architectures") or []
        return num_labels == 1 and any(a.endswith("ForSequenceClassification") for a in architectures)
    except (OSError, ValueError, AttributeError, TypeError):
        return False


def rerank_local(model_path: str, question: str, chunks: list, top_k: int) -> list:
    """Re-sort ``chunks`` with a local cross-encoder; return the best ``top_k``."""
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = CrossEncoder(model_path, max_length=RERANK_MAX_LENGTH, device=device)
    ranked = model.rank(question, chunks, top_k=top_k, batch_size=16, show_progress_bar=False)
    return [chunks[r["corpus_id"]] for r in ranked]


def rerank_litellm(client: OpenAI, model: str, user: str, question: str, chunks: list, top_k: int) -> list:
    """Re-sort ``chunks`` with a proxy-hosted reranker; return the best ``top_k``.

    Posts to the proxy's rerank endpoint (``{LITELLM_BASE_URL}/rerank``)
    through the OpenAI SDK client, so timeouts, retries (honouring the proxy's
    Retry-After) and redirects behave as for the embedding requests. The
    answer is checked before use: a proxy that returns too few results, or
    indices that do not point at a chunk, fails the job instead of silently
    returning wrong or short context.
    """
    top_n = min(top_k, len(chunks))
    try:
        response = client.post(
            "/rerank",
            cast_to=object,
            body={"model": model, "query": question, "documents": chunks, "top_n": top_n, "user": user},
        )
    except APIError as e:
        sys.exit(f"Reranker request to model '{model}' failed: {e}")
    try:
        ranked = sorted(response["results"], key=lambda r: -float(r["relevance_score"]))
        order = list(dict.fromkeys(r["index"] for r in ranked))  # best first, repeats dropped
    except (KeyError, TypeError, ValueError):
        sys.exit(f"Unexpected response from reranker model '{model}': {str(response)[:500]}")
    if len(order) < top_n or not all(type(i) is int and 0 <= i < len(chunks) for i in order):
        sys.exit(
            f"Unexpected response from reranker model '{model}': expected {top_n} "
            f"distinct chunk indices from 0 to {len(chunks) - 1}, got {order}"
        )
    return [chunks[i] for i in order[:top_n]]


def prepare_reranker(rerank_cfg: dict, user: str, top_k: int):
    """Validate the reranker settings before any document is embedded.

    Returns ``(rerank, fetch_k)``: a function ``rerank(question, chunks,
    top_k)``, or None when reranking is off, and how many chunks retrieval
    should return. Bad settings exit here, so a job fails in seconds rather
    than after embedding the whole corpus.
    """
    source = rerank_cfg["source"]
    if source == "none":
        return None, top_k
    if source not in ("litellm", "local"):
        sys.exit(f"Unknown reranker source: {source}")
    candidates = rerank_cfg.get("candidates")
    if type(candidates) is not int or not 1 <= candidates <= MAX_RERANK_CANDIDATES:
        sys.exit(f"Reranker candidates must be a whole number from 1 to {MAX_RERANK_CANDIDATES}.")
    # A reranker needs a wider pool to choose from; never fewer than top_k.
    fetch_k = max(candidates, top_k)
    if source == "litellm":
        model = rerank_cfg.get("model")
        provider = rerank_cfg.get("provider")
        if not model:
            sys.exit("No LiteLLM reranker model selected.")
        if not provider:
            sys.exit("No LiteLLM reranker provider selected.")
        server = resolve_server(load_litellm_config(), provider)
        client = OpenAI(
            api_key=server["LITELLM_API_KEY"],
            base_url=server["LITELLM_BASE_URL"],
            timeout=REQUEST_TIMEOUT,
            max_retries=REQUEST_MAX_RETRIES,
        )
        return partial(rerank_litellm, client, model, user), fetch_k
    model_path = rerank_cfg.get("path")
    if not model_path:
        sys.exit("No reranker model path given.")
    if not os.path.exists(model_path):
        sys.exit(f"Reranker model path does not exist: {model_path}")
    if not is_cross_encoder(model_path):
        sys.exit(
            f"Reranker model is not a single-score cross-encoder: {model_path} "
            "(expected a *ForSequenceClassification model with one output, "
            "e.g. cross-encoder/ms-marco-MiniLM-L6-v2)."
        )
    return partial(rerank_local, model_path), fetch_k


def main():
    # httpx logs every proxy request at INFO; keep the job's stderr for errors.
    logging.getLogger("httpx").setLevel(logging.WARNING)

    context_files = json.loads(sys.argv[1])
    question = (sys.argv[2] or "").strip()
    embed_cfg = json.loads(sys.argv[3])
    top_k = int(sys.argv[4])
    # Galaxy user id + instance URL, for request attribution on the proxy.
    galaxy_user_id = sys.argv[5]
    galaxy_url = sys.argv[6]
    rerank_cfg = json.loads(sys.argv[7])

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
    rerank, fetch_k = prepare_reranker(rerank_cfg, user, top_k)
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
