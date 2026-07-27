import json
import os
import sys
from pathlib import Path

import torch
import yaml
from llama_index.core import Document, SimpleDirectoryReader, VectorStoreIndex
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.embeddings.openai import OpenAIEmbedding
from llama_index.readers.file import PDFReader
from llama_index.readers.json import JSONReader


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


def main():
    context_files = json.loads(sys.argv[1])
    question = (sys.argv[2] or "").strip()
    embed_cfg = json.loads(sys.argv[3])
    top_k = int(sys.argv[4])

    if not question:
        sys.exit("Question is empty.")
    if not context_files:
        sys.exit("No input files given.")
    if top_k <= 0:
        sys.exit("Top K must be a positive integer.")

    if not isinstance(embed_cfg, dict) or "source" not in embed_cfg:
        sys.exit("Invalid embedding configuration: expected a JSON object with a 'source' key.")

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
        embed_model = OpenAIEmbedding(
            model_name=model,
            api_key=server["LITELLM_API_KEY"],
            api_base=server["LITELLM_BASE_URL"],
            timeout=float(os.environ.get("LITELLM_REQUEST_TIMEOUT", "600")),
            max_retries=int(os.environ.get("LITELLM_REQUEST_MAX_RETRIES", "3")),
            embed_batch_size=100,
        )
    else:
        # Local HuggingFace model (preinstalled path or uploaded archive).
        model_path = embed_cfg.get("path")
        if not model_path:
            sys.exit("No embedding model path given.")
        if not os.path.exists(model_path):
            sys.exit(f"Embedding model path does not exist: {model_path}")
        device = "cuda" if torch.cuda.is_available() else "cpu"
        embed_model = HuggingFaceEmbedding(
            model_name=model_path, normalize=True, device=device
        )

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
    retriever = index.as_retriever(similarity_top_k=top_k)
    nodes = retriever.retrieve(question)

    chunks = []
    for n in nodes:
        node = getattr(n, "node", n)
        text = node.get_content()
        if text:
            chunks.append(text.strip())

    context_text = "\n\n---\n\n".join(chunks).strip()

    out = [
        "## Retrieved context\n",
        context_text if context_text else "(No context retrieved)",
    ]
    Path("rag_context.txt").write_text("\n".join(out), encoding="utf-8")


if __name__ == "__main__":
    main()
