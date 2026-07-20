import json
import os
import sys
from pathlib import Path
from typing import Any

import torch
import yaml
from llama_index.core import Document, SimpleDirectoryReader, VectorStoreIndex
from llama_index.core.embeddings import BaseEmbedding
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.readers.file import PDFReader
from llama_index.readers.json import JSONReader


def _load_litellm_config() -> dict:
    """Read the LiteLLM YAML config and validate its presence.

    Mirrors the resolution logic in tools/llm_hub/llm_hub.py so that RAG and
    LLM Hub share the same server configuration.
    """
    config_file = os.environ.get("LITELLM_CONFIG_FILE")
    if not config_file:
        sys.exit("LiteLLM config file is not configured! Please set the LITELLM_CONFIG_FILE environment variable.")
    if not os.path.isfile(config_file):
        sys.exit(f"LiteLLM config file does not exist: {config_file}")
    with open(config_file, "r") as f:
        return yaml.safe_load(f)


def _resolve_server(config: dict, provider: str) -> dict:
    """Resolve the server config for a provider from the YAML config."""
    servers = config.get("servers", {})
    if servers:
        if provider not in servers:
            sys.exit(f"Provider '{provider}' not found in LiteLLM configuration.")
        source = servers[provider]
    else:
        source = config
    if not source.get("LITELLM_API_KEY"):
        sys.exit("LiteLLM API key is not configured! Please set LITELLM_API_KEY in the configuration.")
    if not source.get("LITELLM_BASE_URL"):
        sys.exit("LiteLLM base URL is not configured! Please set LITELLM_BASE_URL in the configuration.")
    return source


class LiteLLMEmbedding(BaseEmbedding):
    """Embedding model backed by an OpenAI-compatible LiteLLM `/v1/embeddings` endpoint.

    The LiteLLM proxy (the same one used by the LLM Hub) exposes embedding
    models such as BGE-M3, Qwen3-Embedding or nomic-embed-text. We call it
    directly through the ``openai`` SDK rather than downloading a local
    HuggingFace model, which keeps the corpus embeddings and the query
    embedding consistent with a single server-side model.
    """

    def __init__(
        self,
        model: str,
        api_key: str,
        base_url: str,
        embed_batch_size: int = 100,
        **kwargs: Any,
    ) -> None:
        from openai import OpenAI

        super().__init__(
            model_name=model,
            embed_batch_size=embed_batch_size,
            **kwargs,
        )
        self._model = model
        self._client = OpenAI(
            api_key=api_key,
            base_url=base_url,
            timeout=float(os.environ.get("LITELLM_REQUEST_TIMEOUT", "600")),
            max_retries=0,
        )

    def _embed(self, texts: list[str]) -> list[list[float]]:
        # The OpenAI SDK normalises input to a list of strings and returns one
        # embedding per input, preserving order.
        response = self._client.embeddings.create(model=self._model, input=texts)
        return [d.embedding for d in sorted(response.data, key=lambda x: x.index)]

    def _get_query_embedding(self, query: str) -> list[float]:
        return self._embed([query])[0]

    def _get_text_embedding(self, text: str) -> list[float]:
        return self._embed([text])[0]

    def _get_text_embeddings(self, texts: list[str]) -> list[list[float]]:
        if not texts:
            return []
        batch = self.embed_batch_size or len(texts)
        out: list[list[float]] = []
        for i in range(0, len(texts), batch):
            out.extend(self._embed(texts[i:i + batch]))
        return out

    async def _aget_query_embedding(self, query: str) -> list[float]:
        return self._get_query_embedding(query)

    async def _aget_text_embedding(self, text: str) -> list[float]:
        return self._get_text_embedding(text)

    async def _aget_text_embeddings(self, texts: list[str]) -> list[list[float]]:
        return self._get_text_embeddings(texts)


def main():
    context_files = json.loads(sys.argv[1])
    question = (sys.argv[2] or "").strip()
    embed_source = sys.argv[3]
    model_path = sys.argv[4]
    model_value = sys.argv[5]
    provider = sys.argv[6]
    top_k = int(sys.argv[7])

    if not question:
        sys.exit("Question is empty.")
    if not context_files:
        sys.exit("No input files given.")
    if top_k <= 0:
        sys.exit("Top K must be a positive integer.")

    embed_model: BaseEmbedding
    if embed_source == "litellm":
        config = _load_litellm_config()
        if not model_value:
            sys.exit("No LiteLLM embedding model selected.")
        if not provider:
            sys.exit("No LiteLLM provider selected.")
        server = _resolve_server(config, provider)
        embed_model = LiteLLMEmbedding(
            model=model_value,
            api_key=server["LITELLM_API_KEY"],
            base_url=server["LITELLM_BASE_URL"],
        )
    else:
        # Local HuggingFace model (preinstalled path or uploaded archive).
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
