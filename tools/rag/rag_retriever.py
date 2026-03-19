import json
import sys
import os
import torch

from pathlib import Path
from llama_index.core import Document, SimpleDirectoryReader, VectorStoreIndex
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.readers.file import PDFReader
from llama_index.readers.json import JSONReader


def main():
    context_files = json.loads(sys.argv[1])
    question = (sys.argv[2] or "").strip()
    embedding_model = sys.argv[3]
    top_k = int(sys.argv[4])

    if not question:
        sys.exit("Question is empty.")
    if not context_files:
        sys.exit("No input files given.")
    if not os.path.exists(embedding_model):
        sys.exit(f"Embedding model path does not exist: {embedding_model}")
    if top_k <= 0:
        sys.exit("Top K must be a positive integer.")

    device = "cuda" if torch.cuda.is_available() else "cpu"

    embed_model = HuggingFaceEmbedding(
        model_name=embedding_model, normalize=True, device=device
    )

    docs: list[Document] = []

    valid_file_types = ["pdf", "jsonl", "json", "txt", "csv", "md"]

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
        elif file_type == "jsonl":
            docs.extend(JSONReader(is_jsonl=True, levels_back=1).load_data(file_path))
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
