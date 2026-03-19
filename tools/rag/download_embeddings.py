"""
Download and store embedding models in a local directory for Galaxy deployments.

Usage:
    export HF_MODEL_DIR=/path/to/huggingface_models
    python download_embeddings.py
"""

import os
import sys
from sentence_transformers import SentenceTransformer
from pathlib import Path

models = [
    "sentence-transformers/all-MiniLM-L6-v2",
    "BAAI/bge-small-en",
]

base_path = os.environ.get("HF_MODEL_DIR")
if not base_path:
    sys.exit("HF_MODEL_DIR environment variable is not set.")

base_path = Path(base_path)

for model in models:
    out_dir = base_path / model
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {model} -> {out_dir}")
    m = SentenceTransformer(model)
    m.save(str(out_dir))
