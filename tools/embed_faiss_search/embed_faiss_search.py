#!/usr/bin/env python3
"""Embed the text written for each entity and look it up in a FAISS index of a knowledge base.

Input and output are JSON Lines, one object per paragraph: {"text": ..., "entities": [{"start", "end", "label", ...}]}.
For every entity the mention itself (the characters from "start" to "end" of its text) is embedded. The --top-k
nearest knowledge-base entries are added to the entity as "index_results", a list of their ids, best first.

The FAISS index and the metadata file go together: row i of the metadata (JSON Lines) describes vector i of the index.
The embedding model must be the one that built the index; the index does not say which one it was, so the vector
size is checked as a guard.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

SLOTS = os.environ.get("GALAXY_SLOTS")
if SLOTS:
    # Galaxy caps the job at GALAXY_SLOTS CPUs; set this before torch/faiss start their thread pools.
    os.environ["OMP_NUM_THREADS"] = os.environ["MKL_NUM_THREADS"] = SLOTS

# The xet downloader of huggingface_hub panics in containers without system CA certificates; plain HTTPS works.
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")

MODEL_RE = re.compile(r"[A-Za-z0-9][\w.-]*(/[\w.-]+)?")
NAME_RE = re.compile(r"[\w.-]+")


def log(message: str) -> None:
    print(f"[embed_faiss_search] {message}", file=sys.stderr, flush=True)


def read_jsonl(path: Path) -> list:
    """A JSON list, or one JSON value per line."""
    raw = path.read_text(encoding="utf-8-sig")
    try:
        data = json.loads(raw)
        return data if isinstance(data, list) else [data]
    except json.JSONDecodeError:
        return [json.loads(line) for line in raw.splitlines() if line.strip()]


def read_rows(path: Path, field: str) -> list[dict]:
    rows = read_jsonl(path)
    for number, row in enumerate(rows, 1):
        if not isinstance(row, dict) or not isinstance(row.get(field, []), list):
            raise SystemExit(f"Input item {number} is not an object with a {field!r} list.")
    if not any(field in row for row in rows):
        raise SystemExit(f"No item has a {field!r} list. Set the field that holds the entities, for example "
                         "'entities', as spaCy NLP Annotators writes it.")
    return rows


def read_ids(path: Path, id_field: str) -> list:
    """The id of every knowledge-base entry, in the order of the index."""
    entries = read_jsonl(path)
    if not id_field:
        return list(range(len(entries)))
    ids = []
    for number, entry in enumerate(entries, 1):
        if not isinstance(entry, dict) or id_field not in entry:
            raise SystemExit(f"Metadata row {number} has no field {id_field!r}. Set the ID field, or leave it empty "
                             "to use the row number.")
        ids.append(entry[id_field])
    return ids


def mention(row: dict, entity: dict):
    """The text of the mention: the characters from start to end of the text. None if the positions are no good."""
    start, end, text = entity.get("start"), entity.get("end"), row.get("text")
    ok = isinstance(text, str) and isinstance(start, int) and isinstance(end, int) and 0 <= start <= end <= len(text)
    return text[start:end] if ok and text[start:end].strip() else None


def load_model(name: str, prompt: str, dtype: str):
    import torch
    from sentence_transformers import SentenceTransformer

    log(f"Loading embedding model {name}.")
    # "auto" keeps the precision the model was saved in, as the notebook does and as the index was built.
    # float32 is a little more precise, but then neighbours at the bottom of the top-10 can swap places.
    kwargs = dict(device="cpu")
    if dtype == "float32":
        kwargs["model_kwargs"] = {"dtype": torch.float32}
    try:
        model = SentenceTransformer(name, **kwargs)
    except Exception:
        # The notebook loads EmbeddingGemma 2 without its vision and audio parts; retry that way before giving up.
        model = SentenceTransformer(name, config_kwargs={"vision_config": None, "audio_config": None}, **kwargs)
    if prompt and prompt not in (model.prompts or {}):
        raise SystemExit(f"{name} has no prompt named {prompt!r}. It has: {', '.join(model.prompts or {}) or 'none'}. "
                         "Leave the query prompt empty to use none.")
    return model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--entities-field", default="entities", help="item field with the list of entities")
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--id-field", default="", help="metadata field with the id; empty = row number")
    parser.add_argument("--model", required=True, help="embedding model on Hugging Face, the one that built the index")
    parser.add_argument("--prompt", default="", help="name of the model's query prompt, empty = none")
    parser.add_argument("--top-k", type=int, default=30)
    parser.add_argument("--dtype", default="auto", choices=["auto", "float32"], help="precision of the embedding model")
    parser.add_argument("--normalize", action="store_true", help="normalize the query vectors, as when the index was built")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    # These values end up in a download request or a dict key, so they are checked here as well as in the form.
    if not MODEL_RE.fullmatch(args.model):
        raise SystemExit(f"Invalid model name {args.model!r}.")
    for option, value in (("--prompt", args.prompt), ("--id-field", args.id_field),
                          ("--entities-field", args.entities_field)):
        if value and not NAME_RE.fullmatch(value):
            raise SystemExit(f"Invalid value for {option}: {value!r}.")

    rows = read_rows(args.input, args.entities_field)
    queries, missing = [], 0
    for row in rows:
        for entity in row.get(args.entities_field, []):
            text = mention(row, entity)
            if text:
                queries.append((entity, text))
            else:
                missing += 1
    if missing:
        raise SystemExit(f"{missing} entity(ies) have empty or out-of-range 'start' and 'end' positions, so there is "
                         "no mention to look up.")
    log(f"{len(rows)} text(s), {len(queries)} entity(ies) to look up.")

    import numpy as np
    import torch  # noqa: F401  torch first: imported after faiss, their OpenMP runtimes clash and crash on macOS
    import faiss

    if SLOTS:
        faiss.omp_set_num_threads(int(SLOTS))
    index = faiss.read_index(str(args.index))
    kb_ids = read_ids(args.metadata, args.id_field)
    if index.ntotal != len(kb_ids):
        raise SystemExit(f"The index has {index.ntotal} vectors but the metadata has {len(kb_ids)} rows. They must "
                         "come from the same knowledge base, in the same order.")

    if queries:
        model = load_model(args.model, args.prompt, args.dtype)
        size = model.get_embedding_dimension()
        if size != index.d:
            raise SystemExit(f"{args.model} makes vectors of size {size} but the index holds vectors of size "
                             f"{index.d}. Use the model that built the index.")
        for start in range(0, len(queries), 256):
            chunk = queries[start:start + 256]
            vectors = model.encode([text for _, text in chunk], prompt_name=args.prompt or None,
                                   normalize_embeddings=args.normalize, batch_size=args.batch_size,
                                   convert_to_numpy=True, show_progress_bar=False)
            _, hits = index.search(np.ascontiguousarray(vectors, dtype="float32"), args.top_k)
            for (entity, _), found in zip(chunk, hits):
                entity["index_results"] = [kb_ids[i] for i in found if i >= 0]  # -1 means fewer entries than top-k
            log(f"{min(start + 256, len(queries))}/{len(queries)} entities done.")

    args.output.write_text("".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows), encoding="utf-8")
    log("Done.")


if __name__ == "__main__":
    main()
