#!/usr/bin/env python3
"""Rerank the candidate knowledge-base entries of every entity with a sentence-transformers CrossEncoder.

Input and output are JSON Lines, one object per paragraph: {"text": ..., "entities": [{"start", "end", "label", ...}]}.
For every entity, its mention (the characters from "start" to "end" of its text) is scored against each candidate listed in
--candidates-field (default "index_results", the ids from the FAISS search). The text of a candidate is
"<title>: <description>", read from the knowledge-base metadata. The entity gets

  reranker_results  the id of the best candidate (null if there were no candidates)
  reranker_scores   [{"id": ..., "score": ...}, ...] for all candidates, best first
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
    # Galaxy caps the job at GALAXY_SLOTS CPUs; set this before torch starts its thread pools.
    os.environ["OMP_NUM_THREADS"] = os.environ["MKL_NUM_THREADS"] = SLOTS

# The xet downloader of huggingface_hub panics in containers without system CA certificates; plain HTTPS works.
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")

MODEL_RE = re.compile(r"[A-Za-z0-9][\w.-]*(/[\w.-]+)?")
NAME_RE = re.compile(r"[\w.-]+")


def log(message: str) -> None:
    print(f"[cross_encoder_rerank] {message}", file=sys.stderr, flush=True)


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


def mention(row: dict, entity: dict):
    """The text of the mention: the characters from start to end of the text. None if the positions are no good."""
    start, end, text = entity.get("start"), entity.get("end"), row.get("text")
    ok = isinstance(text, str) and isinstance(start, int) and isinstance(end, int) and 0 <= start <= end <= len(text)
    return text[start:end] if ok and text[start:end].strip() else None


def read_documents(path: Path, id_field: str, title_field: str, text_field: str) -> dict:
    """The text of every knowledge-base entry, by id. Without an id field the id is the line number."""
    documents = {}
    for number, entry in enumerate(read_jsonl(path)):
        if not isinstance(entry, dict):
            raise SystemExit(f"Metadata row {number + 1} is not a JSON object.")
        if id_field and id_field not in entry:
            raise SystemExit(f"Metadata row {number + 1} has no field {id_field!r}. Set the ID field, or leave it "
                             "empty to use the line number.")
        parts = [str(entry[field]).strip() for field in (title_field, text_field) if field and entry.get(field)]
        if not parts:
            raise SystemExit(f"Metadata row {number + 1} has neither {title_field!r} nor {text_field!r}.")
        documents[entry[id_field] if id_field else number] = ": ".join(parts)  # "title: description", as the notebook
    return documents


def load_model(name: str, dtype: str):
    import torch
    from sentence_transformers import CrossEncoder

    log(f"Loading reranker {name}.")
    kwargs = {"device": "cpu"}
    # float32 was about 4x faster than the bfloat16 most checkpoints are saved in, on the CPU it was measured on.
    # "auto" keeps the saved precision, which needs half the memory.
    if dtype == "float32":
        kwargs["model_kwargs"] = {"dtype": torch.float32}
    try:
        return CrossEncoder(name, **kwargs)  # trust_remote_code stays off: no code from the Hub is run
    except Exception as exc:
        hint = " This model needs custom code from its repository, which this tool does not run." \
            if "trust_remote_code" in str(exc) else ""
        raise SystemExit(f"Could not load {name} as a CrossEncoder: {exc}{hint}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--entities-field", default="entities", help="item field with the list of entities")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--id-field", default="", help="metadata field with the id; empty = line number")
    parser.add_argument("--title-field", default="title")
    parser.add_argument("--text-field", default="description")
    parser.add_argument("--model", required=True, help="CrossEncoder model on Hugging Face")
    parser.add_argument("--candidates-field", default="index_results", help="entity field with the candidate ids")
    parser.add_argument("--dtype", default="float32", choices=["float32", "auto"])
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    # These values end up in a download request or a dict key, so they are checked here as well as in the form.
    if not MODEL_RE.fullmatch(args.model):
        raise SystemExit(f"Invalid model name {args.model!r}.")
    for option, value in (("--id-field", args.id_field), ("--title-field", args.title_field),
                          ("--text-field", args.text_field),
                          ("--candidates-field", args.candidates_field), ("--entities-field", args.entities_field)):
        if value and not NAME_RE.fullmatch(value):
            raise SystemExit(f"Invalid value for {option}: {value!r}.")
    if not (args.title_field or args.text_field):
        raise SystemExit("Give a title field, a text field or both.")

    rows = read_rows(args.input, args.entities_field)
    documents = read_documents(args.metadata, args.id_field, args.title_field, args.text_field)

    work, bad_query, bad_ids = [], 0, set()  # work: (entity, mention, candidate ids)
    for row in rows:
        for entity in row.get(args.entities_field, []):
            candidates = entity.get(args.candidates_field) or []
            if not candidates:
                entity["reranker_results"], entity["reranker_scores"] = None, []
                continue
            query = mention(row, entity)
            if query is None:
                bad_query += 1
            bad_ids.update(c for c in candidates if c not in documents)
            work.append((entity, query, candidates))
    if bad_query:
        raise SystemExit(f"{bad_query} entity(ies) have empty or out-of-range 'start' and 'end' positions, so there is "
                         "no mention to match.")
    if bad_ids:
        raise SystemExit(f"{len(bad_ids)} candidate id(s) are not in the metadata, for example "
                         f"{sorted(map(str, bad_ids))[:3]}. Use the metadata of the knowledge base that was searched.")
    pairs = sum(len(candidates) for _, _, candidates in work)
    log(f"{len(rows)} text(s), {len(work)} entity(ies), {pairs} pair(s) to score.")

    if work:
        model = load_model(args.model, args.dtype)
        done = 0
        for entity, query, candidates in work:
            scores = model.predict([(query, documents[c]) for c in candidates], batch_size=args.batch_size,
                                   show_progress_bar=False, convert_to_numpy=True)
            ranked = sorted(zip(candidates, (float(s) for s in scores)), key=lambda pair: -pair[1])
            entity["reranker_results"] = ranked[0][0]
            entity["reranker_scores"] = [{"id": c, "score": s} for c, s in ranked]
            done += 1
            if done % 25 == 0:
                log(f"{done}/{len(work)} entities done.")

    args.output.write_text("".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows), encoding="utf-8")
    log("Done.")


if __name__ == "__main__":
    main()
