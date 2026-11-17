#!/usr/bin/env python3
"""Run one or more spaCy NER pipelines from the Hugging Face Hub over paragraphs of text; write JSON Lines.

Every output line is {"text": ..., "spans": [{"start", "end", "label"}, ...]}, the format of the T4.5.2 notebook.
The models run one after another (each is freed before the next is loaded) and all their spans end up in the
same list. Spans and other keys already on an input row are kept, so runs can also be chained.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import re
import sys
from pathlib import Path

SLOTS = os.environ.get("GALAXY_SLOTS")
if SLOTS:
    # Galaxy caps the job at GALAXY_SLOTS CPUs; set this before torch/numpy start their thread pools.
    os.environ["OMP_NUM_THREADS"] = os.environ["MKL_NUM_THREADS"] = SLOTS

REPO_RE = re.compile(r"[A-Za-z0-9][\w.-]*/[\w.-]+")
REVISION_RE = re.compile(r"[\w.-]+")
NER_COMPONENTS = ("transformer", "tok2vec", "ner")


def log(message: str) -> None:
    print(f"[spacy_hf_ner] {message}", file=sys.stderr, flush=True)


def read_txt(path: Path, split: str) -> list[str]:
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    blocks = re.split(r"\n\s*\n", text) if split == "blank" else text.splitlines()
    paragraphs = (" ".join(line.strip() for line in block.splitlines() if line.strip()) for block in blocks)
    return [paragraph for paragraph in paragraphs if paragraph]


def read_tei(path: Path) -> list[str]:
    from lxml import etree
    parser = etree.XMLParser(resolve_entities=False, no_network=True, huge_tree=False)
    root = etree.parse(str(path), parser).getroot()
    paragraphs = ("".join(p.itertext()).strip() for p in root.iter("{*}p"))
    return [paragraph for paragraph in paragraphs if paragraph]


def read_json(path: Path) -> list[dict]:
    raw = path.read_text(encoding="utf-8-sig")
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        data = [json.loads(line) for line in raw.splitlines() if line.strip()]
    rows = data if isinstance(data, list) else [data]
    for number, row in enumerate(rows, 1):
        if not isinstance(row, dict) or not isinstance(row.get("text"), str):
            raise SystemExit(f"Input item {number} is not an object with a string 'text' field.")
    return rows


def padded(values: list[str], count: int, option: str) -> list[str | None]:
    """The values of a repeated option, one per model; an option left out entirely means 'default' for all."""
    if not values:
        return [None] * count
    if len(values) != count:
        raise SystemExit(f"{option} must be given once per --model ({count} times), or not at all.")
    return [value or None for value in values]


def requirement_notes(meta: dict) -> str:
    """Say which packages the model asks for that differ from the installed ones."""
    try:
        from importlib.metadata import PackageNotFoundError, version
        from packaging.requirements import Requirement
    except ImportError:
        return ""
    notes = []
    for text in meta.get("requirements", []):
        requirement = Requirement(text)
        try:
            installed = version(requirement.name)
        except PackageNotFoundError:
            continue
        if not requirement.specifier.contains(installed, prereleases=True):
            notes.append(f"the model asks for {text}, this tool has {requirement.name} {installed}")
    return "; ".join(notes)


def load_pipeline(repo: str, revision: str | None):
    import spacy
    from huggingface_hub import snapshot_download

    log(f"Downloading {repo}" + (f" at {revision}" if revision else "") + ".")
    # The repo holds the unpacked pipeline next to a wheel of it; the wheel would double the download.
    folder = Path(snapshot_download(repo, revision=revision, ignore_patterns=["*.whl", "*.md", ".gitattributes"]))
    if not (folder / "meta.json").is_file() or not (folder / "config.cfg").is_file():
        raise SystemExit(f"{repo} is not a spaCy pipeline repository (no meta.json and config.cfg at its root).")
    notes = requirement_notes(json.loads((folder / "meta.json").read_text(encoding="utf-8")))
    if notes:
        log(f"Warning: {notes}.")
    try:
        nlp = spacy.load(folder)
    except Exception as exc:  # the error is usually a package version problem, so show both
        raise SystemExit(f"Could not load {repo}." + (f"\nLikely cause: {notes}." if notes else "")
                         + f"\nOriginal error: {exc}")
    if "ner" not in nlp.pipe_names:
        raise SystemExit(f"{repo} has no 'ner' component (it has: {', '.join(nlp.pipe_names)}).")
    # Only what NER needs
    nlp.select_pipes(enable=[name for name in nlp.pipe_names if name in NER_COMPONENTS])
    return nlp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--format", required=True, choices=["txt", "tei", "json"])
    parser.add_argument("--split", default="blank", choices=["blank", "line"])
    # --model, --revision and --label are repeated once per model, in the same order.
    parser.add_argument("--model", action="append", required=True, help="Hugging Face repository, owner/name")
    parser.add_argument("--revision", action="append", default=[], help="branch, tag or commit hash ('' = latest)")
    parser.add_argument("--label", action="append", default=[], help="label to use instead of the model's own")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    # These values end up in a download request, so they are checked here as well as in the Galaxy form.
    revisions = padded(args.revision, len(args.model), "--revision")
    labels = padded(args.label, len(args.model), "--label")
    for repo, revision in zip(args.model, revisions):
        if not REPO_RE.fullmatch(repo):
            raise SystemExit(f"Invalid model name {repo!r}; expected owner/name.")
        if revision and not REVISION_RE.fullmatch(revision):
            raise SystemExit(f"Invalid revision {revision!r}.")

    if args.format == "json":
        rows = read_json(args.input)
    else:
        paragraphs = read_tei(args.input) if args.format == "tei" else read_txt(args.input, args.split)
        rows = [{"text": paragraph} for paragraph in paragraphs]
    if not rows:
        raise SystemExit("No text was found in the input.")
    log(f"{len(rows)} paragraph(s) to analyse.")

    # Load every model first, so a wrong name or an incompatible model fails now and not after a long first pass.
    pipelines = [(repo, label, load_pipeline(repo, revision)) for repo, revision, label in zip(args.model, revisions, labels)]
    while pipelines:
        repo, label, nlp = pipelines.pop(0)  # popped, so its memory is freed after its pass over the whole text
        texts = (row["text"] for row in rows)
        for number, (row, doc) in enumerate(zip(rows, nlp.pipe(texts, batch_size=args.batch_size)), 1):
            row["spans"] = list(row.get("spans") or []) + [
                {"start": ent.start_char, "end": ent.end_char, "label": label or ent.label_} for ent in doc.ents
            ]
            if number % 50 == 0:
                log(f"{repo}: {number}/{len(rows)} paragraphs done.")
        del nlp  # only one transformer model in memory at a time
        gc.collect()

    args.output.write_text("".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows), encoding="utf-8")
    log(f"Done: {len(rows)} paragraph(s), {sum(len(row['spans']) for row in rows)} span(s).")


if __name__ == "__main__":
    main()
