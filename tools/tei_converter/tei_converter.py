#!/usr/bin/env python3
"""Write the entities of a JSON Lines file as stand-off TEI annotations for Recogito Studio.

Input: JSON Lines, one object per paragraph: {"text": ..., "entities": [{"start", "end", ...}]}.
Every entity that has an id in --id-field (default "reranker_results") becomes one annotation that links its
characters to <base-url><id>. Entities without an id are skipped. The result is one TEI file for the whole input.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from itertools import accumulate
from pathlib import Path

from lxml import etree

NS_TEI = "http://www.tei-c.org/ns/1.0"
NS_XML = "http://www.w3.org/XML/1998/namespace"
NAME_RE = re.compile(r"[\w.-]+")
URL_RE = re.compile(r"https?://\S+")
XML_BAD = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f￾￿]")  # characters XML 1.0 cannot hold


def log(message: str) -> None:
    print(f"[tei_converter] {message}", file=sys.stderr, flush=True)


def read_rows(path: Path, field: str) -> list[dict]:
    raw = path.read_text(encoding="utf-8-sig")
    try:
        data = json.loads(raw)
        rows = data if isinstance(data, list) else [data]
    except json.JSONDecodeError:  # JSON Lines
        rows = [json.loads(line) for line in raw.splitlines() if line.strip()]
    for number, row in enumerate(rows, 1):
        if not isinstance(row, dict) or not isinstance(row.get("text"), str) or not isinstance(row.get(field, []), list):
            raise SystemExit(f"Input item {number} is not an object with a string 'text' and a {field!r} list.")
    if not any(field in row for row in rows):
        raise SystemExit(f"No item has a {field!r} list. Set the field that holds the entities, for example "
                         "'entities' (spaCy NLP Annotators) as spaCy NLP Annotators writes it.")
    return rows


def clean(text: str) -> tuple[str, list[int] | None]:
    """The text without characters XML cannot hold, and for every old position its new one (None if unchanged)."""
    if not XML_BAD.search(text):
        return text, None
    removed = [1 if XML_BAD.fullmatch(char) else 0 for char in text]
    shift = [0, *accumulate(removed)]  # characters removed before each position (and before the end)
    return XML_BAD.sub("", text), [i - shift[i] for i in range(len(text) + 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--entities-field", default="entities", help="item field with the list of entities")
    parser.add_argument("--id-field", default="reranker_results", help="entity field with the id of the chosen entry")
    parser.add_argument("--base-url", default="", help="put in front of the id, e.g. https://topostext.org/place/")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    for option, value in (("--id-field", args.id_field), ("--entities-field", args.entities_field)):
        if not NAME_RE.fullmatch(value):
            raise SystemExit(f"Invalid value for {option}: {value!r}.")
    if args.base_url and not URL_RE.fullmatch(args.base_url):
        raise SystemExit(f"The base URL must start with http:// or https://, not {args.base_url!r}.")
    rows = read_rows(args.input, args.entities_field)

    tei = etree.Element("TEI", nsmap={None: NS_TEI}, version="3.3.0")
    annotations = etree.SubElement(etree.SubElement(tei, "standOff", type="recogito_studio_annotations"), "listAnnotation")
    body = etree.SubElement(etree.SubElement(tei, "text"), "body")

    written = skipped = 0
    for counter, row in enumerate(rows, 1):
        text, new_position = clean(row["text"])
        etree.SubElement(body, "p").text = text
        path = f"/TEI[1]/text[1]/body[1]/p[{counter}]"
        for entity in row.get(args.entities_field, []):
            chosen = entity.get(args.id_field)
            if isinstance(chosen, list):  # e.g. the ranked ids of the search step: take the best
                chosen = chosen[0] if chosen else None
            if chosen is None or chosen == "":
                skipped += 1
                continue
            start, end = entity["start"], entity["end"]
            if not (isinstance(start, int) and isinstance(end, int) and 0 <= start <= end <= len(row["text"])):
                raise SystemExit(f"Paragraph {counter} has an entity with positions {start}:{end} outside its text.")
            if new_position:
                start, end = new_position[start], new_position[end]
            note = etree.SubElement(annotations, "annotation", ana=f"{args.base_url}{chosen}",
                                    target=f"{path}::{start} {path}::{end}")
            note.set(f"{{{NS_XML}}}id", f"UID-FAKE-{written}")  # for recogito studio
            written += 1

    args.output.write_bytes(etree.tostring(tei, xml_declaration=True, encoding="utf-8", pretty_print=True))
    log(f"{len(rows)} paragraph(s), {written} annotation(s) written, {skipped} entity(ies) without an id skipped.")


if __name__ == "__main__":
    main()
