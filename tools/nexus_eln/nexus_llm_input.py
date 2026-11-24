#!/usr/bin/env python
"""Write the LLM input for one raw file: the files one ``llm_hub`` job reads.

Inputs are the raw file's context bundle (``zenodo_nomad_download``), its header
(``spm_raw_header``), and the rules and instructions of the pynxtools plugin. Two stages:

- ``detect``: the input of the job that finds the technique and scan mode;
- ``eln``: the input of the job that fills the rule ELN (``nexus_eln_builder``).

Written to ``--out-dir``:

- ``eln_instructions.txt``: the stage's sections of the instructions file;
- ``base_eln.txt`` (``eln`` only): the rule ELN, with the fields the LLM may fill set to null;
- ``<raw file>_header.json`` (the header's metadata), or ``<raw file>_header.txt`` when the
  header is the file's leading text;
- ``record.txt``: the record's title and description;
- one text file per context file with text (a PDF's ``pdf_text``), named after the file.

Without a definition, the ``eln`` stage writes only ``base_eln.txt``; the check step then keeps
the rule ELN. The report lists the files, the files left out and the ``Example:`` values.
"""

import argparse
import copy
import html
import json
import logging
import re
import sys
from pathlib import Path, PurePosixPath

import yaml

logger = logging.getLogger("nexus_llm_input")

INSTRUCTIONS_NAME, BASE_ELN_NAME, RECORD_NAME = "eln_instructions.txt", "base_eln.txt", "record.txt"
BASE_ELN_NOTE = "# Fill only the null fields; keep every other value.\n"
# Sections of the instructions file; the ELN job also gets the section of its definition.
DETECTION, COMMON, OUTPUT = "Detection", "Common", "Output"


def safe_name(name):
    """A plain file name; same rule as zenodo_nomad_download."""
    return re.sub(r"[^\w.\-]", "_", name)


def unique_name(name, taken):
    """`name`, or `name` with `_2`, `_3`, ... before its extension if it is taken."""
    path, number = PurePosixPath(name), 1
    while name in taken:
        number += 1
        name = f"{path.stem}_{number}{path.suffix}"
    taken.add(name)
    return name


def has_path(eln, dotted):
    node = eln
    for key in dotted.split("."):
        if not isinstance(node, dict) or key not in node:
            return False
        node = node[key]
    return node not in (None, "", [])


def set_path(eln, dotted, value):
    node = eln
    *parents, last = dotted.split(".")
    for key in parents:
        node = node.setdefault(key, {})
    node[last] = value


def base_eln(eln, rules):
    """The rule ELN with every empty `fill` field and every `replace` field set to null."""
    base = copy.deepcopy(eln)
    for field in rules["fill"]:
        if not has_path(base, field):
            set_path(base, field, None)
    for field in rules["replace"]:
        set_path(base, field, None)
    return base


def definition_of(eln, rules):
    """The ELN's definition, else the one of its technique; None if neither is known."""
    definitions = rules["definitions"]
    if eln.get("definition") in definitions:
        return eln["definition"]
    by_technique = {value["technique"]: name for name, value in definitions.items()}
    return by_technique.get(eln.get("experiment_technique"))


def plain_text(description_html):
    """Text of Zenodo's HTML description."""
    text = re.sub(r"<br\s*/?>|</p>", "\n", description_html or "", flags=re.IGNORECASE)
    text = html.unescape(re.sub(r"<[^>]+>", "", text))
    return re.sub(r"\n\s*\n+", "\n\n", text).strip()


def read_sections(path):
    """{name: text} of the `# <name>` sections; the text before the first one is a note, not sent.

    Lines inside fenced code blocks never start a section.
    """
    sections, name, fenced = {}, None, False
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("```"):
            fenced = not fenced
        if not fenced and line.startswith("# "):
            name = line[2:].strip()
            sections[name] = []
        elif name:
            sections[name].append(line)
    return {key: "\n".join(lines).strip() + "\n" for key, lines in sections.items()}


def checked_text(sections, names, needed, what):
    """The named sections joined; stops if one is missing or misses a value of the rules."""
    absent = [name for name in names if name not in sections]
    if absent:
        sys.exit(f"No section {', '.join(absent)} in the instructions file")
    text = "\n".join(sections[name] for name in names)
    missing = [value for value in needed if value not in text]
    if missing:
        sys.exit(f"Instructions for {what} do not match the rules; missing: {missing}")
    return text


def eln_instructions(definition, rules, path):
    """Instructions of the ELN job: Common, the definition's section, Output."""
    needed = rules["definitions"][definition]["scan_mode"] + rules["fill"] + rules["replace"]
    return checked_text(read_sections(path), [COMMON, definition, OUTPUT], needed, definition)


def detection_instructions(rules, path):
    """Instructions of the detection job; they name every technique and scan mode of the rules."""
    needed = [value for item in rules["definitions"].values()
              for value in [item["technique"], *item["scan_mode"]]]
    return checked_text(read_sections(path), [DETECTION], needed, "the detection")


def example_values(instructions):
    """Quoted values on `Example:` lines; the check step drops an ELN value that copies one."""
    lines = re.findall(r"^\s*Example:(.*)$", instructions, flags=re.MULTILINE)
    return [value for line in lines for value in re.findall(r'"([^"]+)"', line)]


def context_texts(bundle):
    """(path, text) per context file with text, and the files left out with a reason."""
    texts, left_out = [], []
    for item in bundle.get("context", []):
        text = item.get("pdf_text") if "pdf_text" in item else item.get("text")
        if text and text.strip():
            texts.append((item["path"], text))
        else:
            reason = item.get("pdf_text_error") or "no_text"
            left_out.append({"path": item["path"], "reason": reason})
    return texts, left_out


def write_input(bundle, eln, header, rules, instructions_path, out_dir):
    """Write the LLM input files to `out_dir`; return the report. Without `eln`: detection stage."""
    raw_path = bundle["raw"]["path"]
    definition = definition_of(eln, rules) if eln is not None else None
    files = []

    def write(name, text, kind, source=None):
        (out_dir / name).write_text(text, encoding="utf-8")
        files.append({"name": name, "kind": kind, "source": source, "chars": len(text)})

    report = {
        "schema_version": 1,
        "sha256": bundle["raw"]["sha256"],
        "raw_path": raw_path,
        "record_id": bundle["record"]["record_id"],
        "stage": "detect" if eln is None else "eln",
        "definition": definition,
        "status": "ok",
        "instructions": None,
        "files": files,
        "left_out": [],
        "example_values": [],
    }
    if eln is None:
        sections = [DETECTION]
        instructions = detection_instructions(rules, instructions_path)
    else:
        yaml_text = yaml.safe_dump(base_eln(eln, rules), allow_unicode=True, sort_keys=False)
        write(BASE_ELN_NAME, BASE_ELN_NOTE + yaml_text, "base_eln")
        if not definition:
            logger.warning("No definition for %s; only %s is written", raw_path, BASE_ELN_NAME)
            report["status"] = "no_definition"
            return report
        sections = [COMMON, definition, OUTPUT]
        instructions = eln_instructions(definition, rules, instructions_path)
        report["example_values"] = example_values(instructions)
    write(INSTRUCTIONS_NAME, instructions, "instructions", " + ".join(sections))
    report["instructions"] = sections

    sources, report["left_out"] = source_files(bundle, header)
    for source in sources:
        write(source["name"], source["text"], source["kind"], source["path"])
    return report


def header_source(header_path):
    """(suffix, text) of the header file of spm_raw_header; None if it holds no header.

    Metadata is given as the JSON file itself; a leading-text header as its plain text.
    """
    text = Path(header_path).read_text(encoding="utf-8")
    header = json.loads(text)
    if "text" in header:
        return (".txt", header["text"]) if header["text"].strip() else None
    return ".json", text


def record_text(record):
    """The record's title and description as text; empty if it has neither."""
    title, description = record.get("title") or "", plain_text(record.get("description_html"))
    parts = [f"Title: {title}" if title else "", description]
    return "\n\n".join(part for part in parts if part) + "\n" if title or description else ""


def source_files(bundle, header_path):
    """The source files the LLM gets, as dicts (name, kind, path, text), and the files left out.

    Also used by nexus_eln_llm_merge to check the quotes against the same texts.
    """
    taken, sources = {INSTRUCTIONS_NAME, BASE_ELN_NAME, RECORD_NAME}, []
    raw_path = bundle["raw"]["path"]
    header = header_source(header_path)
    if header:
        suffix, text = header
        name = unique_name(safe_name(f"{PurePosixPath(raw_path).stem}_header{suffix}"), taken)
        sources.append({"name": name, "kind": "header", "path": raw_path, "text": text})
    record = record_text(bundle["record"])
    if record:
        sources.append({"name": RECORD_NAME, "kind": "record", "path": None, "text": record})
    texts, left_out = context_texts(bundle)
    for path, text in texts:
        name = unique_name(safe_name(PurePosixPath(path).name), taken)
        sources.append({"name": name, "kind": "context", "path": path, "text": text})
    return sources, left_out


def main():
    parser = argparse.ArgumentParser(description="Write the LLM input for one raw file.")
    parser.add_argument("bundle", help="Context bundle JSON of one raw file.")
    parser.add_argument("--eln", help="ELN YAML of the rule builder; without it: the detection stage.")
    parser.add_argument("--header", required=True, help="Header JSON of the raw file (spm_raw_header).")
    parser.add_argument("--rules", required=True, help="Rules YAML of the pynxtools plugin.")
    parser.add_argument("--instructions", required=True, help="Instructions file of the pynxtools plugin.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    bundle = json.loads(Path(args.bundle).read_text(encoding="utf-8"))
    eln = None
    if args.eln:
        eln = yaml.safe_load(Path(args.eln).read_text(encoding="utf-8")) or {}
    rules = yaml.safe_load(Path(args.rules).read_text(encoding="utf-8"))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report = write_input(bundle, eln, args.header, rules, args.instructions, out_dir)
    report["total_chars"] = sum(item["chars"] for item in report["files"])
    Path(args.report).write_text(json.dumps(report, indent=2), encoding="utf-8")
    logger.info("LLM input for %s: %d files, %d left out",
                report["raw_path"], len(report["files"]), len(report["left_out"]))


if __name__ == "__main__":
    main()
