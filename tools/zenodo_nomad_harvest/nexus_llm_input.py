#!/usr/bin/env python
"""Write the LLM input for one raw SPM file: the files one ``llm_hub`` job reads.

Inputs are the raw file's context bundle (``zenodo_download_unzip``), its rule ELN
(``nexus_eln_builder``) and, optionally, the raw file. Written to ``--out-dir``:

- ``eln_instructions.txt``: the fixed instructions for the ELN's definition
  (``eln_instructions/<definition>.txt``, picked through ``eln_rules.yaml``);
- ``base_eln.txt``: the rule ELN, with the fields the LLM may fill set to null;
- one text file per context file with text (a PDF's ``pdf_text``), named after the file;
- ``<raw file>_header.txt`` for a Nanonis .sxm.

Without a definition only ``base_eln.txt`` is written; the check step then keeps the
rule ELN. The report lists the files, the files left out and the ``Example:`` values.
"""

import argparse
import copy
import json
import logging
import re
import sys
from pathlib import Path, PurePosixPath

import yaml

from nexus_eln_builder import sxm_header_text

logger = logging.getLogger("nexus_llm_input")

TOOL_DIR = Path(__file__).resolve().parent
INSTRUCTIONS_NAME, BASE_ELN_NAME = "eln_instructions.txt", "base_eln.txt"
BASE_ELN_NOTE = "# Fill only the null fields; keep every other value.\n"


def safe_name(name):
    """A plain file name; same rule as zenodo_download_unzip."""
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


def load_instructions(definition, rules):
    """Instructions text; stops if it misses an allowed value or a field of `eln_rules.yaml`."""
    text = (TOOL_DIR / rules["definitions"][definition]["instructions"]).read_text(encoding="utf-8")
    needed = rules["definitions"][definition]["scan_mode"] + rules["fill"] + rules["replace"]
    missing = [value for value in needed if value not in text]
    if missing:
        sys.exit(f"Instructions for {definition} do not match eln_rules.yaml; missing: {missing}")
    return text


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


def write_input(bundle, eln, raw_file, rules, out_dir):
    """Write the LLM input files to `out_dir`; return the report."""
    raw_path = bundle["raw"]["path"]
    definition = definition_of(eln, rules)
    files = []

    def write(name, text, kind, source=None):
        (out_dir / name).write_text(text, encoding="utf-8")
        files.append({"name": name, "kind": kind, "source": source, "chars": len(text)})

    report = {
        "schema_version": 1,
        "sha256": bundle["raw"]["sha256"],
        "raw_path": raw_path,
        "record_id": bundle["record"]["record_id"],
        "definition": definition,
        "status": "ok" if definition else "no_definition",
        "instructions": None,
        "files": files,
        "left_out": [],
        "example_values": [],
    }
    yaml_text = yaml.safe_dump(base_eln(eln, rules), allow_unicode=True, sort_keys=False)
    write(BASE_ELN_NAME, BASE_ELN_NOTE + yaml_text, "base_eln")
    if not definition:
        logger.warning("No definition for %s; only %s is written", raw_path, BASE_ELN_NAME)
        return report
    instructions = load_instructions(definition, rules)
    write(INSTRUCTIONS_NAME, instructions, "instructions",
          rules["definitions"][definition]["instructions"])
    report["instructions"] = rules["definitions"][definition]["instructions"]
    report["example_values"] = example_values(instructions)

    sources, report["left_out"] = source_files(bundle, raw_file)
    for source in sources:
        write(source["name"], source["text"], source["kind"], source["path"])
    return report


def source_files(bundle, raw_file):
    """The source files the LLM gets, as dicts (name, kind, path, text), and the files left out.

    Also used by nexus_eln_llm_merge to check the quotes against the same texts.
    """
    taken, sources = {INSTRUCTIONS_NAME, BASE_ELN_NAME}, []
    raw_path = bundle["raw"]["path"]
    header = sxm_header_text(raw_file) if raw_file else None
    if header:
        name = unique_name(safe_name(f"{PurePosixPath(raw_path).stem}_header.txt"), taken)
        sources.append({"name": name, "kind": "header", "path": raw_path, "text": header})
    texts, left_out = context_texts(bundle)
    for path, text in texts:
        name = unique_name(safe_name(PurePosixPath(path).name), taken)
        sources.append({"name": name, "kind": "context", "path": path, "text": text})
    return sources, left_out


def main():
    parser = argparse.ArgumentParser(description="Write the LLM input for one raw SPM file.")
    parser.add_argument("bundle", help="Context bundle JSON of one raw file.")
    parser.add_argument("eln", help="ELN YAML of the rule builder.")
    parser.add_argument("--raw", help="The raw file; a Nanonis .sxm header is added as a text file.")
    parser.add_argument("--rules", default=str(TOOL_DIR / "eln_rules.yaml"))
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    bundle = json.loads(Path(args.bundle).read_text(encoding="utf-8"))
    eln = yaml.safe_load(Path(args.eln).read_text(encoding="utf-8")) or {}
    rules = yaml.safe_load(Path(args.rules).read_text(encoding="utf-8"))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report = write_input(bundle, eln, args.raw, rules, out_dir)
    report["total_chars"] = sum(item["chars"] for item in report["files"])
    Path(args.report).write_text(json.dumps(report, indent=2), encoding="utf-8")
    logger.info("LLM input for %s: %d files, %d left out",
                report["raw_path"], len(report["files"]), len(report["left_out"]))


if __name__ == "__main__":
    main()
