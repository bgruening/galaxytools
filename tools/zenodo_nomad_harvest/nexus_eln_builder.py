#!/usr/bin/env python
"""Build an ELN YAML for one raw SPM file from its context bundle.

The bundle (written by the ``zenodo_download_unzip`` tool) holds the Zenodo
record metadata and the text of the raw file's context files. The ELN fills
the fields that ``pynxtools-spm`` needs besides the raw file, following the
ELN examples of pynxtools-spm: https://github.com/FAIRmat-NFDI/pynxtools-spm

Values are only taken from the record metadata or found in the context text;
nothing is guessed. Fields that were not found are listed in the report.
"""

import argparse
import html
import json
import logging
import re
import sys
from pathlib import PurePosixPath

import yaml

logger = logging.getLogger("nexus_eln_builder")

# Technique per raw-file extension, as pynxtools-spm reads them.
TECHNIQUE_BY_EXTENSION = {
    ".sxm": "STM", ".sm4": "STM", ".dat": "STS",
    ".spm": "AFM", ".txt": "AFM", ".flt": "AFM",
}
DEFINITION_BY_TECHNIQUE = {"STM": "NXstm", "STS": "NXsts", "AFM": "NXafm"}

# "Label: value" lines searched in the context text; the first match wins.
CONTEXT_FIELDS = {
    "Sample.name": r"sample(?:\s+name)?",
    "Sample.description": r"sample\s+description",
    "scan_mode": r"scan\s*mode",
    "User.name": r"(?:user|operator|measured\s+by)",
}
# Fields the SPM reader expects and that only an ELN can give.
EXPECTED_FIELDS = ["User.name", "Sample.name", "scan_mode"]


def plain_text(description_html):
    """Text of Zenodo's HTML description."""
    text = re.sub(r"<br\s*/?>|</p>", "\n", description_html or "", flags=re.IGNORECASE)
    text = html.unescape(re.sub(r"<[^>]+>", "", text))
    return re.sub(r"\n\s*\n+", "\n\n", text).strip()


def find_in_context(context, label):
    """First `label: value` line in the context files, with the file it came from."""
    pattern = re.compile(rf"^\s*{label}\s*[:=]\s*(\S.*?)\s*$", re.IGNORECASE | re.MULTILINE)
    for item in context:
        match = pattern.search(item.get("text") or "")
        if match:
            return match.group(1), item["path"]
    return None, None


def set_path(eln, dotted, value):
    node = eln
    *parents, last = dotted.split(".")
    for key in parents:
        node = node.setdefault(key, {})
    node[last] = value


def build_eln(bundle, technique):
    """The ELN content and a report of where each value came from."""
    record, raw = bundle["record"], bundle["raw"]
    extension = PurePosixPath(raw["path"]).suffix.lower()
    if technique == "auto":
        technique = TECHNIQUE_BY_EXTENSION.get(extension)
        technique_source = f"extension {extension}"
    else:
        technique_source = "tool parameter"
    sources = {}
    eln = {}
    if technique:
        eln["experiment_technique"] = technique
        eln["definition"] = DEFINITION_BY_TECHNIQUE[technique]
        sources["experiment_technique"] = technique_source
    record_url = f"https://zenodo.org/records/{record['record_id']}"
    eln["identifier_collection"] = (
        f"https://doi.org/{record['version_doi']}" if record.get("version_doi") else record_url
    )
    eln["identifier_experiment"] = PurePosixPath(raw["path"]).name
    description = plain_text(record.get("description_html"))
    title = record.get("title") or ""
    if title or description:
        eln["experiment_description"] = "\n\n".join(part for part in (title, description) if part)
        sources["experiment_description"] = "record title and description"
    cite = {"url": [record_url]}
    if record.get("creators"):
        cite["author"] = list(record["creators"])
    if record.get("version_doi"):
        cite["doi"] = [record["version_doi"]]
    if description:
        cite["description"] = [description]
    eln["citeID"] = cite
    sources["citeID"] = "record metadata"
    for field, label in CONTEXT_FIELDS.items():
        value, path = find_in_context(bundle.get("context", []), label)
        if value is not None:
            set_path(eln, field, [value] if field == "User.name" else value)
            sources[field] = path
    missing = [field for field in EXPECTED_FIELDS if field not in sources]
    if not technique:
        missing.insert(0, "experiment_technique")
    report = {
        "schema_version": 1,
        "sha256": raw["sha256"],
        "raw_path": raw["path"],
        "record_id": record["record_id"],
        "technique": technique,
        "sources": sources,
        "missing": missing,
    }
    return eln, report


def main():
    parser = argparse.ArgumentParser(description="Build an ELN YAML from a context bundle.")
    parser.add_argument("bundle", help="Context bundle JSON of one raw file.")
    parser.add_argument("--technique", default="auto", choices=["auto", *DEFINITION_BY_TECHNIQUE])
    parser.add_argument("--eln", required=True, help="Output ELN YAML.")
    parser.add_argument("--report", required=True, help="Output JSON: value sources and missing fields.")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    with open(args.bundle, encoding="utf-8") as handle:
        bundle = json.load(handle)
    eln, report = build_eln(bundle, args.technique)
    header = (
        f"# ELN built from Zenodo record {report['record_id']} for {report['raw_path']}.\n"
        f"# Not found in the record or its context files: {', '.join(report['missing']) or 'none'}.\n"
    )
    with open(args.eln, "w", encoding="utf-8") as handle:
        handle.write(header)
        yaml.safe_dump(eln, handle, allow_unicode=True, sort_keys=False)
    with open(args.report, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    logger.info("ELN for %s; missing: %s", report["raw_path"], report["missing"])


if __name__ == "__main__":
    main()
