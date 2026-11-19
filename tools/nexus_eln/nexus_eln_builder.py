#!/usr/bin/env python
"""Build an ELN YAML for one raw file from its context bundle.

The bundle (written by the ``zenodo_download_unzip`` tool) holds the Zenodo
record metadata and the text of the raw file's context files. The technique and
NeXus definition come from the raw file's extension, as listed in the rules file
(``definitions.<name>.extensions``) of the pynxtools plugin, e.g.
https://github.com/FAIRmat-NFDI/pynxtools-spm

Values are only taken from the record metadata and the context text; nothing is
guessed. Fields that were not found are listed in the report.
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

# "Label: value" lines searched in the context text; the first match wins.
CONTEXT_FIELDS = {
    "Sample.name": r"sample(?:\s+name)?",
    "Sample.description": r"sample\s+description",
    "scan_mode": r"scan\s*mode",
}
# Fields the SPM reader expects and that only an ELN can give.
EXPECTED_FIELDS = ["User", "Sample.name", "scan_mode"]
# The ELN's only user is the person who runs the harvest (from the work list), never a
# person named in the dataset files; the Zenodo creators are in citeID.
USER_KEYS = ("name", "email", "affiliation", "address")
USER_ROLE = "harvested and converted the data"


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


def definition_by_technique(rules):
    return {value["technique"]: name for name, value in rules["definitions"].items()}


def technique_of(extension, rules):
    """The technique whose `extensions` in the rules hold `extension`; None if none does."""
    for value in rules["definitions"].values():
        if extension in value.get("extensions", []):
            return value["technique"]
    return None


def build_eln(bundle, technique, rules):
    """The ELN content and a report of where each value came from.

    The technique parameter, unless "auto", wins over the extension rule.
    """
    record, raw = bundle["record"], bundle["raw"]
    extension = PurePosixPath(raw["path"]).suffix.lower()
    if technique != "auto":
        technique_source = "tool parameter"
    else:
        technique = technique_of(extension, rules)
        technique_source = f"extension {extension}"
    sources = {}
    eln = {}
    if technique:
        eln["experiment_technique"] = technique
        eln["definition"] = definition_by_technique(rules)[technique]
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
            set_path(eln, field, value)
            sources[field] = path
    user = {key: bundle["user"][key] for key in USER_KEYS if (bundle.get("user") or {}).get(key)}
    if "name" in user:
        # A list of users: pynxtools-spm writes USER[user_1].
        eln["User"] = [user | {"role": USER_ROLE}]
        sources["User"] = "the harvest's user (work list)"
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
    parser.add_argument("--rules", required=True, help="Rules YAML of the pynxtools plugin.")
    parser.add_argument("--technique", default="auto", help="auto, or a technique of the rules file.")
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
    with open(args.rules, encoding="utf-8") as handle:
        rules = yaml.safe_load(handle)
    techniques = definition_by_technique(rules)
    if args.technique not in ("auto", *techniques):
        sys.exit(f"Unknown technique {args.technique!r}; the rules file has: {', '.join(techniques)}")
    eln, report = build_eln(bundle, args.technique, rules)
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
