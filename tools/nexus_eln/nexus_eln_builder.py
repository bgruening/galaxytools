#!/usr/bin/env python
"""Build an ELN YAML for one raw file from its context bundle.

The bundle (written by the ``zenodo_download_unzip`` tool) holds the Zenodo
record metadata and the text of the raw file's context files. The technique, and so
the NeXus definition, comes from the first of: the tool parameter; the LLM detection
answer, if its quote is found in a source file; the header markers of the rules file;
the raw file's extension. The rules file is the pynxtools plugin's, e.g.
https://github.com/FAIRmat-NFDI/pynxtools-spm

Values are only taken from the record metadata, the context text and checked
answers; nothing is guessed. Fields that were not found are listed in the report.
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path, PurePosixPath

import yaml
from nexus_eln_llm_merge import allowed_value, Dropped, find_quote, norm, parse_answer
from nexus_llm_input import plain_text, set_path, source_files

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


def find_in_context(context, label):
    """First `label: value` line in the context files, with the file it came from."""
    pattern = re.compile(rf"^\s*{label}\s*[:=]\s*(\S.*?)\s*$", re.IGNORECASE | re.MULTILINE)
    for item in context:
        match = pattern.search(item.get("text") or "")
        if match:
            return match.group(1), item["path"]
    return None, None


def definition_by_technique(rules):
    return {value["technique"]: name for name, value in rules["definitions"].items()}


def technique_of(extension, rules):
    """The technique whose `extensions` in the rules hold `extension`; None if none does."""
    for value in rules["definitions"].values():
        if extension in value.get("extensions", []):
            return value["technique"]
    return None


def flat_value(value):
    """A header value as text: `2.1 V` for a value with a unit, JSON for other objects."""
    if isinstance(value, dict) and set(value) == {"value", "unit"}:
        return f"{value['value']} {value['unit']}"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def header_lines(header_path):
    """The header as the markers read it: `key: value` and `channel: <name>` lines, or its text."""
    header = json.loads(Path(header_path).read_text(encoding="utf-8"))
    if "text" in header:
        return header["text"]
    lines = [f"{key}: {flat_value(value)}" for key, value in (header.get("meta") or {}).items()]
    for name, channel in (header.get("channels") or {}).items():
        lines.append(f"channel: {name}")
        meta = channel.get("meta") or {}
        lines += [f"{key}: {flat_value(value)}" for key, value in meta.items()]
    return "\n".join(lines)


def marker_hits(text, rules):
    """{technique: the matched text} for each definition whose `header_markers` match."""
    hits = {}
    for value in rules["definitions"].values():
        for marker in value.get("header_markers") or []:
            try:
                match = re.search(marker, text, flags=re.MULTILINE)
            except re.error as error:
                sys.exit(f"Header marker {marker!r} of {value['technique']} is no regular expression: {error}")
            if match:
                hits[value["technique"]] = match.group(0)
                break
    return hits


def llm_detection(answer_text, sources, rules):
    """The checked `technique` and `scan_mode` of the detection answer, with a status.

    A value is kept when its quote is found in a source file; a technique must be one of the rules.
    """
    result = {"status": "no_answer", "fields": {}}
    if answer_text is None:
        return result
    values, evidence = parse_answer(answer_text)
    if values is None:
        result["status"] = "unparsable"
        return result
    result["status"] = "checked"
    techniques = list(definition_by_technique(rules))
    for field in ("technique", "scan_mode"):
        value = values.get(field)
        if value in (None, "", "null"):
            continue
        entry = {"value": value, "notes": []}
        try:
            if not isinstance(value, str):
                raise Dropped("not_text")
            find_quote(evidence.get(field), sources, entry)
            if field == "technique":
                known = [name for name in techniques if name.lower() == norm(value).lower()]
                if not known:
                    raise Dropped("unknown_technique")
                entry["value"] = known[0]
            entry["status"] = "kept"
        except Dropped as reason:
            entry.update(status="dropped", reason=str(reason))
        result["fields"][field] = entry
    return result


def pick_technique(parameter, llm, markers, extension, rules):
    """(technique, source, conflict) by the order: parameter, LLM, one header marker, extension."""
    llm_technique = llm["fields"].get("technique", {})
    candidates = []
    if parameter != "auto":
        candidates.append((parameter, "tool parameter"))
    if llm_technique.get("status") == "kept":
        candidates.append((llm_technique["value"], f"LLM detection ({llm_technique['file']})"))
    if len(markers) == 1:
        (technique, match), = markers.items()
        candidates.append((technique, f"header marker {match!r}"))
    by_extension = technique_of(extension, rules)
    if by_extension:
        candidates.append((by_extension, f"extension {extension}"))
    if not candidates:
        return None, None, False
    technique, source = candidates[0]
    conflict = any(other != technique for other, _ in candidates[1:])
    return technique, source, conflict


def build_eln(bundle, technique, rules, header_path, answer_text=None):
    """The ELN content and a report of where each value came from."""
    record, raw = bundle["record"], bundle["raw"]
    extension = PurePosixPath(raw["path"]).suffix.lower()
    llm = llm_detection(answer_text, source_files(bundle, header_path)[0], rules)
    markers = marker_hits(header_lines(header_path), rules)
    parameter = technique
    technique, technique_source, conflict = pick_technique(parameter, llm, markers, extension, rules)
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
    llm_scan_mode = llm["fields"].get("scan_mode", {})
    if llm_scan_mode.get("status") == "kept" and technique and "scan_mode" not in sources:
        allowed = rules["definitions"][eln["definition"]]["scan_mode"]
        try:
            eln["scan_mode"] = allowed_value(llm_scan_mode["value"], allowed, llm_scan_mode["notes"])
            sources["scan_mode"] = f"LLM detection ({llm_scan_mode['file']})"
        except Dropped as reason:
            llm_scan_mode.update(status="dropped", reason=f"{reason} for {technique}")
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
        "detection": {
            "parameter": None if parameter == "auto" else parameter,
            "llm": llm,
            "header_markers": markers,
            "extension": technique_of(extension, rules),
            "conflict": conflict,
        },
        "sources": sources,
        "missing": missing,
    }
    return eln, report


def main():
    parser = argparse.ArgumentParser(description="Build an ELN YAML from a context bundle.")
    parser.add_argument("bundle", help="Context bundle JSON of one raw file.")
    parser.add_argument("--rules", required=True, help="Rules YAML of the pynxtools plugin.")
    parser.add_argument("--technique", default="auto", help="auto, or a technique of the rules file.")
    parser.add_argument("--header", required=True, help="Header JSON of the raw file (raw_file_header).")
    parser.add_argument("--answer", help="The llm_hub answer of the detection job (Markdown).")
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
    answer = Path(args.answer).read_text(encoding="utf-8", errors="replace") if args.answer else None
    # A skipped llm_hub step (LLM off in the workflow) gives a dataset that holds only "null".
    if answer is not None and answer.strip() == "null":
        answer = None
    eln, report = build_eln(bundle, args.technique, rules, args.header, answer)
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
