#!/usr/bin/env python
"""Collect the conversion results of a harvest run and build the NOMAD projects.

Reads the manifest of the ``zenodo_download_unzip`` tool and the collections of
the later steps, all keyed by the same element identifiers, and writes:

- ``results.json``: conversion status per raw file (by SHA-256), with
  ``nexus_path``, the NeXus file's path in its project zip;
- one project zip per record with at least one converted raw file, laid out as
  ``<version DOI>/<raw file stem>/`` (raw file, NeXus file, ``eln_data.yaml``)
  and ``<version DOI>/Complementaries/`` (the record's other files);
- one log per record version (named after its DOI) and an HTML summary.

``results.json`` and the zips are read by the NOMAD side,
nomad-eosc-galaxy-actions (``manifest.py``, ``nomad_projects.py``):
https://github.com/FAIRmat-NFDI/nomad-eosc-galaxy-actions
"""

import argparse
import html
import json
import logging
import re
import sys
import zipfile
from pathlib import Path, PurePosixPath

logger = logging.getLogger("nexus_report_builder")

# Folder for the record files that are not converted (README, notes, ...).
COMPLEMENTARIES = "Complementaries"
ELN_NAME = "eln_data.yaml"
ZIP_DATE = (2020, 1, 1, 0, 0, 0)  # fixed, so the same content gives the same zip
MAX_WARNINGS_IN_LOG = 20


def safe_name(name):
    return re.sub(r"[^\w.\-]", "_", name)


def element_id(record_id, sha256, path):
    """Same rule as in zenodo_download_unzip.py (the tools ship separately)."""
    return f"{record_id}-{sha256[:12]}-{safe_name(PurePosixPath(path).name)}"


def read_dir(folder):
    """Files of a folder by name (the element identifier); empty if not given."""
    if not folder or not Path(folder).is_dir():
        return {}
    return {path.name: path for path in Path(folder).iterdir()}


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def version_name(record_id, version_doi):
    """`10.5281/zenodo.1` → `10.5281_zenodo.1`; `zenodo_<record ID>` without a DOI."""
    return version_doi.replace("/", "_") if version_doi else f"zenodo_{record_id}"


def log_name(record_id, version_doi):
    return f"{version_name(record_id, version_doi)}.log"


def folder_names(files):
    """A unique folder name per raw file: its stem, else `<stem>_<ext>`, else `<stem>_<ext>_<sha6>`."""

    def candidates(file):
        path = PurePosixPath(file["path"])
        stem, ext = safe_name(path.stem), safe_name(path.suffix.lstrip("."))
        with_ext = f"{stem}_{ext}" if ext else stem
        return [stem, with_ext, f"{with_ext}_{file['sha256'][:6]}"]

    names = {}
    for level in range(3):
        counts = {}
        for file in files:
            if file["sha256"] not in names:
                name = candidates(file)[level]
                counts[name] = counts.get(name, 0) + 1
        for file in files:
            name = candidates(file)[level]
            if file["sha256"] not in names and (counts[name] == 1 or level == 2):
                names[file["sha256"]] = name
    return names


def missing_warnings(converter_log):
    if converter_log is None:
        return []
    return [
        line.strip()
        for line in converter_log.read_text(encoding="utf-8", errors="replace").splitlines()
        if "hasn't been supplied" in line
    ]


class Run:
    def __init__(self, args):
        self.manifest = read_json(args.manifest)
        self.statuses = {name: read_json(path) for name, path in read_dir(args.statuses_dir).items()}
        self.nexus = read_dir(args.nexus_dir)
        self.elns = read_dir(args.eln_dir)
        self.eln_reports = {name: read_json(p) for name, p in read_dir(args.eln_report_dir).items()}
        self.raw = read_dir(args.raw_dir)
        self.context = read_dir(args.context_dir)
        self.converter_logs = read_dir(args.converter_log_dir)
        self.records = self._records()

    def _records(self):
        """Record ID → its archives, files and download failures, in manifest order."""
        records = {}

        def record(record_id, version_doi=None):
            entry = records.setdefault(record_id, {
                "record_id": record_id, "version_doi": None,
                "archives": [], "files": [], "download_failures": [],
            })
            entry["version_doi"] = entry["version_doi"] or version_doi
            return entry

        for archive in self.manifest.get("archives", []):
            record(archive["record_id"], archive.get("version_doi"))["archives"].append(archive)
        for file in self.manifest.get("files", []):
            record(file["record_id"], file.get("version_doi"))["files"].append(file)
        for failure in self.manifest.get("download_failures", []):
            record(failure["record_id"])["download_failures"].append(failure)
        return records

    def status_of(self, file):
        return self.statuses.get(element_id(file["record_id"], file["sha256"], file["path"]))

    def build_project(self, record, projects_dir):
        """Write the record's project zip; return {sha256: nexus_path} of its NeXus files."""
        converted = [
            file for file in record["files"]
            if file["role"] == "raw" and file["action"] == "process"
            and (self.status_of(file) or {}).get("status") == "success"
        ]
        if not converted:
            return {}
        root = version_name(record["record_id"], record["version_doi"])
        folders = folder_names(converted)
        members, nexus_paths = {}, {}
        for file in converted:
            name = element_id(file["record_id"], file["sha256"], file["path"])
            folder = f"{root}/{folders[file['sha256']]}"
            raw_name = PurePosixPath(file["path"]).name
            members[f"{folder}/{raw_name}"] = self.raw[name]
            nexus_path = f"{folder}/{PurePosixPath(raw_name).stem}.nxs"
            members[nexus_path] = self.nexus[name]
            nexus_paths[file["sha256"]] = nexus_path
            if name in self.elns:
                members[f"{folder}/{ELN_NAME}"] = self.elns[name]
        # Every other file the download step kept (context files, images, other
        # data formats); files matched by an ignore pattern are not in self.context.
        for file in record["files"]:
            if file["role"] == "raw":
                continue
            source = self.context.get(element_id(file["record_id"], file["sha256"], file["path"]))
            if source:
                members[f"{root}/{COMPLEMENTARIES}/{file['path']}"] = source
        zip_path = projects_dir / f"{record['record_id']}.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
            for arcname in sorted(members):
                info = zipfile.ZipInfo(arcname, ZIP_DATE)
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, members[arcname].read_bytes())
        return nexus_paths

    def record_log(self, record, nexus_paths):
        lines = [
            f"Zenodo record {record['record_id']}",
            f"DOI: {record['version_doi'] or '-'}",
            f"Run: {self.manifest.get('run_id')}  Pipeline version: {self.manifest.get('tool_version')}",
            f"NOMAD project zip: {'yes' if nexus_paths else 'no (no converted raw file)'}",
            "",
            "Archives:",
        ]
        lines += [
            f"  {a['key']}: {a['status']}" + (f" ({a['reason']})" if a.get("reason") else "")
            for a in record["archives"]
        ] or ["  none"]
        if record["download_failures"]:
            lines.append("Download failures:")
            lines += [f"  {f['key']}: {f['reason']}" for f in record["download_failures"]]
        lines += ["", "Raw files:"]
        raw_files = [file for file in record["files"] if file["role"] == "raw"]
        for file in raw_files:
            name = element_id(file["record_id"], file["sha256"], file["path"])
            if file["action"] != "process":
                lines.append(f"  {file['path']}: {file['action']}")
                continue
            status = self.status_of(file)
            if status is None:
                lines.append(f"  {file['path']}: failed (no_result)")
                continue
            outcome = status["status"] + (f": {status['error']}" if status.get("error") else "")
            lines.append(f"  {file['path']}: {outcome}")
            lines.append(f"    sha256 {file['sha256']}, application definition {status.get('nxdl')}")
            if file["sha256"] in nexus_paths:
                lines.append(f"    NeXus file: {nexus_paths[file['sha256']]}")
            missing = self.eln_reports.get(name, {}).get("missing")
            if missing:
                lines.append(f"    ELN fields not found: {', '.join(missing)}")
            warnings = missing_warnings(self.converter_logs.get(name))
            if warnings:
                lines.append(f"    Validation: {len(warnings)} required NeXus concepts not supplied")
                lines += [f"      {w}" for w in warnings[:MAX_WARNINGS_IN_LOG]]
            if file["context_paths"]:
                lines.append(f"    Context files: {', '.join(file['context_paths'])}")
        if not raw_files:
            lines.append("  none")
        return "\n".join(lines) + "\n"

    def counts(self, record):
        raw = [f for f in record["files"] if f["role"] == "raw"]
        processed = [f for f in raw if f["action"] == "process"]
        success = sum((self.status_of(f) or {}).get("status") == "success" for f in processed)
        return {
            "raw files": len(raw),
            "converted": success,
            "failed": len(processed) - success,
            "already known": sum(f["action"] == "already_known" for f in raw),
            "archives not extracted": sum(a["status"] != "extracted" for a in record["archives"]),
            "download failures": len(record["download_failures"]),
        }


def summary_html(run, rows):
    columns = ["record", "DOI", *rows[0][1].keys(), "project"] if rows else []
    body = "".join(
        "<tr>" + "".join(
            f"<td>{html.escape(str(value))}</td>"
            for value in [record["record_id"], record["version_doi"] or "-", *counts.values(), project]
        ) + "</tr>"
        for record, counts, project in rows
    )
    head = "".join(f"<th>{html.escape(c)}</th>" for c in columns)
    table = f"<table><tr>{head}</tr>{body}</table>" if rows else "<p>No records in this run.</p>"
    return (
        "<!doctype html><html><head><meta charset='utf-8'><title>Harvest summary</title>"
        "<style>body{font-family:sans-serif}table{border-collapse:collapse}"
        "td,th{border:1px solid #999;padding:4px 8px;text-align:left}</style></head><body>"
        f"<h1>Harvest summary</h1><p>Run {html.escape(str(run.manifest.get('run_id')))}, "
        f"pipeline version {html.escape(str(run.manifest.get('tool_version')))}</p>"
        f"{table}</body></html>\n"
    )


def main():
    parser = argparse.ArgumentParser(description="Build results.json, project zips and logs.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--statuses-dir", required=True)
    parser.add_argument("--nexus-dir", required=True)
    parser.add_argument("--eln-dir", required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--context-dir", required=True)
    parser.add_argument("--eln-report-dir", default=None)
    parser.add_argument("--converter-log-dir", default=None)
    parser.add_argument("--results", required=True)
    parser.add_argument("--projects-dir", required=True)
    parser.add_argument("--logs-dir", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    run = Run(args)
    projects_dir, logs_dir = Path(args.projects_dir), Path(args.logs_dir)
    projects_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    nexus_paths, rows = {}, []
    for record in run.records.values():
        paths = run.build_project(record, projects_dir)
        nexus_paths.update(paths)
        (logs_dir / log_name(record["record_id"], record["version_doi"])).write_text(
            run.record_log(record, paths), encoding="utf-8"
        )
        rows.append((record, run.counts(record), "yes" if paths else "no"))

    results = {
        "schema_version": 1,
        "run_id": run.manifest["run_id"],
        "files": [
            {
                "sha256": status["sha256"],
                "status": status["status"],
                "error": status.get("error"),
                "nexus_path": nexus_paths.get(status["sha256"]) if status["status"] == "success" else None,
            }
            for _, status in sorted(run.statuses.items())
        ],
    }
    with open(args.results, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    Path(args.summary).write_text(summary_html(run, rows), encoding="utf-8")
    logger.info("%d records, %d results, %d project zips",
                len(run.records), len(results["files"]), sum(1 for _ in projects_dir.iterdir()))


if __name__ == "__main__":
    main()
