#!/usr/bin/env python
"""Download the files of a work list from Zenodo, extract them safely and
write ``manifest.json``.

The work list (``work_list.json``) and the manifest are the contract with the
NOMAD side, nomad-eosc-galaxy-actions (``work_list.py``, ``manifest.py``):
https://github.com/FAIRmat-NFDI/nomad-eosc-galaxy-actions

For each record the script downloads every listed file with ``zenodo_get``
(a separate program, not imported), checks its MD5 against the work list,
extracts ``.zip`` archives, computes SHA-256 checksums and classifies every
file as ``raw``, ``context`` or ``ignored``. Raw files to convert are written
to ``--raw-dir``, one context bundle per raw file (same file name) to
``--bundle-dir``, and every other file not matched by ``ignore_patterns`` to
``--context-dir``, for the project's ``Complementaries`` folder.
The text of a PDF context file is read with ``pdftotext`` (Poppler) into the
bundle's ``pdf_text``; a PDF without text gets ``pdf_text_error`` instead.

Failures (download, MD5, unsafe archive) are recorded in the manifest and the
script still exits with 0, so the run goes on with the other files.
"""

import argparse
import functools
import glob
import hashlib
import json
import logging
import re
import shutil
import subprocess
import sys
import zipfile
from fnmatch import fnmatch
from pathlib import Path, PurePosixPath

logger = logging.getLogger("zenodo_download_unzip")

ARCHIVE_EXTENSIONS = (".zip",)
# Archives inside an archive are not extracted (yet).
NESTED_ARCHIVE_EXTENSIONS = (".zip", ".tar", ".tar.gz", ".tgz", ".gz", ".bz2", ".xz", ".7z", ".rar")
# Context files whose text goes into the bundle, for the ELN builder.
TEXT_EXTENSIONS = (".txt", ".md", ".csv", ".json", ".yaml", ".yml", ".rtf")
MAX_TEXT_BYTES = 1_000_000
# PDF context files: their text goes into the bundle as `pdf_text`, read by pdftotext (Poppler).
MAX_PDF_BYTES = 100_000_000
PDF_TIMEOUT_SECONDS = 120
SUPPORTED_SCHEMA_VERSION = 1


class UncompressedLimitError(Exception):
    pass


def file_hashes(path):
    """SHA-256 and MD5 of a file, read once."""
    sha256, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(2**20), b""):
            sha256.update(chunk)
            md5.update(chunk)
    return sha256.hexdigest(), md5.hexdigest()


def safe_name(name):
    """A plain file name; also used by the tool XML for prefetched files."""
    return re.sub(r"[^\w.\-]", "_", name)


def element_id(record_id, sha256, path):
    """Collection element identifier; ends with the file name, for its extension."""
    return f"{record_id}-{sha256[:12]}-{safe_name(PurePosixPath(path).name)}"


def has_extension(name, extensions):
    return name.lower().endswith(tuple(ext.lower() for ext in extensions))


def is_ignored(path, patterns):
    """True if a pattern matches the path or the path below any of its folders."""
    parts = PurePosixPath(path).parts
    tails = ["/".join(parts[i:]) for i in range(len(parts))]
    return any(fnmatch(tail, pattern) for tail in tails for pattern in patterns)


def classify(path, work_list):
    """Role and skip reason of a file.

    An ignored file that no ignore pattern matches is still kept for the project
    (reason ``complementary`` or ``unsupported_archive``), see ``is_complementary``.
    """
    if is_ignored(path, work_list.get("ignore_patterns", [])):
        return "ignored", "ignore_pattern"
    if has_extension(path, work_list.get("raw_extensions", [])):
        return "raw", None
    if has_extension(path, work_list.get("context_extensions", [])):
        return "context", None
    if has_extension(path, NESTED_ARCHIVE_EXTENSIONS):
        return "ignored", "unsupported_archive"
    return "ignored", "complementary"


def is_complementary(row):
    """True for files that go to the project's Complementaries folder: all non-raw
    files except those matched by an ignore pattern."""
    return row["role"] != "raw" and row["reason"] != "ignore_pattern"


def unsafe_member(info):
    """True for absolute paths, `..`, drive letters and symbolic links ("zip slip")."""
    name = info.filename.replace("\\", "/")
    is_symlink = (info.external_attr >> 16) & 0o170000 == 0o120000
    return (
        name.startswith("/")
        or re.match(r"^[A-Za-z]:", name) is not None
        or ".." in PurePosixPath(name).parts
        or is_symlink
    )


def check_archive(archive, max_uncompressed_bytes):
    """Why the archive must not be extracted, or None."""
    infos = archive.infolist()
    if any(unsafe_member(info) for info in infos):
        return "unsafe_path"
    if any(info.flag_bits & 0x1 for info in infos):  # encrypted member
        return "unsupported_archive"
    # Checked on the zip's own index first, against "zip bombs".
    if sum(info.file_size for info in infos) > max_uncompressed_bytes:
        return "uncompressed_limit"
    return None


def extract(archive, target, max_uncompressed_bytes):
    """Extract all files; also count the real bytes, in case the index lies."""
    target = target.resolve()
    total = 0
    members = []
    for info in archive.infolist():
        if info.is_dir():
            continue
        destination = (target / info.filename).resolve()
        if not destination.is_relative_to(target):
            raise ValueError(f"unsafe path {info.filename!r}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with archive.open(info) as source, open(destination, "wb") as sink:
            for chunk in iter(lambda: source.read(2**20), b""):
                total += len(chunk)
                if total > max_uncompressed_bytes:
                    raise UncompressedLimitError()
                sink.write(chunk)
        members.append((info.filename, destination))
    return members


def fetch(record_id, entry, prefetched, download_dir, zenodo_get):
    """Local path of a listed file and the failure reason (or None)."""
    key = entry["key"]
    path = prefetched.get(safe_name(key))
    if path is None:
        out = download_dir / record_id
        out.mkdir(parents=True, exist_ok=True)
        # One exact file name per call; `glob.escape` quotes `*`, `?` and `[`.
        command = [
            zenodo_get, "-r", record_id, "-g", glob.escape(key),
            "-o", str(out), "-k", "-v", "1",
        ]
        logger.info("Downloading %s from record %s", key, record_id)
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            logger.warning("zenodo_get failed for %s: %s", key, completed.stderr[-2000:])
        path = out / key
        if not path.is_file():
            return None, "download_failed"
    sha256, md5 = file_hashes(path)
    if entry.get("md5") and entry["md5"].lower() != md5:
        return path, "md5_mismatch"
    return path, None


def file_row(record, source, path, local_path, archive_key=None):
    sha256, md5 = file_hashes(local_path)
    return {
        "concept_record_id": record["concept_record_id"],
        "version_doi": record.get("version_doi"),
        "record_id": record["record_id"],
        "source": source,
        "archive_key": archive_key,
        "path": path,
        "size": local_path.stat().st_size,
        "sha256": sha256,
        "md5": md5 if source == "loose" else None,
        "role": None,
        "action": None,
        "reason": None,
        "context_paths": [],
        "_local": local_path,
    }


def assign_actions(rows, record, work_list):
    """Role and action per file; a raw file is processed once per record."""
    known = set(record.get("known_sha256", []))
    seen = set()
    for row in rows:
        row["role"], row["reason"] = classify(row["path"], work_list)
        if row["role"] == "ignored":
            row["action"] = "skipped"
        elif row["role"] == "context":
            row["action"] = "process"
        elif row["sha256"] in known or row["sha256"] in seen:
            row["action"] = "already_known"
        else:
            row["action"] = "process"
            seen.add(row["sha256"])


def assign_context(rows):
    """Folder rule: context files in the raw file's folder and every parent folder
    of the same archive, plus all loose context files of the record."""
    context = [row for row in rows if row["role"] == "context"]
    for row in rows:
        if row["role"] != "raw":
            continue
        folders = PurePosixPath(row["path"]).parents
        row["context_paths"] = sorted(
            other["path"]
            for other in context
            if other["source"] == "loose"
            or (
                other["archive_key"] == row["archive_key"]
                and PurePosixPath(other["path"]).parent in folders
            )
        )


def drop_rtf_groups(text):
    """Remove header groups (font and colour tables, `{\\*...}`), nested braces included."""
    pattern = re.compile(r"\{\\(?:\*|fonttbl|colortbl|stylesheet|info)")
    while match := pattern.search(text):
        depth, end = 0, match.start()
        for end in range(match.start(), len(text)):
            depth += {"{": 1, "}": -1}.get(text[end], 0)
            if depth == 0:
                break
        text = text[:match.start()] + text[end + 1:]
    return text


def strip_rtf(text):
    """Rough plain text of an RTF file: drops header groups, control words and braces."""
    text = drop_rtf_groups(text)
    text = re.sub(r"\\'[0-9a-fA-F]{2}", "", text)
    text = re.sub(r"\\[a-zA-Z]+-?\d* ?", "", text)
    return re.sub(r"[{}]", "", text).strip()


def context_text(row):
    """Text of a text-like context file, or None (for example a PDF)."""
    if not has_extension(row["path"], TEXT_EXTENSIONS) or row["size"] > MAX_TEXT_BYTES:
        return None
    text = row["_local"].read_text(encoding="utf-8", errors="replace")
    return strip_rtf(text) if has_extension(row["path"], (".rtf",)) else text


def pdf_text(row, pdftotext="pdftotext"):
    """Text of a PDF file and the failure reason (one of them is None).

    No `-layout`: it mixes the two columns of a paper. The LLM prompts were tested on this text.
    """
    if row["size"] > MAX_PDF_BYTES:
        return None, "too_large"
    command = [pdftotext, "-enc", "UTF-8", str(row["_local"]), "-"]
    try:
        completed = subprocess.run(
            command, capture_output=True, timeout=PDF_TIMEOUT_SECONDS, check=False
        )
    except subprocess.TimeoutExpired:
        return None, "timeout"
    if completed.returncode != 0:
        stderr = completed.stderr.decode("utf-8", errors="replace")[-2000:]
        logger.warning("pdftotext failed for %s: %s", row["path"], stderr)
        return None, "unreadable"
    text = completed.stdout.decode("utf-8", errors="replace")
    # Only page breaks or spaces: a scanned PDF without a text layer.
    if not text.strip():
        return None, "no_text"
    return text, None


def context_item(row):
    """Bundle entry of one context file; a PDF also gets `pdf_text` and `pdf_text_error`."""
    item = {"path": row["path"], "sha256": row["sha256"], "text": context_text(row)}
    if has_extension(row["path"], (".pdf",)):
        item["pdf_text"], item["pdf_text_error"] = pdf_text(row)
        if item["pdf_text_error"]:
            logger.warning("No text from PDF %s: %s", row["path"], item["pdf_text_error"])
    return item


def write_outputs(rows, record, work_list, raw_dir, bundle_dir, context_dir):
    by_path = {row["path"]: row for row in rows if row["role"] == "context"}
    # Many raw files share a context file; read each one (a PDF can be slow) only once.
    context_items = functools.cache(lambda path: context_item(by_path[path]))
    record_meta = {
        key: record.get(key)
        for key in (
            "record_id", "version_doi", "concept_record_id", "concept_doi",
            "license", "title", "creators", "description_html",
        )
    }
    for row in rows:
        if is_complementary(row):
            # Same name means same content, so a repeated file is copied once.
            name = element_id(record["record_id"], row["sha256"], row["path"])
            if not (context_dir / name).exists():
                shutil.copyfile(row["_local"], context_dir / name)
        if row["role"] != "raw" or row["action"] != "process":
            continue
        name = element_id(record["record_id"], row["sha256"], row["path"])
        shutil.copyfile(row["_local"], raw_dir / name)
        bundle = {
            "schema_version": 1,
            "run_id": work_list["run_id"],
            # The person who runs the harvest; the ELN's only user. Optional.
            "user": work_list.get("user"),
            "record": record_meta,
            "raw": {key: row[key] for key in ("path", "sha256", "size", "source", "archive_key")},
            "context": [context_items(path) for path in row["context_paths"]],
        }
        (bundle_dir / name).write_text(json.dumps(bundle, indent=2), encoding="utf-8")


def process_record(record, work_list, prefetched, dirs, zenodo_get):
    """Download, extract and classify one record; return its manifest parts."""
    limits = work_list["limits"]
    rows, archives, failures = [], [], []
    for entry in record.get("download", []):
        key = entry["key"]
        archive_row = {
            "concept_record_id": record["concept_record_id"],
            "record_id": record["record_id"],
            "version_doi": record.get("version_doi"),
            "key": key,
            "md5": entry.get("md5"),
            "status": "extracted",
            "reason": None,
        }
        is_archive = has_extension(key, ARCHIVE_EXTENSIONS)
        if is_archive and entry.get("size", 0) > limits["max_archive_bytes"]:
            archives.append(archive_row | {"status": "skipped", "reason": "too_large"})
            continue
        path, reason = fetch(record["record_id"], entry, prefetched, dirs["download"], zenodo_get)
        if not is_archive:
            if reason:
                # No row in `files`: without trusted content there is no SHA-256.
                failures.append({"record_id": record["record_id"], "key": key, "reason": reason})
            else:
                rows.append(file_row(record, "loose", key, path))
            continue
        if reason:
            archives.append(archive_row | {"status": "failed", "reason": reason})
            continue
        target = dirs["extract"] / record["record_id"] / safe_name(key)
        try:
            with zipfile.ZipFile(path) as archive:
                reason = check_archive(archive, limits["max_uncompressed_bytes"])
                members = [] if reason else extract(
                    archive, target, limits["max_uncompressed_bytes"]
                )
        except UncompressedLimitError:
            reason = "uncompressed_limit"
        except (zipfile.BadZipFile, NotImplementedError, ValueError) as error:
            logger.warning("Cannot extract %s: %s", key, error)
            reason = "unsupported_archive"
        if reason:
            shutil.rmtree(target, ignore_errors=True)  # nothing extracted is kept
            archives.append(archive_row | {"status": "failed", "reason": reason})
            continue
        archives.append(archive_row)
        rows.extend(
            file_row(record, "archive", name, local, archive_key=key)
            for name, local in members
        )
    assign_actions(rows, record, work_list)
    assign_context(rows)
    write_outputs(rows, record, work_list, dirs["raw"], dirs["bundle"], dirs["context"])
    files = [{k: v for k, v in row.items() if k != "_local"} for row in rows]
    return archives, files, failures


def main():
    parser = argparse.ArgumentParser(description="Download, extract and classify a work list.")
    parser.add_argument("work_list", help="work_list.json from the NOMAD side.")
    parser.add_argument("--prefetched-dir", default=None,
                        help="Files already in Galaxy, named after their Zenodo key.")
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--context-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--tool-version", required=True)
    parser.add_argument("--pipeline-version", default=None,
                        help="Version written to the manifest; defaults to --tool-version.")
    parser.add_argument("--zenodo-get", default="zenodo_get")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    work_list = json.loads(Path(args.work_list).read_text(encoding="utf-8"))
    if work_list.get("schema_version") != SUPPORTED_SCHEMA_VERSION:
        sys.exit(f"Unsupported work list schema_version: {work_list.get('schema_version')!r}")

    prefetched = {}
    if args.prefetched_dir and Path(args.prefetched_dir).is_dir():
        prefetched = {path.name: path for path in Path(args.prefetched_dir).iterdir()}
    work = Path(args.work_dir)
    dirs = {"download": work / "download", "extract": work / "extract",
            "raw": Path(args.raw_dir), "bundle": Path(args.bundle_dir),
            "context": Path(args.context_dir)}
    for folder in dirs.values():
        folder.mkdir(parents=True, exist_ok=True)

    manifest = {
        "schema_version": 1,
        "run_id": work_list["run_id"],
        "tool_version": args.pipeline_version or args.tool_version,
        "archives": [],
        "files": [],
        # Not read by the NOMAD side yet; listed for the report.
        "download_failures": [],
    }
    for record in work_list.get("records", []):
        archives, files, failures = process_record(
            record, work_list, prefetched, dirs, args.zenodo_get
        )
        manifest["archives"].extend(archives)
        manifest["files"].extend(files)
        manifest["download_failures"].extend(failures)
        logger.info("Record %s: %d archives, %d files, %d download failures",
                    record["record_id"], len(archives), len(files), len(failures))
    shutil.rmtree(work, ignore_errors=True)
    with open(args.manifest, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


if __name__ == "__main__":
    main()
