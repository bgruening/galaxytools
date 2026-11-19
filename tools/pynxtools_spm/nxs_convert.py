#!/usr/bin/env python
"""Convert one raw file to NeXus with ``pynx convert`` and always write a status.

A failed conversion is recorded in ``status.json`` (``status: failed`` plus the
reason) and the script still exits with 0. When the tool is mapped over a
collection, one bad file then does not stop the other files, and the report step
can list every failure. Only a broken tool (for example an unwritable output)
exits with an error.

``status.json`` carries the fields of ``FileResult`` in nomad-eosc-galaxy-actions
(``sha256``, ``status``, ``error``), so results can be joined with the manifest:
https://github.com/FAIRmat-NFDI/nomad-eosc-galaxy-actions
"""

import argparse
import hashlib
import json
import logging
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

logger = logging.getLogger("nxs_convert")

# Longest error message kept in status.json; the full output is in the log.
MAX_ERROR_LENGTH = 500


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(2**20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_name(name):
    """A plain file name that keeps the extension the reader looks at."""
    name = re.sub(r"[^\w.\-]", "_", Path(name).name).lstrip(".")
    return name or "raw_file"


def last_error(log_text):
    """The most telling line of the converter output: the last one naming an error."""
    lines = [line.strip() for line in log_text.splitlines() if line.strip()]
    for line in reversed(lines):
        if re.search(r"error|exception", line, re.IGNORECASE):
            return line[:MAX_ERROR_LENGTH]
    return lines[-1][:MAX_ERROR_LENGTH] if lines else "converter produced no output"


def nxdl_from_eln(eln):
    """The ELN's `definition` (for example NXstm), or None."""
    try:
        content = yaml.safe_load(Path(eln).read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, yaml.YAMLError):
        return None
    definition = content.get("definition") if isinstance(content, dict) else None
    return definition if isinstance(definition, str) and definition.strip() else None


def convert(raw, raw_name, eln, config, reader, nxdl, output_nxs, log_path):
    """Run ``pynx convert``; return the error message, or None on success."""
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        # Galaxy names every dataset `dataset_N.dat`, but the reader picks inputs
        # by extension (`.dat` even counts as an STS raw file), so link real names.
        inputs = [(raw, safe_name(raw_name)), (eln, "eln_data.yaml")]
        if config:
            # Passed as an input file: pynxtools-spm 0.4.1 and 0.4.2 do not accept `-c`.
            inputs.append((config, "config.json"))
        for source, name in inputs:
            (work / name).symlink_to(Path(source).resolve())
        output = work / "output.nxs"
        command = [
            "pynx", "convert", "--reader", reader, "--nxdl", nxdl,
            "--output", str(output), *(name for _, name in inputs),
        ]
        logger.info("Running: %s", " ".join(command))
        completed = subprocess.run(
            command, cwd=work, capture_output=True, text=True, check=False
        )
        log_text = completed.stdout + completed.stderr
        Path(log_path).write_text(log_text, encoding="utf-8")
        if completed.returncode != 0:
            return last_error(log_text)
        if not output.is_file():
            return "converter finished but wrote no NeXus file"
        Path(output_nxs).write_bytes(output.read_bytes())
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Convert one raw file to NeXus and write status.json."
    )
    parser.add_argument("raw", help="Path to the raw file.")
    parser.add_argument(
        "--raw-name", required=True, help="Original file name, for its extension."
    )
    parser.add_argument("--eln", required=True, help="ELN YAML file.")
    parser.add_argument("--config", default=None, help="Optional reader config JSON.")
    parser.add_argument("--reader", default="spm")
    parser.add_argument(
        "--nxdl", required=True, help="Application definition, e.g. NXstm, or `auto` for the ELN's."
    )
    parser.add_argument("--output-nxs", required=True)
    parser.add_argument("--status-json", required=True)
    parser.add_argument("--log", required=True, help="Converter output log.")
    parser.add_argument("--tool-version", required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    # Galaxy expects every output to exist; after a failed conversion it stays empty.
    Path(args.output_nxs).write_bytes(b"")
    nxdl = nxdl_from_eln(args.eln) if args.nxdl == "auto" else args.nxdl
    if nxdl is None:
        error = "the ELN has no `definition`; choose an application definition"
        Path(args.log).write_text(error + "\n", encoding="utf-8")
    else:
        error = convert(
            args.raw, args.raw_name, args.eln, args.config,
            args.reader, nxdl, args.output_nxs, args.log,
        )
    status = {
        "sha256": sha256_of(args.raw),
        "status": "failed" if error else "success",
        "error": error,
        "raw_name": args.raw_name,
        "reader": args.reader,
        "nxdl": nxdl,
        "tool_version": args.tool_version,
    }
    with open(args.status_json, "w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2)
    logger.info("Conversion of %s: %s", args.raw_name, status["status"])


if __name__ == "__main__":
    main()
