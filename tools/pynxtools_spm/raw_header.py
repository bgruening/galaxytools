#!/usr/bin/env python
"""Write the header (metadata) of a raw file as JSON, for the ELN tools.

The header is read with gwyddionpy (https://github.com/FAIRmat-NFDI/gwyddionPy), which reads many
scanning probe microscopy formats, also binary ones. A file it cannot read gives its leading text.
"""

import argparse
import json
import logging
import sys

import gwyddionpy
from gwyddionpy.export.json import to_metadata_dict

logger = logging.getLogger("raw_header")

MAX_TEXT_BYTES = 64 * 1024
# Lines that start the data block of a text format (Nanonis .dat); the header ends there.
DATA_MARKERS = ("[DATA]",)


def gwyddion_header(path):
    """(gwyddionpy's metadata JSON payload, format), or (None, reason)."""
    try:
        data = gwyddionpy.load(path)
    except (gwyddionpy.GwyddionPyError, OSError) as error:
        return None, f"gwyddionpy: {type(error).__name__}"
    if not data.metadata:
        return None, f"gwyddionpy: no metadata ({data.source_format})"
    return shared_meta(to_metadata_dict(data)), data.source_format


def shared_meta(payload):
    """Move the metadata that all channels share to a top-level `meta`, so each fact is written once."""
    channels = list(payload["channels"].values())
    metas = [channel.get("meta") or {} for channel in channels]
    shared = {key: value for key, value in metas[0].items()
              if all(key in meta and meta[key] == value for meta in metas[1:])} if metas else {}
    for channel in channels:
        own = {key: value for key, value in (channel.get("meta") or {}).items() if key not in shared}
        channel.pop("meta", None)
        if own:
            channel["meta"] = own
    return {"source_format": payload["source_format"], "meta": shared, "channels": payload["channels"]}


def leading_text(path):
    """The text before the first NUL byte or data marker, or None if it is not text."""
    with open(path, "rb") as handle:
        data = handle.read(MAX_TEXT_BYTES).split(b"\0", 1)[0]
    lines = []
    for line in data.decode("latin-1").splitlines():
        if line.strip() in DATA_MARKERS:
            break
        lines.append(line.rstrip())
    text = "\n".join(lines).strip()
    printable = sum(char.isprintable() or char.isspace() for char in text)
    return text if text and printable >= 0.95 * len(text) else None


def raw_header(path, name):
    """The header JSON: `file`, `read_with`, then gwyddionpy's metadata or the leading `text`."""
    header = {"file": name}
    payload, source = gwyddion_header(path)
    if payload:
        header["read_with"] = f"gwyddionpy ({source})"
        return header | payload
    text = leading_text(path)
    header["read_with"] = f"leading text; {source}" if text else f"no header; {source}"
    header["text"] = text or ""
    return header


def main():
    parser = argparse.ArgumentParser(description="Write the header of a raw file as JSON.")
    parser.add_argument("raw", help="The raw file.")
    parser.add_argument("--name", required=True, help="The raw file's name, written as `file`.")
    parser.add_argument("--header", required=True, help="Output JSON file.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr)])

    header = raw_header(args.raw, args.name)
    # Vendor metadata holds µ, ° and Å; non-finite numbers are already strings (gwyddionpy).
    with open(args.header, "w", encoding="utf-8") as handle:
        json.dump(header, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
    logger.info("%s: header %s", args.name, header["read_with"])


if __name__ == "__main__":
    main()
