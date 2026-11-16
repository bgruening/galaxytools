#!/usr/bin/env python
"""Write the header of a raw file as text, for the ELN tools.

The header is read with gwyddionpy (https://github.com/FAIRmat-NFDI/gwyddionPy), which reads many
scanning probe microscopy formats, also binary ones. A file it cannot read gives its leading text.
"""

import argparse
import logging
import sys

import gwyddionpy

logger = logging.getLogger("raw_header")

MAX_TEXT_BYTES = 64 * 1024
# Lines that start the data block of a text format (Nanonis .dat); the header ends there.
DATA_MARKERS = ("[DATA]",)


def gwyddion_header(path):
    """(header lines, format) from gwyddionpy, or (None, reason)."""
    try:
        data = gwyddionpy.load(path)
    except (gwyddionpy.GwyddionPyError, OSError) as error:
        return None, f"gwyddionpy: {type(error).__name__}"
    lines = [f"{key}: {value}" for key, value in data.metadata.items()]
    if not lines:
        return None, f"gwyddionpy: no metadata ({data.source_format})"
    return lines, data.source_format


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
    return lines if text and printable >= 0.95 * len(text) else None


def raw_header(path):
    """(header lines, how it was read); empty lines if nothing was found."""
    lines, source = gwyddion_header(path)
    if lines:
        return lines, f"gwyddionpy ({source})"
    text = leading_text(path)
    if text:
        return text, f"leading text; {source}"
    return [], f"no header; {source}"


def main():
    parser = argparse.ArgumentParser(description="Write the header of a raw file as text.")
    parser.add_argument("raw", help="The raw file.")
    parser.add_argument("--name", required=True, help="The raw file's name, for the first line.")
    parser.add_argument("--header", required=True, help="Output text file.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr)])

    lines, method = raw_header(args.raw)
    with open(args.header, "w", encoding="utf-8") as handle:
        handle.write(f"# Header of {args.name}, read with: {method}\n")
        handle.writelines(f"{line}\n" for line in lines)
    logger.info("%s: %d header lines, %s", args.name, len(lines), method)


if __name__ == "__main__":
    main()
