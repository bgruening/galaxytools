"""Merge sylph-tax reports while retaining samples with no detected genomes."""

import argparse
import csv
import subprocess
from pathlib import Path


def merge_taxonomy(profile, sketches, taxonomy, output, column, paired=False):
    suffix = ".paired.sylsp" if paired else ".sylsp"
    samples = sorted(path.name[:-len(suffix)] for path in Path(sketches).glob(f"*{suffix}"))
    if not samples:
        raise ValueError("No read sketches found; cannot identify samples")

    with open(profile) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"Sample_file", "Genome_file"} <= set(reader.fieldnames):
            raise ValueError("Missing Sylph profile header; this is not a valid zero-hit result")
        positive_samples = {row["Sample_file"] for row in reader}

    files = sorted(Path(taxonomy).glob("*.sylphmpa"))
    if files:
        subprocess.run(
            ["sylph-tax", "merge", *map(str, files), "--column", column, "-o", str(output)],
            check=True,
        )
        with open(output) as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        header = rows[0]
        if header[0] != "clade_name" or set(header[1:]) != positive_samples:
            raise ValueError("Taxonomy report sample columns do not match the Sylph profile")
        if not positive_samples <= set(samples):
            raise ValueError("Sylph profile samples do not match the read sketches")
        if set(header[1:]) == set(samples):
            return
    else:
        if positive_samples:
            raise ValueError("No taxonomy reports were produced for a nonempty Sylph profile")
        header, rows = ["clade_name"], [["clade_name"]]

    # Preserve the upstream values verbatim, adding zeros only for absent samples.
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["clade_name", *samples])
        for row in rows[1:]:
            values = dict(zip(header[1:], row[1:]))
            writer.writerow([row[0], *(values.get(sample, "0") for sample in samples)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--sketches", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--column", required=True)
    parser.add_argument("--paired", action="store_true")
    merge_taxonomy(**vars(parser.parse_args()))
