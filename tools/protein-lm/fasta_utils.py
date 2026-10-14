"""Plain-Python FASTA parsing and FASTA/labels validation (no biopython dependency)."""

from pathlib import Path

_RARE_AA_TRANS = str.maketrans("UZOB", "XXXX")
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(fasta_path: str) -> tuple[list[str], list[str]]:
    """Read a FASTA file and return (ids, sequences).

    ``id`` is the first whitespace-separated token after '>', matching
    Biopython's ``record.id`` behaviour. Sequence lines are concatenated
    and whitespace-stripped until the next '>' or end of file.
    """
    ids: list[str] = []
    seqs: list[str] = []
    current_seq: list[str] = []
    with open(fasta_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if ids:
                    seqs.append("".join(current_seq))
                ids.append(line[1:].split()[0] if line[1:].split() else "")
                current_seq = []
            else:
                current_seq.append(line.strip())
        if ids:
            seqs.append("".join(current_seq))

    return ids, seqs


def preprocess_seq(seq: str) -> str:
    """Replace rare residues (U/Z/O/B -> X) and space-separate residues.

    Required by ProtBERT/ProtT5 tokenisers. Do NOT apply to ESM2 inputs,
    which tokenise raw sequences directly.
    """
    return " ".join(seq.upper().translate(_RARE_AA_TRANS))


def warn_nonstandard_residues(ids: list[str], seqs: list[str]) -> list[str]:
    """Return human-readable warnings for sequences containing non-standard residues."""
    warnings = []
    for seq_id, seq in zip(ids, seqs):
        nonstandard = sorted(set(seq.upper()) - _STANDARD_AA)
        if nonstandard:
            warnings.append(
                f"Sequence '{seq_id}' contains non-standard residue(s) {nonstandard}; "
                "these will be mapped to X by ProtBERT/ProtT5, or passed through as-is to ESM2."
            )
    return warnings


def load_labels(labels_path: str, id_col: str, label_cols: list[str]) -> dict[str, list[str]]:
    """Read a tab-separated labels table into {id: [label_values...]}."""
    import csv

    labels: dict[str, list[str]] = {}
    with open(labels_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in [id_col, *label_cols] if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                f"Labels table {labels_path} is missing column(s): {missing}. "
                f"Columns found: {reader.fieldnames}. "
                "If this looks like one column was split or merged unexpectedly, check that the "
                "file is tab-separated (not comma-separated) and that the header matches exactly "
                "(case-sensitive, no extra whitespace)."
            )
        for row in reader:
            labels[row[id_col]] = [row[c] for c in label_cols]
    return labels


def load_labels_or_fail(labels_path: str, id_col: str, label_cols: list[str]) -> dict[str, list[str]]:
    """load_labels, failing with exit code 7 (invalid label values) instead of raising."""
    from errors import fail

    try:
        return load_labels(labels_path, id_col=id_col, label_cols=label_cols)
    except ValueError as exc:
        fail(str(exc), exit_code=7)


def check_id_overlap(fasta_ids: list[str], label_ids) -> tuple[list[str], list[str]]:
    """Return (fasta_only, labels_only): ids present in one input but not the other."""
    fasta_set = set(fasta_ids)
    label_set = set(label_ids)
    return sorted(fasta_set - label_set), sorted(label_set - fasta_set)


def sanitize_h5_key(seq_id: str) -> str:
    """Sanitise a sequence id for use as an HDF5 dataset key.

    HDF5 keys may contain '/' (interpreted as a group separator) and other
    characters that are awkward in practice; replace anything outside a safe
    set with '_'. The original id is preserved separately as an attribute,
    so this only affects the key, never the recorded identity of the record.
    """
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in seq_id)


def ensure_file_exists(path: str, description: str) -> None:
    if not path or not Path(path).exists():
        raise FileNotFoundError(f"{description} does not exist: {path}")


def parse_fasta_or_fail(fasta_path: str) -> tuple[list[str], list[str]]:
    """parse_fasta, failing with exit code 5 if the file contains zero sequences."""
    from errors import fail

    ids, seqs = parse_fasta(fasta_path)
    if not ids:
        fail(f"No sequences found in FASTA file: {fasta_path}", exit_code=5)
    return ids, seqs
