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


def load_sequence_label_table(path: str) -> tuple[list[str], list[str]]:
    """Read a 2-column tab-separated [sequence, label] table.

    Matches tools/dna_foundation_models' convention exactly: a header row is
    present but skipped positionally (its actual column names are never
    checked), and there is no id column -- the sequence string itself is the
    only identifier carried through to the output.
    """
    import csv

    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))[1:]

    if not rows:
        raise ValueError(f"{path} has no data rows (expected a header row followed by sequence/label rows).")

    malformed = [i + 2 for i, row in enumerate(rows) if len(row) != 2]  # +2: 1-indexed, header already skipped
    if malformed:
        raise ValueError(
            f"{path}: line(s) {malformed[:10]} do not have exactly 2 tab-separated columns "
            "(sequence, label)."
        )

    seqs = [row[0] for row in rows]
    labels = [row[1] for row in rows]
    return seqs, labels


def load_sequence_label_table_or_fail(path: str) -> tuple[list[str], list[str]]:
    """load_sequence_label_table, failing with exit code 7 instead of raising."""
    from errors import fail

    try:
        return load_sequence_label_table(path)
    except ValueError as exc:
        fail(str(exc), exit_code=7)


def load_and_filter_sequence_label_table(path: str, max_length: int, on_long_sequence: str) -> tuple[list[str], list[str]]:
    """load_sequence_label_table_or_fail, then warn/length-filter like a FASTA input.

    Rows have no id, so warnings reference 1-indexed row position instead.
    Shared by predict.py and finetune.py.
    """
    import sys

    from embed import filter_by_length
    from errors import fail

    ensure_file_exists(path, "Dataset")
    all_seqs, all_labels = load_sequence_label_table_or_fail(path)
    row_ids = [str(i + 1) for i in range(len(all_seqs))]

    for warning in warn_nonstandard_residues(row_ids, all_seqs):
        print(f"WARNING: {warning}", file=sys.stderr)
    kept_ids, seqs, length_warnings = filter_by_length(row_ids, all_seqs, max_length, on_long_sequence)
    for warning in length_warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if not seqs:
        fail(f"No sequences remain in {path} after length filtering.", exit_code=5)

    kept_indices = [int(row_id) - 1 for row_id in kept_ids]
    labels = [all_labels[i] for i in kept_indices]
    return seqs, labels


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
