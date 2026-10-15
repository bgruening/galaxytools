"""Extract per-residue (H5) and/or per-protein (TSV) embeddings from a FASTA file.

Exit codes: see errors.py.
"""

import argparse
import sys

import h5py
import numpy as np
import pandas as pd
import torch
from errors import fail, run_main
from fasta_utils import (
    ensure_file_exists,
    parse_fasta_or_fail,
    preprocess_seq,
    sanitize_h5_key,
    warn_nonstandard_residues,
)
from model_utils import load_model, resolve_device, resolve_max_length, resolve_model_source
from tqdm import tqdm

# ESM2 (HF) and ProtBERT both prepend one special token (CLS/BOS) and append
# one (EOS/SEP); ProtT5 only appends a trailing EOS. Matches src/pooling.py's
# mean_pool_cls_eos / mean_pool_eos_only conventions in the main pipeline.
_LEADING_SPECIAL_TOKENS = {"esm": 1, "protbert": 1, "prott5": 0}


def tokenize_batch(tokenizer, family: str, seqs: list[str], max_length: int):
    batch = seqs if family == "esm" else [preprocess_seq(s) for s in seqs]
    return tokenizer(
        batch, return_tensors="pt", padding=True, truncation=True, max_length=max_length
    )


def strip_special_tokens(hidden: torch.Tensor, attention_mask: torch.Tensor, family: str) -> list[np.ndarray]:
    """Return one [seq_len, hidden_size] array per sequence, special tokens and padding removed."""
    leading = _LEADING_SPECIAL_TOKENS[family]
    lengths = attention_mask.sum(dim=1)
    out = []
    for i in range(hidden.size(0)):
        real_len = int(lengths[i].item())
        out.append(hidden[i, leading:real_len - 1].cpu().numpy())
    return out


def filter_by_length(ids, seqs, max_length: int, on_long: str):
    """Apply --on_long_sequence to sequences over max_length - 2 (special-token budget).

    Returns (kept_ids, kept_seqs, warnings).
    """
    budget = max_length - 2
    kept_ids, kept_seqs, warnings = [], [], []
    for seq_id, seq in zip(ids, seqs):
        if len(seq) <= budget:
            kept_ids.append(seq_id)
            kept_seqs.append(seq)
            continue
        if on_long == "skip":
            warnings.append(f"Skipped '{seq_id}': length {len(seq)} exceeds max_length budget {budget}.")
        else:
            kept_ids.append(seq_id)
            kept_seqs.append(seq[:budget])
            warnings.append(f"Truncated '{seq_id}' from length {len(seq)} to {budget}.")
    return kept_ids, kept_seqs, warnings


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract protein language model embeddings")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--model_path", default=None, help="Preinstalled model folder (from the huggingface data table)")
    parser.add_argument("--model_archive", default=None, help="Uploaded model archive (tar/tgz)")
    parser.add_argument("--output_format", choices=["per_residue", "per_protein", "both"], default="per_protein")
    parser.add_argument("--output_h5", default=None)
    parser.add_argument("--output_tsv", default=None)
    parser.add_argument("--max_length", type=int, default=1024)
    parser.add_argument("--on_long_sequence", choices=["skip", "truncate"], default="truncate")
    parser.add_argument("--batch_size", type=int, default=8)
    args = parser.parse_args()
    args.model_dir = resolve_model_source(args.model_path, args.model_archive)

    ensure_file_exists(args.input, "Input FASTA file")
    ids, seqs = parse_fasta_or_fail(args.input)

    for warning in warn_nonstandard_residues(ids, seqs):
        print(f"WARNING: {warning}", file=sys.stderr)

    family, tokenizer, model = load_model(args.model_dir)
    device = resolve_device()
    model.eval().to(device)
    print(f"Loaded {family} model from {args.model_dir} on {device}.", file=sys.stderr)
    args.max_length = resolve_max_length(model, args.max_length)

    ids, seqs, length_warnings = filter_by_length(ids, seqs, args.max_length, args.on_long_sequence)
    for warning in length_warnings:
        print(f"WARNING: {warning}", file=sys.stderr)

    if not seqs:
        fail("No sequences remain after length filtering; nothing to embed.", exit_code=5)

    want_h5 = args.output_format in ("per_residue", "both")
    want_tsv = args.output_format in ("per_protein", "both")

    h5_file = h5py.File(args.output_h5, "w") if want_h5 else None
    mean_pooled_rows = [] if want_tsv else None

    try:
        for i in tqdm(range(0, len(seqs), args.batch_size), desc="Embedding"):
            batch_ids = ids[i:i + args.batch_size]
            batch_seqs = seqs[i:i + args.batch_size]
            encoded = tokenize_batch(tokenizer, family, batch_seqs, args.max_length)
            encoded = {k: v.to(device) for k, v in encoded.items()}

            with torch.no_grad():
                out = model(**encoded)

            per_seq = strip_special_tokens(out.last_hidden_state, encoded["attention_mask"], family)
            for seq_id, residues in zip(batch_ids, per_seq):
                if h5_file is not None:
                    key = sanitize_h5_key(seq_id)
                    dataset = h5_file.create_dataset(key, data=residues.astype(np.float32))
                    dataset.attrs["original_id"] = seq_id
                if mean_pooled_rows is not None:
                    mean_pooled_rows.append((seq_id, len(residues), residues.mean(axis=0)))
    finally:
        if h5_file is not None:
            h5_file.close()

    if mean_pooled_rows is not None:
        emb_cols = [f"emb_{j}" for j in range(mean_pooled_rows[0][2].shape[0])]
        df = pd.DataFrame([row[2] for row in mean_pooled_rows], columns=emb_cols)
        df.insert(0, "id", [row[0] for row in mean_pooled_rows])
        df.insert(1, "seq_len", [row[1] for row in mean_pooled_rows])
        df.to_csv(args.output_tsv, sep="\t", index=False)

    print(f"Embedded {len(seqs)} sequences.", file=sys.stderr)


if __name__ == "__main__":
    run_main(main)
