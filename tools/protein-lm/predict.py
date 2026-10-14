"""Train a classifier/regressor on protein language model embeddings.

Follows the PRIMARY reference's (tabicl) generic single-target layout: one
feature matrix (here, mean-pooled per-protein embeddings) plus one target
column, evaluated via an internal train/test split.

Exit codes: see errors.py.
"""

import argparse
import sys

import numpy as np
import pandas as pd
import torch
from embed import filter_by_length, strip_special_tokens, tokenize_batch
from errors import fail, run_main
from fasta_utils import (
    check_id_overlap,
    ensure_file_exists,
    load_labels_or_fail,
    parse_fasta_or_fail,
    warn_nonstandard_residues,
)
from model_utils import load_model, resolve_device, resolve_max_length, resolve_model_source
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_recall_curve,
    r2_score,
    root_mean_squared_error,
)
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.preprocessing import label_binarize
from sklearn.svm import SVC, SVR

SEED = 42

CLASSIFIERS = {
    "lr": lambda args: LogisticRegression(max_iter=1000, random_state=SEED),
    "rf": lambda args: RandomForestClassifier(n_estimators=args.n_estimators, random_state=SEED),
    "svm": lambda args: SVC(kernel=args.svm_kernel, probability=True, random_state=SEED),
    "mlp": lambda args: MLPClassifier(hidden_layer_sizes=tuple(args.hidden_sizes), random_state=SEED),
}
REGRESSORS = {
    "linear": lambda args: LinearRegression(),
    "rf": lambda args: RandomForestRegressor(n_estimators=args.n_estimators, random_state=SEED),
    "svm": lambda args: SVR(kernel=args.svm_kernel),
    "mlp": lambda args: MLPRegressor(hidden_layer_sizes=tuple(args.hidden_sizes), random_state=SEED),
}


def embed_sequences(seqs: list[str], family, tokenizer, model, max_length: int, batch_size: int, device) -> np.ndarray:
    """Embed an already-decided list of sequences (parsing/filtering done by the caller)."""
    rows = []
    for i in range(0, len(seqs), batch_size):
        batch_seqs = seqs[i:i + batch_size]
        encoded = tokenize_batch(tokenizer, family, batch_seqs, max_length)
        encoded = {k: v.to(device) for k, v in encoded.items()}
        with torch.no_grad():
            out = model(**encoded)
        for residues in strip_special_tokens(out.last_hidden_state, encoded["attention_mask"], family):
            rows.append(residues.mean(axis=0))
    return np.vstack(rows)


def build_xy(fasta_path, labels_path, label_col, model_dir, max_length, on_long_sequence, batch_size, device, task):
    """Parse+validate the FULL fasta/labels pair (so a genuine ID mismatch is never
    confused with a sequence that --on_long_sequence skip legitimately dropped),
    then length-filter and embed only the matched survivors.
    """
    ensure_file_exists(fasta_path, "FASTA file")
    ensure_file_exists(labels_path, "Labels table")
    all_ids, all_seqs = parse_fasta_or_fail(fasta_path)
    labels = load_labels_or_fail(labels_path, id_col="id", label_cols=[label_col])

    fasta_only, labels_only = check_id_overlap(all_ids, labels.keys())
    if fasta_only or labels_only:
        fail(
            "Sequence IDs in the FASTA file and labels table do not match. "
            f"In FASTA but not labels: {fasta_only[:10]}{'...' if len(fasta_only) > 10 else ''}. "
            f"In labels but not FASTA: {labels_only[:10]}{'...' if len(labels_only) > 10 else ''}.",
            exit_code=6,
        )

    family, tokenizer, model = load_model(model_dir)
    model.eval().to(device)
    max_length = resolve_max_length(model, max_length)

    for warning in warn_nonstandard_residues(all_ids, all_seqs):
        print(f"WARNING: {warning}", file=sys.stderr)
    ids, seqs, length_warnings = filter_by_length(all_ids, all_seqs, max_length, on_long_sequence)
    for warning in length_warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if not seqs:
        fail(f"No sequences remain in {fasta_path} after length filtering.", exit_code=5)

    raw_values = [labels[seq_id][0] for seq_id in ids]
    if task == "regression":
        try:
            y = np.array([float(v) for v in raw_values], dtype=np.float64)
        except ValueError as exc:
            fail(f"Labels table column '{label_col}' contains a non-numeric value for regression: {exc}", exit_code=7)
    else:
        y = np.array(raw_values)

    embeddings = embed_sequences(seqs, family, tokenizer, model, max_length, batch_size, device)
    return ids, embeddings, y


def performance_metrics(y_true, y_pred, task: str) -> dict:
    if task == "classification":
        return {
            "accuracy": accuracy_score(y_true, y_pred),
            "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
            "f1_weighted": f1_score(y_true, y_pred, average="weighted", zero_division=0),
        }
    return {
        "rmse": root_mean_squared_error(y_true, y_pred),
        "r2": r2_score(y_true, y_pred),
    }


def prediction_plot(y_true, y_pred, task, y_scores=None):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6))
    if task == "classification":
        classes = np.unique(y_true)
        if len(classes) == 2:
            y_binary = label_binarize(y_true, classes=classes).ravel()
            precision, recall, _ = precision_recall_curve(y_binary, y_scores[:, 1])
            ap = average_precision_score(y_binary, y_scores[:, 1])
            plt.plot(recall, precision, label=f"Precision-recall (AP={ap:.2f})")
            plt.title("Precision-recall curve (binary classification)")
        else:
            y_binarized = label_binarize(y_true, classes=classes)
            for index, class_name in enumerate(classes):
                precision, recall, _ = precision_recall_curve(y_binarized[:, index], y_scores[:, index])
                ap = average_precision_score(y_binarized[:, index], y_scores[:, index])
                plt.plot(recall, precision, label=f"{class_name} (AP={ap:.2f})")
            plt.title("Precision-recall curve (multiclass classification)")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.legend(loc="lower left")
    else:
        rmse, r2 = root_mean_squared_error(y_true, y_pred), r2_score(y_true, y_pred)
        plt.scatter(y_true, y_pred, alpha=0.8)
        lo, hi = min(np.min(y_true), np.min(y_pred)), max(np.max(y_true), np.max(y_pred))
        plt.plot([lo, hi], [lo, hi], "r--", label="y = x")
        plt.xlabel("True values")
        plt.ylabel("Predicted values")
        plt.title(f"True vs predicted (RMSE={rmse:.2f}, R2={r2:.2f})")
        plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("output_plot.png")
    plt.close()


def make_estimator(args):
    table = CLASSIFIERS if args.task == "classification" else REGRESSORS
    return table[args.estimator](args)


def train_test(args, device):
    """Embed the whole dataset once, then split it into train/test by --test_fraction."""
    ids, X, y = build_xy(
        args.data_fasta, args.data_labels, args.label_col, args.model_dir,
        args.max_length, args.on_long_sequence, args.batch_size, device, args.task,
    )

    stratify = y if args.task == "classification" else None
    if args.task == "classification":
        _, counts = np.unique(y, return_counts=True)
        if counts.min() < 2:
            fail(
                "Cannot split into train/test: every class must have at least 2 samples "
                f"(smallest class has {counts.min()}). Use more data, or a larger test "
                "fraction so every class can appear on both sides of the split.",
                exit_code=7,
            )

    indices = np.arange(len(ids))
    try:
        train_idx, test_idx = train_test_split(
            indices, test_size=args.test_fraction, random_state=args.random_seed, stratify=stratify
        )
    except ValueError as exc:
        fail(f"Cannot split into train/test with test_fraction={args.test_fraction}: {exc}", exit_code=7)

    estimator = make_estimator(args)
    estimator.fit(X[train_idx], y[train_idx])
    predicted = estimator.predict(X[test_idx])
    y_test = y[test_idx]
    test_ids = [ids[i] for i in test_idx]

    pd.DataFrame([performance_metrics(y_test, predicted, args.task)]).to_csv(
        "test_metrics.tsv", sep="\t", index=False
    )
    scores = estimator.predict_proba(X[test_idx]) if args.task == "classification" else None
    prediction_plot(y_test, predicted, args.task, scores)

    output = pd.DataFrame({"id": test_ids, "true_label": y_test, "predicted_label": predicted})
    output.to_csv("output_predicted_data", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a classifier/regressor on protein LM embeddings")
    parser.add_argument("--model_path", default=None, help="Preinstalled model folder (from the huggingface data table)")
    parser.add_argument("--model_archive", default=None, help="Uploaded model archive (tar/tgz)")
    parser.add_argument("--task", choices=["classification", "regression"], required=True)
    parser.add_argument("--data_fasta", required=True)
    parser.add_argument("--data_labels", required=True)
    parser.add_argument("--label_col", required=True)
    parser.add_argument("--test_fraction", type=float, default=0.2)
    parser.add_argument("--random_seed", type=int, default=SEED)
    parser.add_argument("--estimator", choices=["lr", "rf", "svm", "mlp", "linear"], default="rf")
    parser.add_argument("--n_estimators", type=int, default=100)
    parser.add_argument("--svm_kernel", default="rbf")
    parser.add_argument("--hidden_sizes", type=int, nargs="+", default=[128])
    parser.add_argument("--max_length", type=int, default=1024)
    parser.add_argument("--on_long_sequence", choices=["skip", "truncate"], default="truncate")
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = parser.parse_args()
    args.model_dir = resolve_model_source(args.model_path, args.model_archive)

    device = resolve_device(args.device)
    train_test(args, device)


if __name__ == "__main__":
    run_main(main)
