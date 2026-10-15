"""Train an MLP classifier/regressor on protein language model embeddings.

Input convention matches tools/dna_foundation_models: a 2-column tab-separated
table (sequence, label), separate pre-split train/test files -- no FASTA, no
internal splitting. Users who need a split can use the existing
sklearn_train_test_split Galaxy tool first.

Exit codes: see errors.py.
"""

import argparse

import numpy as np
import pandas as pd
import torch
from embed import strip_special_tokens, tokenize_batch
from errors import fail, run_main
from fasta_utils import load_and_filter_sequence_label_table
from model_utils import load_model, resolve_device, resolve_max_length, resolve_model_source
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_recall_curve,
    r2_score,
    root_mean_squared_error,
)
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.preprocessing import label_binarize


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


def build_xy(data_path, family, tokenizer, model, max_length, on_long_sequence, batch_size, device, task):
    """Load a [sequence, label] tabular file, length-filter, and embed it."""
    seqs, raw_values = load_and_filter_sequence_label_table(data_path, max_length, on_long_sequence)
    if task == "regression":
        try:
            y = np.array([float(v) for v in raw_values], dtype=np.float64)
        except ValueError as exc:
            fail(f"{data_path}: label column contains a non-numeric value for regression: {exc}", exit_code=7)
    else:
        y = np.array(raw_values)

    embeddings = embed_sequences(seqs, family, tokenizer, model, max_length, batch_size, device)
    return seqs, embeddings, y


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
    if args.task == "classification":
        return MLPClassifier(hidden_layer_sizes=tuple(args.hidden_sizes), random_state=42)
    return MLPRegressor(hidden_layer_sizes=tuple(args.hidden_sizes), random_state=42)


def train_test(args, device):
    family, tokenizer, model = load_model(args.model_dir)
    model.eval().to(device)
    max_length = resolve_max_length(model, args.max_length)

    train_seqs, X_train, y_train = build_xy(
        args.train_data, family, tokenizer, model, max_length, args.on_long_sequence, args.batch_size, device, args.task,
    )
    test_seqs, X_test, y_test = build_xy(
        args.test_data, family, tokenizer, model, max_length, args.on_long_sequence, args.batch_size, device, args.task,
    )

    estimator = make_estimator(args)
    estimator.fit(X_train, y_train)
    predicted = estimator.predict(X_test)

    pd.DataFrame([performance_metrics(y_test, predicted, args.task)]).to_csv(
        "test_metrics.tsv", sep="\t", index=False
    )
    scores = estimator.predict_proba(X_test) if args.task == "classification" else None
    prediction_plot(y_test, predicted, args.task, scores)

    output = pd.DataFrame({"sequence": test_seqs, "label": y_test, "prediction": predicted})
    output.to_csv("output_predicted_data", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train an MLP classifier/regressor on protein LM embeddings")
    parser.add_argument("--model_path", default=None, help="Preinstalled model folder (from the huggingface data table)")
    parser.add_argument("--model_archive", default=None, help="Uploaded model archive (tar/tgz)")
    parser.add_argument("--task", choices=["classification", "regression"], required=True)
    parser.add_argument("--train_data", required=True, help="Tabular [sequence, label] training file")
    parser.add_argument("--test_data", required=True, help="Tabular [sequence, label] test file")
    parser.add_argument("--hidden_sizes", type=int, nargs="+", default=[128])
    parser.add_argument("--max_length", type=int, default=1024)
    parser.add_argument("--on_long_sequence", choices=["skip", "truncate"], default="truncate")
    parser.add_argument("--batch_size", type=int, default=8)
    args = parser.parse_args()
    args.model_dir = resolve_model_source(args.model_path, args.model_archive)

    device = resolve_device()
    train_test(args, device)


if __name__ == "__main__":
    run_main(main)
