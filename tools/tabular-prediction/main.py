"""TabICL in-context classification, regression and SHAP runner."""

import argparse
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_recall_curve,
    r2_score,
    root_mean_squared_error,
)
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import label_binarize

SEED = 42


def read_table(path, header):
    return pd.read_csv(path, sep="\t", header=0 if header == "true" else None)


def split_xy(path, header):
    data = read_table(path, header)
    return data.iloc[:, :-1].copy(), data.iloc[:, -1].copy()


def parse_auto_boolean(value):
    if value == "auto":
        return "auto"
    return value == "true"


def optional_int(value):
    if value.lower() in {"", "none"}:
        return None
    return int(value)


def model_config(args):
    config = {"random_state": args.random_state, "model_path": args.model_path}
    if args.advanced_icl != "true":
        return config
    norm_methods = args.norm_methods.split(",")
    if norm_methods == ["none_power"]:
        norm_methods = ["none", "power"]
    allowed_norm_methods = {"none", "power", "quantile", "quantile_rtdl", "robust"}
    if not set(norm_methods) <= allowed_norm_methods:
        raise ValueError("Unknown normalization method.")
    config.update(
        {
            "n_estimators": args.n_estimators,
            "norm_methods": norm_methods,
            "feat_shuffle_method": args.feat_shuffle_method,
            "outlier_threshold": args.outlier_threshold,
            "batch_size": args.batch_size,
            "kv_cache": {"false": False, "true": True, "kv": "kv", "repr": "repr"}[
                args.kv_cache
            ],
            "device": None,
            "use_amp": "auto",
            "use_fa3": "auto",
            "offload_mode": "auto",
            "disk_offload_dir": None,
            "verbose": args.verbose == "true",
            "inference_config": None,
        }
    )
    if args.selected_task == "Classification":
        config.update(
            {
                "class_shuffle_method": args.class_shuffle_method,
                "softmax_temperature": args.softmax_temperature,
                "average_logits": args.average_logits == "true",
                "support_many_classes": args.support_many_classes == "true",
            }
        )
    return config


def make_estimator(args):
    from tabicl import TabICLClassifier, TabICLRegressor

    cls = (
        TabICLClassifier if args.selected_task == "Classification" else TabICLRegressor
    )
    return cls(**model_config(args))


def prediction_plot(y_true, y_pred, task, y_scores=None):
    plt.figure(figsize=(8, 6))
    if task == "Classification":
        classes = np.unique(y_true)
        if y_scores is None:
            raise ValueError(
                "Classification plotting requires prediction probabilities."
            )
        if len(classes) == 2:
            y_binary = label_binarize(y_true, classes=classes).ravel()
            precision, recall, _ = precision_recall_curve(y_binary, y_scores[:, 1])
            average_precision = average_precision_score(y_binary, y_scores[:, 1])
            plt.plot(
                recall,
                precision,
                label=f"Precision-recall (AP={average_precision:.2f})",
            )
            plt.title("Precision-recall curve (binary classification)")
        else:
            y_binarized = label_binarize(y_true, classes=classes)
            for index, class_name in enumerate(classes):
                precision, recall, _ = precision_recall_curve(
                    y_binarized[:, index], y_scores[:, index]
                )
                average_precision = average_precision_score(
                    y_binarized[:, index], y_scores[:, index]
                )
                plt.plot(
                    recall,
                    precision,
                    label=f"{class_name} (AP={average_precision:.2f})",
                )
            precision, recall, _ = precision_recall_curve(
                y_binarized.ravel(), y_scores.ravel()
            )
            plt.plot(recall, precision, "--", color="black", label="Micro-average")
            plt.title("Precision-recall curve (multiclass classification)")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.legend(loc="lower left")
    else:
        rmse, r2 = root_mean_squared_error(y_true, y_pred), r2_score(y_true, y_pred)
        plt.scatter(y_true, y_pred, alpha=0.8)
        low, high = min(np.min(y_true), np.min(y_pred)), max(
            np.max(y_true), np.max(y_pred)
        )
        plt.plot([low, high], [low, high], "r--", label="y = x")
        plt.xlabel("True values")
        plt.ylabel("Predicted values")
        plt.title(f"True vs predicted (RMSE={rmse:.2f}, R2={r2:.2f})")
        plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("output_plot.png")
    plt.close()


def make_shap_plot(estimator, features, limit, plot_type="beeswarm"):
    from tabicl.shap import get_shap_explainer, plot_shap

    explain_data = features.iloc[:limit]
    explain_array = np.asarray(explain_data, dtype=np.float64)
    predict_method = (
        "predict_proba" if hasattr(estimator, "predict_proba") else "predict"
    )
    explainer = get_shap_explainer(
        estimator, explain_array, predict_fn=predict_method, algorithm="permutation"
    )
    # Permutation SHAP requires at least two evaluations per feature plus one.
    values = explainer(explain_array, max_evals=2 * explain_array.shape[1] + 1)
    if hasattr(values, "feature_names"):
        values.feature_names = [str(column) for column in features.columns]
    plot_shap(values, kind=plot_type)
    plt.tight_layout()
    plt.savefig("shap_plot.png", bbox_inches="tight")
    plt.close()


def performance_metrics(y_true, y_pred, task):
    if task == "Classification":
        return {
            "accuracy": accuracy_score(y_true, y_pred),
            "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
            "f1_weighted": f1_score(
                y_true, y_pred, average="weighted", zero_division=0
            ),
        }
    return {
        "rmse": root_mean_squared_error(y_true, y_pred),
        "r2": r2_score(y_true, y_pred),
    }


def train_test(args):
    x_train, y_train = split_xy(args.train_data, args.train_header)
    if args.testhaslabels == "true":
        x_test, y_test = split_xy(args.test_data, args.test_header)
    else:
        x_test, y_test = read_table(args.test_data, args.test_header), None
    estimator = make_estimator(args)
    estimator.fit(x_train, y_train)
    predicted = estimator.predict(x_test)
    if y_test is not None:
        pd.DataFrame(
            [performance_metrics(y_test, predicted, args.selected_task)]
        ).to_csv("test_metrics.tsv", sep="\t", index=False)
        scores = (
            estimator.predict_proba(x_test)
            if args.selected_task == "Classification"
            else None
        )
        prediction_plot(y_test, predicted, args.selected_task, scores)
    if args.shap == "true":
        make_shap_plot(estimator, x_test, args.shap_max_samples, args.shap_plot_type)
    output = x_test.copy()
    if y_test is not None:
        output["true_labels"] = y_test.to_numpy()
    output["predicted_labels"] = predicted
    output.to_csv("output_predicted_data", sep="\t", index=False)


def cross_validate(args):
    features, labels = split_xy(args.train_data, args.train_header)
    if args.selected_task == "Classification" and args.cv_strategy == "stratified":
        too_small = labels.value_counts()[lambda counts: counts < args.n_splits]
        if not too_small.empty:
            raise ValueError(
                "Cannot run stratified cross validation: each class must contain at "
                f"least {args.n_splits} samples. Classes below that limit: "
                + ", ".join(map(str, too_small.index))
            )
        splits = StratifiedKFold(
            args.n_splits, shuffle=True, random_state=args.random_state
        ).split(features, labels)
    else:
        if args.selected_task == "Regression" and args.cv_strategy == "stratified":
            raise ValueError(
                "Stratified cross validation is only available for classification."
            )
        splits = KFold(
            args.n_splits, shuffle=True, random_state=args.random_state
        ).split(features)
    predictions = pd.Series(index=features.index, dtype=object)
    fold_numbers = pd.Series(index=features.index, dtype="Int64")
    metrics = []
    for fold_number, (train_index, test_index) in enumerate(splits, 1):
        estimator = make_estimator(args)
        estimator.fit(features.iloc[train_index], labels.iloc[train_index])
        predicted = estimator.predict(features.iloc[test_index])
        predictions.iloc[test_index], fold_numbers.iloc[test_index] = (
            predicted,
            fold_number,
        )
        fold_metrics = performance_metrics(
            labels.iloc[test_index], predicted, args.selected_task
        )
        metrics.append({"fold": fold_number, **fold_metrics})
        metric_columns = list(fold_metrics)
    output = features.copy()
    output["true_labels"], output["fold"], output["predicted_labels"] = (
        labels,
        fold_numbers,
        predictions,
    )
    output.to_csv("output_predicted_data", sep="\t", index=False)
    metrics_df = pd.DataFrame(metrics)
    summary = [
        {
            "fold": name,
            **{
                column: getattr(metrics_df[column], name)() for column in metric_columns
            },
        }
        for name in ("mean", "std")
    ]
    pd.concat([metrics_df, pd.DataFrame(summary)], ignore_index=True).to_csv(
        "cv_metrics.tsv", sep="\t", index=False
    )


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected_task", required=True)
    parser.add_argument("--evaluation_method", default="train_test")
    parser.add_argument("--train_data", required=True)
    parser.add_argument("--train_header", required=True)
    parser.add_argument("--test_data")
    parser.add_argument("--test_header")
    parser.add_argument("--testhaslabels", default="false")
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--advanced_icl", default="false")
    parser.add_argument("--n_estimators", type=int, default=8)
    parser.add_argument(
        "--norm_methods",
        default="none_power",
        help="Comma-separated normalization methods",
    )
    parser.add_argument(
        "--feat_shuffle_method",
        choices=["none", "shift", "random", "latin"],
        default="latin",
    )
    parser.add_argument(
        "--class_shuffle_method",
        choices=["none", "shift", "random", "latin"],
        default="shift",
    )
    parser.add_argument("--outlier_threshold", type=float, default=4.0)
    parser.add_argument("--softmax_temperature", type=float, default=0.9)
    parser.add_argument("--average_logits", default="true")
    parser.add_argument("--support_many_classes", default="true")
    parser.add_argument("--batch_size", type=optional_int, default=8)
    parser.add_argument(
        "--kv_cache", choices=["false", "true", "kv", "repr"], default="false"
    )
    parser.add_argument("--use_amp", choices=["auto", "true", "false"], default="auto")
    parser.add_argument("--random_state", type=optional_int, default=SEED)
    parser.add_argument("--n_jobs", type=optional_int, default=0)
    parser.add_argument("--verbose", default="false")
    parser.add_argument("--inference_config", default="")
    parser.add_argument("--n_splits", type=int, default=5)
    parser.add_argument("--cv_strategy", default="stratified")
    parser.add_argument("--shap", default="false")
    parser.add_argument("--shap_max_samples", type=int, default=10)
    parser.add_argument(
        "--shap_plot_type", choices=["bar", "scatter", "beeswarm"], default="beeswarm"
    )
    return parser


def main():
    args = make_parser().parse_args()
    start = time.time()
    if args.evaluation_method == "cross_validation":
        cross_validate(args)
    else:
        train_test(args)
    print(f"Time taken: {time.time() - start:.2f} seconds")


if __name__ == "__main__":
    main()
