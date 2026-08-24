import os
import csv
import json
import logging
from dataclasses import dataclass, field
from typing import Any, Optional, Dict, Sequence, Tuple, List, Union

import torch
import transformers
import sklearn
import numpy as np
from torch.utils.data import Dataset
from scipy.stats import pearsonr, spearmanr

from peft import (
    LoraConfig,
    get_peft_model,
)

# -----------------------------
# Arguments
# -----------------------------
@dataclass
class ModelArguments:
    model_name_or_path: Optional[str] = field(default="zhihan1996/DNABERT-2-117M")

    problem_type: str = field(
        default="classification",
        metadata={"help": 'Task type: one of {"classification", "regression"}'}
    )

    use_lora: bool = field(
        default=False,
        metadata={"help": "Whether to enable LoRA fine-tuning"}
    )

    # LoRA args
    lora_r: int = field(default=8, metadata={"help": "hidden dimension for LoRA"})
    lora_alpha: int = field(default=32, metadata={"help": "alpha for LoRA"})
    lora_dropout: float = field(default=0.05, metadata={"help": "dropout rate for LoRA"})
    lora_target_modules: str = field(default="Wqkv,wo", metadata={"help": "where to perform LoRA"})


@dataclass
class DataArguments:
    train_file: str = field(default=None, metadata={"help": "Path to train data"})
    eval_file: str = field(default=None, metadata={"help": "Path to eval/dev data"})
    test_file: str = field(default=None, metadata={"help": "Path to test data"})


@dataclass
class TrainingArguments(transformers.TrainingArguments):
    cache_dir: Optional[str] = field(default=None)
    run_name: str = field(default="run")
    optim: str = field(default="adamw_torch")
    model_max_length: int = field(default=512, metadata={"help": "Maximum sequence length."})
    gradient_accumulation_steps: int = field(default=1)
    per_device_train_batch_size: int = field(default=1)
    per_device_eval_batch_size: int = field(default=1)
    num_train_epochs: int = field(default=1)
    fp16: bool = field(default=False)
    logging_steps: int = field(default=100)
    save_steps: int = field(default=100)
    eval_steps: int = field(default=100)
    eval_strategy: str = field(default="steps")
    save_strategy: str = field(default="steps")
    warmup_steps: int = field(default=50)
    weight_decay: float = field(default=0.01)
    learning_rate: float = field(default=1e-4)
    save_total_limit: int = field(default=3)
    load_best_model_at_end: bool = field(default=True)
    output_dir: str = field(default="output")
    find_unused_parameters: bool = field(default=False)
    checkpointing: bool = field(default=False)
    dataloader_pin_memory: bool = field(default=False)
    eval_and_save_results: bool = field(default=True)
    save_model: bool = field(default=False)
    seed: int = field(default=42)

# -----------------------------
# Dataset
# -----------------------------
class SupervisedDataset(Dataset):
    """
    Dataset for supervised fine-tuning.
    """
    def __init__(
        self,
        data_path: str,
        tokenizer: transformers.PreTrainedTokenizer,
        problem_type: str = "classification"
    ):
        super(SupervisedDataset, self).__init__()

        with open(data_path, "r") as f:
            data = list(csv.reader(f, delimiter="\t"))[1:]  # skip header

        if len(data[0]) == 2:
            # data is in the format of [text, label]
            logging.warning("Perform single sequence task...")
            texts = [d[0] for d in data]
            labels = [int(d[1]) for d in data] if problem_type == "classification" else [float(d[1]) for d in data]
        else:
            raise ValueError("Data format not supported.")

        output = tokenizer(
            texts,
            return_tensors="pt",
            padding="longest",
            max_length=tokenizer.model_max_length,
            truncation=True,
        )

        self.input_ids = output["input_ids"]
        self.attention_mask = output["attention_mask"]
        self.labels = labels
        self.problem_type = problem_type
        self.num_labels = len(set(labels)) if problem_type == "classification" else 1

    def __len__(self):
        return len(self.input_ids)

    def __getitem__(self, i) -> Dict[str, torch.Tensor]:
        return dict(input_ids=self.input_ids[i], labels=self.labels[i])


@dataclass
class DataCollatorForSupervisedDataset:
    """Collate examples for supervised fine-tuning."""

    tokenizer: transformers.PreTrainedTokenizer
    problem_type: str = "classification"

    def __call__(self, instances: Sequence[Dict]) -> Dict[str, torch.Tensor]:
        input_ids, labels = tuple([instance[key] for instance in instances] for key in ("input_ids", "labels"))
        input_ids = torch.nn.utils.rnn.pad_sequence(
            input_ids, batch_first=True, padding_value=self.tokenizer.pad_token_id
        )
        if self.problem_type == "classification":
            labels = torch.tensor(labels).long()
        else:
            labels = torch.tensor(labels).float()
        return dict(
            input_ids=input_ids,
            labels=labels,
            attention_mask=input_ids.ne(self.tokenizer.pad_token_id),
        )


# -----------------------------
# Metrics
# -----------------------------
def calculate_classification_metrics(predictions: np.ndarray, labels: np.ndarray):
    valid_mask = labels != -100
    valid_predictions = predictions[valid_mask]
    valid_labels = labels[valid_mask]

    results = {
        "accuracy": sklearn.metrics.accuracy_score(valid_labels, valid_predictions),
        "f1": sklearn.metrics.f1_score(valid_labels, valid_predictions, average="macro", zero_division=0),
        "matthews_correlation": sklearn.metrics.matthews_corrcoef(valid_labels, valid_predictions),
        "precision": sklearn.metrics.precision_score(valid_labels, valid_predictions, average="macro", zero_division=0),
        "recall": sklearn.metrics.recall_score(valid_labels, valid_predictions, average="macro", zero_division=0),
    }
    return results


def calculate_regression_metrics(predictions: np.ndarray, labels: np.ndarray):
    valid_mask = labels != -100
    preds = predictions[valid_mask]
    targets = labels[valid_mask]

    mse = sklearn.metrics.mean_squared_error(targets, preds)
    rmse = np.sqrt(mse)
    mae = sklearn.metrics.mean_absolute_error(targets, preds)
    mape = sklearn.metrics.mean_absolute_percentage_error(targets, preds)
    r2 = sklearn.metrics.r2_score(targets, preds)

    try:
        pearson, pearson_p = pearsonr(targets, preds)
    except Exception:
        pearson = 0.0
        pearson_p = 1.0

    try:
        spearman, spearman_p = spearmanr(targets, preds)
    except Exception:
        spearman = 0.0
        spearman_p = 1.0

    results =  {
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "mape": mape,
        "r2": r2,
        "pearson": pearson,
        "pearson_p": pearson_p,
        "spearman": spearman,
        "spearman_p": spearman_p,
    }
    return results


def make_preprocess_logits_for_metrics(problem_type: str):
    def _preprocess_logits_for_metrics(logits: Union[torch.Tensor, Tuple[torch.Tensor, Any]], _):
        if isinstance(logits, tuple):
            logits = logits[0]
        if problem_type == "classification":
            if logits.ndim == 3:
                logits = logits.reshape(-1, logits.shape[-1])
            return logits
        else:
            if logits.ndim > 1:
                logits = logits.squeeze(-1)
            return logits
    return _preprocess_logits_for_metrics


def make_compute_metrics(problem_type: str):
    def _compute_metrics(eval_pred):
        logits, labels = eval_pred
        labels = np.squeeze(labels)

        if problem_type == "classification":
            probs = torch.softmax(torch.tensor(logits), dim=-1).numpy()
            predictions = np.argmax(probs, axis=-1)
            results = calculate_classification_metrics(predictions, labels)

            valid_mask = labels != -100
            valid_labels = labels[valid_mask]

            if probs.shape[-1] == 2 and len(np.unique(valid_labels)) > 1: # for binary classification
                results["roc_auc"] = sklearn.metrics.roc_auc_score(valid_labels, probs[valid_mask, 1])
                results["pr_auc"] = sklearn.metrics.average_precision_score(valid_labels, probs[valid_mask, 1])
            return results
        else:
            predictions = np.squeeze(logits)
            return calculate_regression_metrics(predictions, labels)
    return _compute_metrics


# -----------------------------
# Prediction dump
# -----------------------------
def get_texts(data_path: str):
    texts = []
    with open(data_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)
        for row in reader:
            texts.append(row[0])
    return texts

def dump_test_predictions(trainer: transformers.Trainer, test_dataset: Dataset, sequences: list[str],
                          output_dir: str, problem_type: str):
    pred_output = trainer.predict(test_dataset=test_dataset)
    labels = pred_output.label_ids
    preds_raw = pred_output.predictions

    if problem_type == "classification":
        probs = torch.softmax(torch.tensor(preds_raw), dim=-1).numpy()
        preds = np.argmax(probs, axis=-1)
    else:
        preds = np.squeeze(preds_raw)

    labels = np.squeeze(labels)
    preds = np.squeeze(preds)

    pred_path = os.path.join(output_dir, "test_predictions.csv")
    with open(pred_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence", "label", "prediction"])
        for seq, y_true, y_pred in zip(sequences, labels.tolist(), preds.tolist()):
            writer.writerow([seq, y_true, y_pred])

# -----------------------------
# Train
# -----------------------------
def train():
    parser = transformers.HfArgumentParser((ModelArguments, DataArguments, TrainingArguments))
    model_args, data_args, training_args = parser.parse_args_into_dataclasses()

    problem_type = model_args.problem_type.lower().strip()
    if problem_type not in {"classification", "regression"}:
        raise ValueError(f'Invalid problem_type={model_args.problem_type}. Use one of: "classification", "regression"')

    if not data_args.train_file or not data_args.eval_file or not data_args.test_file:
        raise ValueError("You must provide --train_file, --eval_file, and --test_file.")

    os.makedirs(training_args.output_dir, exist_ok=True)

    tokenizer = transformers.AutoTokenizer.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=training_args.cache_dir,
        model_max_length=training_args.model_max_length,
        padding_side="right",
        use_fast=True,
        trust_remote_code=True,
    )

    train_dataset = SupervisedDataset(
        tokenizer=tokenizer,
        data_path=data_args.train_file,
        problem_type=problem_type,
    )
    val_dataset = SupervisedDataset(
        tokenizer=tokenizer,
        data_path=data_args.eval_file,
        problem_type=problem_type,
    )
    test_dataset = SupervisedDataset(
        tokenizer=tokenizer,
        data_path=data_args.test_file,
        problem_type=problem_type,
    )
    data_collator = DataCollatorForSupervisedDataset(tokenizer=tokenizer, problem_type=problem_type)

    config = transformers.AutoConfig.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=training_args.cache_dir,
        trust_remote_code=True,
    )
    config.pad_token_id = tokenizer.pad_token_id
    if problem_type == "classification":
        config.num_labels = train_dataset.num_labels
        config.problem_type = "single_label_classification"
    else:
        config.num_labels = 1
        config.problem_type = "regression"

    model = transformers.AutoModelForSequenceClassification.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=training_args.cache_dir,
        config=config,
        trust_remote_code=True,
    )
    model.config.pad_token_id = tokenizer.pad_token_id

    # Apply LoRA
    if model_args.use_lora:
        lora_config = LoraConfig(
            r=model_args.lora_r,
            lora_alpha=model_args.lora_alpha,
            target_modules=list(model_args.lora_target_modules.split(",")),
            lora_dropout=model_args.lora_dropout,
            bias="none",
            task_type="SEQ_CLS",
            inference_mode=False,
        )
        model = get_peft_model(model, lora_config)
        model.print_trainable_parameters()

    trainer = transformers.Trainer(
        model=model,
        args=training_args,
        preprocess_logits_for_metrics=make_preprocess_logits_for_metrics(problem_type),
        compute_metrics=make_compute_metrics(problem_type),
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        data_collator=data_collator,
    )

    trainer.train()

    # i) save full fine-tuned model if requested
    if training_args.save_model:
        model_dir = os.path.join(training_args.output_dir, "model")
        os.makedirs(model_dir, exist_ok=True)

        if model_args.use_lora:
            merged_model = trainer.model.merge_and_unload()
            merged_model.save_pretrained(model_dir)
            tokenizer.save_pretrained(model_dir)
        else:
            trainer.save_model(model_dir)

    # ii) save trainer state
    trainer.save_state()

    # iii) output test evaluation metrics
    if training_args.eval_and_save_results:
        test_metrics = trainer.evaluate(eval_dataset=test_dataset)
        metrics_path = os.path.join(training_args.output_dir, "evaluation.json")
        with open(metrics_path, "w") as f:
            json.dump(test_metrics, f, indent=2)

    # iv) output test labels and predictions
    texts = get_texts(data_path=data_args.test_file)  
    dump_test_predictions(
        trainer=trainer,
        test_dataset=test_dataset,
        sequences=texts,
        output_dir=training_args.output_dir,
        problem_type=problem_type,
    )


if __name__ == "__main__":
    train()