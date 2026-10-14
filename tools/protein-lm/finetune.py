"""LoRA/AdaLoRA fine-tune a protein language model encoder for classification or regression.

Generalizes the pipeline's three LoRA scripts (ESM full-layer, ESM last-layer,
ProtT5) into one family-agnostic tool: model family is auto-detected (see
model_utils.detect_family) and determines the target_modules and pooling
function; which transformer layers the adapter is applied to (--lora_layer_scope)
is an orthogonal choice that applies to any family. Evaluation always uses a
single FASTA + labels table, split internally by --val_fraction (same pattern
as predict.py's train/test split).

Exit codes: see errors.py.
"""

import argparse
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from embed import filter_by_length, tokenize_batch
from errors import fail, run_main
from fasta_utils import check_id_overlap, ensure_file_exists, load_labels_or_fail, parse_fasta_or_fail, warn_nonstandard_residues
from model_utils import detect_family, hidden_size, load_model, resolve_device, resolve_max_length, resolve_model_source
from peft import AdaLoraConfig, get_peft_model, LoraConfig, TaskType
from predict import performance_metrics, prediction_plot
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader, Dataset

SEED = 42

_TARGET_MODULES = {"esm": ["query", "value"], "protbert": ["query", "value"], "prott5": ["q", "v"]}
# peft's layers_to_transform needs an explicit layers_pattern for T5: its modules are
# named "encoder.block.<N>.layer.<0|1>...", and peft's default pattern list (which
# matches a literal "layer" segment) picks the inner 0/1 sub-block index instead of
# the real transformer-block index N. ESM/BERT use "encoder.layer.<N>...", which
# already matches peft's default, so they need no override.
_LAYERS_PATTERN = {"esm": None, "protbert": None, "prott5": "block"}
_LEADING_SPECIAL_TOKENS = {"esm": 1, "protbert": 1, "prott5": 0}


def num_transformer_layers(encoder) -> int:
    return getattr(encoder.config, "num_hidden_layers", None) or getattr(encoder.config, "num_layers")


def resolve_layers_to_transform(scope: str, last_n: int, custom_layers: str, num_layers: int):
    """Turn --lora_layer_scope (+ its companion param) into peft's layers_to_transform.

    Returns None for "all" (every layer); a list of 0-indexed layer numbers otherwise.
    """
    if scope == "all":
        return None
    if scope == "last_n":
        if last_n < 1 or last_n > num_layers:
            fail(
                f"--lora_last_n must be between 1 and this model's layer count ({num_layers}); got {last_n}.",
                exit_code=7,
            )
        return list(range(num_layers - last_n, num_layers))
    try:
        layers = sorted({int(x.strip()) for x in custom_layers.split(",") if x.strip()})
    except ValueError:
        fail(f"--lora_custom_layers must be a comma-separated list of integers; got '{custom_layers}'.", exit_code=7)
    if not layers:
        fail("--lora_custom_layers must list at least one layer index.", exit_code=7)
    out_of_range = [i for i in layers if i < 0 or i >= num_layers]
    if out_of_range:
        fail(
            f"--lora_custom_layers contains out-of-range index(es) {out_of_range}; "
            f"this model has layers 0..{num_layers - 1}.",
            exit_code=7,
        )
    return layers


class PooledHead(nn.Module):
    """Encoder (LoRA-wrapped) -> mean-pool over real residues -> linear head."""

    def __init__(self, encoder, family: str, hidden_dim: int, output_dim: int):
        super().__init__()
        self.encoder = encoder
        self.family = family
        self.head = nn.Linear(hidden_dim, output_dim)

    def forward(self, input_ids, attention_mask):
        out = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        leading = _LEADING_SPECIAL_TOKENS[self.family]
        mask = attention_mask.clone().float()
        mask[:, :leading] = 0.0
        lengths = attention_mask.sum(dim=1)
        for i in range(mask.size(0)):
            eos_pos = int(lengths[i].item()) - 1
            if eos_pos >= leading:
                mask[i, eos_pos] = 0.0
        mask = mask.unsqueeze(-1)
        pooled = (out.last_hidden_state * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1e-9)
        return self.head(pooled)


class SeqDataset(Dataset):
    def __init__(self, seqs, targets, tokenizer, family, max_length):
        self.seqs = seqs
        self.targets = targets
        self.tokenizer = tokenizer
        self.family = family
        self.max_length = max_length

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        encoded = tokenize_batch(self.tokenizer, self.family, [self.seqs[idx]], self.max_length)
        return {
            "input_ids": encoded["input_ids"].squeeze(0),
            "attention_mask": encoded["attention_mask"].squeeze(0),
            "target": self.targets[idx],
        }


def make_collate_fn(pad_token_id: int):
    def _collate(batch):
        input_ids = nn.utils.rnn.pad_sequence(
            [b["input_ids"] for b in batch], batch_first=True, padding_value=pad_token_id
        )
        attention_mask = nn.utils.rnn.pad_sequence(
            [b["attention_mask"] for b in batch], batch_first=True, padding_value=0
        )
        targets = torch.stack([torch.as_tensor(b["target"]) for b in batch])
        return {"input_ids": input_ids, "attention_mask": attention_mask, "targets": targets}

    return _collate


def load_fasta_labels(fasta_path, labels_path, label_col, max_length, on_long_sequence):
    """Parse+validate the FULL fasta/labels pair, then length-filter the survivors.

    Mirrors predict.py's build_xy: the ID-overlap check runs on the full FASTA
    before filtering, so a sequence --on_long_sequence skip legitimately drops
    is never mistaken for a genuine FASTA/labels mismatch.
    """
    ensure_file_exists(fasta_path, "FASTA file")
    ensure_file_exists(labels_path, "Labels table")
    all_ids, all_seqs = parse_fasta_or_fail(fasta_path)
    labels = load_labels_or_fail(labels_path, id_col="id", label_cols=[label_col])
    fasta_only, labels_only = check_id_overlap(all_ids, labels.keys())
    if fasta_only or labels_only:
        fail(
            "Sequence IDs in the FASTA file and labels table do not match. "
            f"In FASTA but not labels: {fasta_only[:10]}. In labels but not FASTA: {labels_only[:10]}.",
            exit_code=6,
        )
    for warning in warn_nonstandard_residues(all_ids, all_seqs):
        print(f"WARNING: {warning}", file=sys.stderr)
    ids, seqs, length_warnings = filter_by_length(all_ids, all_seqs, max_length, on_long_sequence)
    for warning in length_warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if not seqs:
        fail(f"No sequences remain in {fasta_path} after length filtering.", exit_code=5)
    raw_values = [labels[seq_id][0] for seq_id in ids]
    return ids, seqs, raw_values


def train_epoch(model, loader, optimizer, loss_fn, device, adapter_method, global_step):
    """Returns (avg_loss, new_global_step). adapter_method/global_step drive AdaLoRA's
    rank reallocation, which peft requires to be called every step: backward -> step ->
    update_and_allocate -> zero_grad (see AdaLoraModel.update_and_allocate's docstring).
    """
    model.train()
    total_loss = 0.0
    for batch in loader:
        logits = model(batch["input_ids"].to(device), batch["attention_mask"].to(device))
        loss = loss_fn(logits, batch["targets"].to(device))
        loss.backward()
        optimizer.step()
        if adapter_method == "adalora":
            model.encoder.base_model.update_and_allocate(global_step)
        optimizer.zero_grad()
        total_loss += loss.item() * batch["input_ids"].size(0)
        global_step += 1
    return total_loss / len(loader.dataset), global_step


@torch.no_grad()
def predict_epoch(model, loader, device, task):
    """Return (predictions, scores). scores is the softmax probability matrix
    for classification (needed by prediction_plot's PR curve), None for regression.
    """
    model.eval()
    all_preds, all_scores = [], []
    for batch in loader:
        logits = model(batch["input_ids"].to(device), batch["attention_mask"].to(device))
        if task == "classification":
            probs = torch.softmax(logits, dim=1)
            all_preds.append(torch.argmax(probs, dim=1).cpu().numpy())
            all_scores.append(probs.cpu().numpy())
        else:
            all_preds.append(logits.squeeze(-1).cpu().numpy())
    predictions = np.concatenate(all_preds)
    scores = np.concatenate(all_scores) if task == "classification" else None
    return predictions, scores


def main() -> None:
    parser = argparse.ArgumentParser(description="LoRA/AdaLoRA fine-tune a protein LM encoder")
    parser.add_argument("--model_path", default=None, help="Preinstalled model folder (from the huggingface data table)")
    parser.add_argument("--model_archive", default=None, help="Uploaded model archive (tar/tgz)")
    parser.add_argument("--task", choices=["classification", "regression"], required=True)
    parser.add_argument("--data_fasta", required=True)
    parser.add_argument("--data_labels", required=True)
    parser.add_argument("--label_col", required=True)
    parser.add_argument("--val_fraction", type=float, default=0.2)
    parser.add_argument("--random_seed", type=int, default=SEED)
    parser.add_argument("--adapter_method", choices=["lora", "adalora"], default="lora")
    parser.add_argument("--lora_layer_scope", choices=["all", "last_n", "custom"], default="all")
    parser.add_argument("--lora_last_n", type=int, default=1)
    parser.add_argument("--lora_custom_layers", default="")
    parser.add_argument("--lora_r", type=int, default=8)
    parser.add_argument("--lora_alpha", type=int, default=16)
    parser.add_argument("--lora_dropout", type=float, default=0.05)
    parser.add_argument("--adalora_init_r", type=int, default=12)
    parser.add_argument("--adalora_target_r", type=int, default=8)
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--weight_decay", type=float, default=0.01)
    parser.add_argument("--max_length", type=int, default=1024)
    parser.add_argument("--on_long_sequence", choices=["skip", "truncate"], default="truncate")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = parser.parse_args()
    args.model_dir = resolve_model_source(args.model_path, args.model_archive)

    device = resolve_device(args.device)
    family = detect_family(args.model_dir)
    _, tokenizer, encoder = load_model(args.model_dir)
    args.max_length = resolve_max_length(encoder, args.max_length)

    all_ids, all_seqs, all_raw = load_fasta_labels(
        args.data_fasta, args.data_labels, args.label_col, args.max_length, args.on_long_sequence
    )

    label_encoder = None
    if args.task == "classification":
        label_encoder = LabelEncoder()
        all_targets = label_encoder.fit_transform(all_raw)
        output_dim = len(label_encoder.classes_)
        loss_fn = nn.CrossEntropyLoss()
        stratify = all_targets
        _, counts = np.unique(all_targets, return_counts=True)
        if counts.min() < 2:
            fail(
                "Cannot split into train/validation: every class must have at least 2 samples "
                f"(smallest class has {counts.min()}). Use more data, or a larger validation "
                "fraction so every class can appear on both sides of the split.",
                exit_code=7,
            )
    else:
        try:
            all_targets = np.array([float(v) for v in all_raw], dtype=np.float32)
        except ValueError as exc:
            fail(f"Labels table column '{args.label_col}' contains a non-numeric value: {exc}", exit_code=7)
        output_dim = 1
        loss_fn = nn.MSELoss()
        stratify = None

    indices = np.arange(len(all_ids))
    try:
        train_idx, val_idx = train_test_split(
            indices, test_size=args.val_fraction, random_state=args.random_seed, stratify=stratify
        )
    except ValueError as exc:
        fail(f"Cannot split into train/validation with val_fraction={args.val_fraction}: {exc}", exit_code=7)

    train_seqs = [all_seqs[i] for i in train_idx]
    val_seqs = [all_seqs[i] for i in val_idx]
    val_ids = [all_ids[i] for i in val_idx]
    val_raw = [all_raw[i] for i in val_idx]

    if args.task == "classification":
        train_targets_t = torch.as_tensor(all_targets[train_idx], dtype=torch.long)
    else:
        train_targets_t = torch.as_tensor(all_targets[train_idx], dtype=torch.float32).unsqueeze(-1)

    num_layers = num_transformer_layers(encoder)
    layers_to_transform = resolve_layers_to_transform(
        args.lora_layer_scope, args.lora_last_n, args.lora_custom_layers, num_layers
    )
    layers_pattern = _LAYERS_PATTERN[family] if layers_to_transform is not None else None
    common_config = dict(
        task_type=TaskType.FEATURE_EXTRACTION,
        target_modules=_TARGET_MODULES[family],
        layers_to_transform=layers_to_transform,
        layers_pattern=layers_pattern,
        lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout,
        bias="none",
    )
    if args.adapter_method == "adalora":
        steps_per_epoch = max(1, -(-len(train_seqs) // args.batch_size))  # ceil division
        total_step = steps_per_epoch * args.epochs
        adapter_config = AdaLoraConfig(
            init_r=args.adalora_init_r, target_r=args.adalora_target_r, total_step=total_step, **common_config
        )
    else:
        adapter_config = LoraConfig(r=args.lora_r, **common_config)
    encoder = get_peft_model(encoder, adapter_config)
    encoder.print_trainable_parameters()

    model = PooledHead(encoder, family, hidden_size(encoder.get_base_model()), output_dim).to(device)

    pad_token_id = tokenizer.pad_token_id or 0
    collate = make_collate_fn(pad_token_id)
    train_loader = DataLoader(
        SeqDataset(train_seqs, train_targets_t, tokenizer, family, args.max_length),
        batch_size=args.batch_size, shuffle=True, collate_fn=collate,
    )

    optimizer = torch.optim.AdamW(
        [p for p in model.parameters() if p.requires_grad], lr=args.lr, weight_decay=args.weight_decay
    )

    global_step = 0
    for epoch in range(1, args.epochs + 1):
        loss, global_step = train_epoch(model, train_loader, optimizer, loss_fn, device, args.adapter_method, global_step)
        print(f"epoch {epoch}/{args.epochs}  loss={loss:.4f}", file=sys.stderr)

    val_loader = DataLoader(
        SeqDataset(val_seqs, [0] * len(val_seqs), tokenizer, family, args.max_length),
        batch_size=args.batch_size, shuffle=False, collate_fn=collate,
    )
    predicted, scores = predict_epoch(model, val_loader, device, args.task)

    if args.task == "classification":
        y_true = label_encoder.transform(val_raw)
        pred_labels = label_encoder.inverse_transform(predicted)
    else:
        y_true = all_targets[val_idx]
        pred_labels = predicted

    pd.DataFrame([performance_metrics(y_true, predicted, args.task)]).to_csv(
        "val_metrics.tsv", sep="\t", index=False
    )
    prediction_plot(y_true, predicted, args.task, scores)
    pd.DataFrame({"id": val_ids, "true_label": val_raw, "predicted_label": pred_labels}).to_csv(
        "output_predicted_data", sep="\t", index=False
    )

    model.encoder.save_pretrained("output_adapter", safe_serialization=True)
    print("Saved LoRA adapter to output_adapter/", file=sys.stderr)


if __name__ == "__main__":
    run_main(main)
