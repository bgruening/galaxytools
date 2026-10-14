"""Model-family detection, offline loading, and archive extraction.

All three model families are loaded through Hugging Face ``transformers``
from a local folder only -- no network access. ``HF_HUB_OFFLINE`` /
``TRANSFORMERS_OFFLINE`` are set at import time, and every ``from_pretrained``
call also passes ``local_files_only=True``. Every call also passes
``use_safetensors=True``: the tool must never deserialise a ``pytorch_model.bin``
(a pickle) from an uploaded archive, since that is an arbitrary-code-execution
risk for inputs we do not trust (see PR description for the CVE this guards
against on older torch/transformers combinations).
"""

import json
import os
import sys
import tarfile
from pathlib import Path

os.environ["HF_HUB_OFFLINE"] = "1"
os.environ["TRANSFORMERS_OFFLINE"] = "1"

import torch  # noqa: E402  (offline env vars must be set before importing transformers)
from errors import fail  # noqa: E402
from transformers import (  # noqa: E402
    AutoTokenizer,
    BertModel,
    BertTokenizer,
    EsmModel,
    PreTrainedModel,
    PreTrainedTokenizerBase,
    T5EncoderModel,
    T5Tokenizer,
)

_FAMILY_BY_MODEL_TYPE = {"esm": "esm", "bert": "protbert", "t5": "prott5"}


def resolve_device(requested: str = "auto") -> torch.device:
    if requested == "cpu":
        return torch.device("cpu")
    if requested == "cuda" and not torch.cuda.is_available():
        fail("GPU requested but no CUDA device is available.", exit_code=2)
    if requested == "cuda":
        return torch.device("cuda")
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def extract_archive(archive_path: str, dest_dir: str) -> str:
    """Extract a tar/tgz model archive to dest_dir, handling an optional
    single top-level wrapper directory. Returns the folder that directly
    contains config.json.
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive_path) as tar:
            members = tar.getmembers()
            if not members:
                fail(f"Uploaded model archive is empty: {archive_path}", exit_code=3)
            tar.extractall(dest, filter="data")
    except tarfile.TarError as exc:
        fail(f"Uploaded file is not a valid tar/tgz archive: {exc}", exit_code=3)

    if (dest / "config.json").exists():
        return str(dest)

    # No config.json at the top level -- look for exactly one wrapper directory.
    entries = [p for p in dest.iterdir() if p.name != "__MACOSX"]
    subdirs = [p for p in entries if p.is_dir()]
    if len(subdirs) == 1 and (subdirs[0] / "config.json").exists():
        return str(subdirs[0])

    fail(
        f"Extracted archive does not contain a config.json (checked {dest} "
        "and its single top-level subdirectory, if any).",
        exit_code=2,
    )


def resolve_model_source(model_path: str | None, model_archive: str | None) -> str:
    """Resolve the --model_path (preinstalled) / --model_archive (upload) pair
    that every entrypoint script's argparse exposes into one local model_dir.
    """
    if model_archive:
        return extract_archive(model_archive, "./extracted_model")
    if not model_path:
        fail("No model was selected: provide either a preinstalled model path or an uploaded archive.", exit_code=2)
    if not Path(model_path).exists():
        fail(f"Model path does not exist: {model_path}", exit_code=2)
    return model_path


def detect_family(model_dir: str) -> str:
    """Read config.json's model_type and map it to esm / protbert / prott5."""
    config_path = Path(model_dir) / "config.json"
    if not config_path.exists():
        fail(f"Model path does not exist or has no config.json: {model_dir}", exit_code=2)
    try:
        config = json.loads(config_path.read_text())
    except json.JSONDecodeError as exc:
        fail(f"config.json at {model_dir} is not valid JSON: {exc}", exit_code=2)

    model_type = config.get("model_type")
    family = _FAMILY_BY_MODEL_TYPE.get(model_type)
    if family is None:
        fail(
            f"Unsupported model_type '{model_type}' in {config_path} "
            f"(supported: {sorted(_FAMILY_BY_MODEL_TYPE)}).",
            exit_code=2,
        )
    return family


def load_model(model_dir: str) -> tuple[str, PreTrainedTokenizerBase, PreTrainedModel]:
    """Detect the model family and load (family, tokenizer, model) from model_dir.

    Every from_pretrained call is offline-only and safetensors-only.
    """
    family = detect_family(model_dir)
    common = {"local_files_only": True, "use_safetensors": True}

    if family == "esm":
        tokenizer = AutoTokenizer.from_pretrained(model_dir, **common)
        model = EsmModel.from_pretrained(model_dir, **common)
    elif family == "protbert":
        tokenizer = BertTokenizer.from_pretrained(model_dir, do_lower_case=False, **common)
        model = BertModel.from_pretrained(model_dir, **common)
    else:
        tokenizer = T5Tokenizer.from_pretrained(model_dir, do_lower_case=False, **common)
        model = T5EncoderModel.from_pretrained(model_dir, **common)

    return family, tokenizer, model


def hidden_size(model: PreTrainedModel) -> int:
    config = model.config
    return getattr(config, "hidden_size", None) or getattr(config, "d_model")


def resolve_max_length(model: PreTrainedModel, requested_max_length: int) -> int:
    """Clamp requested_max_length to the model's absolute position-embedding
    capacity, if it has one (BERT/ESM; T5 uses relative position biases and
    has no such limit). Without this, a sequence that survives --on_long_sequence
    truncation but still exceeds the model's actual capacity would crash deep in
    the forward pass with a raw PyTorch shape-mismatch error instead of a clear
    tool message.
    """
    model_max = getattr(model.config, "max_position_embeddings", None)
    if model_max is None or requested_max_length <= model_max:
        return requested_max_length
    print(
        f"WARNING: requested max_length={requested_max_length} exceeds this model's "
        f"max_position_embeddings={model_max}; clamping to {model_max}.",
        file=sys.stderr,
    )
    return model_max
