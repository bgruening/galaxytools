"""Regenerate the tiny test-data checkpoints (deterministic, real architectures).

Each checkpoint pairs the REAL tokenizer files for its family (copied from a
local Hugging Face model folder -- already tiny, a few KB to ~230KB for
ProtT5's sentencepiece model) with a freshly random-initialized model at
minimal dimensions (real transformers config class, not a mock). This keeps
every committed file under the repo CI's 1MB file_sizes limit while still
exercising the real tokenizer vocabulary and the real model_type detection
path (config.json's model_type is read by model_utils.detect_family).

Run from the repo root:
    conda run -n protein_lm_galaxy python tools/protein_lm/test-data/make_tiny_models.py

Source tokenizer folders (already downloaded locally by this project):
    models/esm/esm2_t33_650M        (vocab_size=33)
    models/protbert/prot_bert_bfd   (vocab_size=30)
    models/prott5/prot_t5_xl_uniref50 (vocab_size=128)
"""

import shutil
from pathlib import Path

import torch
from transformers import (
    AutoTokenizer,
    BertConfig,
    BertModel,
    BertTokenizer,
    EsmConfig,
    EsmModel,
    T5Config,
    T5EncoderModel,
    T5Tokenizer,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent

SEED = 42

TOKENIZER_FILES = {
    "esm": ["vocab.txt", "tokenizer_config.json", "special_tokens_map.json"],
    "protbert": ["vocab.txt", "tokenizer_config.json", "special_tokens_map.json"],
    "prott5": ["spiece.model", "tokenizer_config.json", "special_tokens_map.json"],
}
SOURCE_DIRS = {
    "esm": REPO_ROOT / "models/esm/esm2_t33_650M",
    "protbert": REPO_ROOT / "models/protbert/prot_bert_bfd",
    "prott5": REPO_ROOT / "models/prott5/prot_t5_xl_uniref50",
}


def copy_tokenizer_files(family: str, dest: Path) -> None:
    src = SOURCE_DIRS[family]
    for filename in TOKENIZER_FILES[family]:
        source_path = src / filename
        if source_path.exists():
            shutil.copy(source_path, dest / filename)


def make_tiny_esm(dest: Path) -> None:
    torch.manual_seed(SEED)
    config = EsmConfig(
        vocab_size=33,
        hidden_size=16,
        num_hidden_layers=2,
        num_attention_heads=2,
        intermediate_size=32,
        max_position_embeddings=1026,
        pad_token_id=1,
    )
    model = EsmModel(config)
    model.save_pretrained(dest, safe_serialization=True)
    copy_tokenizer_files("esm", dest)
    AutoTokenizer.from_pretrained(dest)  # sanity check: loads cleanly


def make_tiny_protbert(dest: Path) -> None:
    torch.manual_seed(SEED)
    config = BertConfig(
        vocab_size=30,
        hidden_size=16,
        num_hidden_layers=2,
        num_attention_heads=2,
        intermediate_size=32,
        max_position_embeddings=1026,
        pad_token_id=0,
    )
    model = BertModel(config, add_pooling_layer=False)
    model.save_pretrained(dest, safe_serialization=True)
    copy_tokenizer_files("protbert", dest)
    BertTokenizer.from_pretrained(dest, do_lower_case=False)


def make_tiny_prott5(dest: Path) -> None:
    torch.manual_seed(SEED)
    config = T5Config(
        vocab_size=128,
        d_model=16,
        d_ff=32,
        num_layers=2,
        num_heads=2,
        d_kv=8,
    )
    model = T5EncoderModel(config)
    model.save_pretrained(dest, safe_serialization=True)
    copy_tokenizer_files("prott5", dest)
    T5Tokenizer.from_pretrained(dest, do_lower_case=False)


def report_sizes(dest: Path) -> None:
    for path in sorted(dest.rglob("*")):
        if path.is_file():
            print(f"  {path.relative_to(OUT_DIR)}: {path.stat().st_size:,} bytes")


def main() -> None:
    builders = {"tiny_esm": make_tiny_esm, "tiny_protbert": make_tiny_protbert, "tiny_prott5": make_tiny_prott5}
    for name, builder in builders.items():
        dest = OUT_DIR / name
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        builder(dest)
        print(f"\n{name}/")
        report_sizes(dest)


if __name__ == "__main__":
    main()
