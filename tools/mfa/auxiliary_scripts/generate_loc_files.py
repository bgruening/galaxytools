"""Generate Galaxy data tables from downloaded MFA models."""

import argparse
from pathlib import Path


MODEL_EXTENSIONS = {
    "acoustic": ".zip",
    "dictionary": ".dict",
    "g2p": ".zip",
    "language_model": ".zip",
    "tokenizer": ".zip",
}


def generate_loc_files(models_root, output_dir):
    models_root = models_root.resolve()
    if not models_root.is_dir():
        raise ValueError(f"Model directory does not exist: {models_root}")

    tables = {}
    for model_type, extension in MODEL_EXTENSIONS.items():
        rows = ["#value\tdbkey\tname\tpath\n"]
        for model in sorted((models_root / model_type).glob(f"*{extension}")):
            if not model.is_file():
                continue
            name = model.stem.replace("_", " ").title()
            fields = (model.stem, "default", name, str(model))
            if any(any(char in field for char in "\t\r\n") for field in fields):
                raise ValueError(f"Model path contains a table delimiter: {model!s}")
            rows.append("\t".join(fields) + "\n")
        tables[model_type] = "".join(rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    for model_type, content in tables.items():
        destination = output_dir / f"mfa_{model_type}.loc"
        destination.write_text(content, encoding="utf-8")
        print(f"Created: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-root", type=Path, required=True,
                        help="Downloaded MFA pretrained_models directory")
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "tool-data",
                        help="Destination for .loc files (default: the suite's tool-data directory)")
    args = parser.parse_args()
    try:
        generate_loc_files(args.models_root, args.output_dir)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
