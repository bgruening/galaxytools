import argparse
import os
import sys
from pathlib import Path
from sentence_transformers import SentenceTransformer


def main():
    parser = argparse.ArgumentParser(
        description="Download multiple SentenceTransformer models."
    )
    parser.add_argument(
        "--models",
        nargs="+",
        required=True,
        help="Hugging Face model IDs, separated by spaces",
    )

    args = parser.parse_args()

    base_path = os.environ.get("HF_MODEL_DIR")
    if not base_path:
        print("Error: HF_MODEL_DIR environment variable is not set.")
        sys.exit(1)

    base_path = Path(base_path)

    for model_id in args.models:
        out_dir = base_path / model_id
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n[START] Downloading: {model_id}")
        print(f"[PATH] Destination: {out_dir}")

        try:
            model = SentenceTransformer(model_id)
            model.save(str(out_dir))
            print(f"[SUCCESS] Saved {model_id}")
        except Exception as e:
            print(f"[ERROR] Could not download {model_id}: {e}")

    print("\nAll tasks completed.")


if __name__ == "__main__":
    main()
