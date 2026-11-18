"""Download the models of the Galaxy WhisperX tool into one directory.

Run it inside the tool's container, so it fetches exactly the models that
WhisperX version loads (its alignment model per language and its diarization
model), from the same loaders the tool uses:

    docker run --rm --user "$(id -u):$(id -g)" \\
        -v /data/db/whisperx_models:/models \\
        -v "$PWD/download_models.py:/download_models.py:ro" \\
        -e HF_TOKEN=<token> \\
        quay.io/galaxy/whisperx:3.8.6 python3 /download_models.py --model-dir /models

With rootless Docker, leave out --user. The models need about 60 GB, mostly
the 35 Hugging Face alignment models. The token needs the terms of the
pyannote diarization model accepted. The directory is the `path` of the
whisperx rows in the huggingface data table.
"""

import argparse
import gc
import os
import sys
from pathlib import Path

# The value column of the rows in tool-data/huggingface.loc.sample
WHISPER_MODELS = [
    "mobiuslabsgmbh/faster-whisper-large-v3-turbo",
    "Systran/faster-whisper-tiny",
    "Systran/faster-whisper-base",
    "Systran/faster-whisper-small",
    "Systran/faster-whisper-medium",
    "Systran/faster-whisper-large-v3",
]


def main():
    parser = argparse.ArgumentParser(
        description="Download the models of the Galaxy WhisperX tool"
    )
    parser.add_argument(
        "--model-dir", required=True, help="Directory to store the models"
    )
    parser.add_argument(
        "--hf-token",
        default=os.getenv("HF_TOKEN") or os.getenv("HF_AUTH_TOKEN"),
        help="Hugging Face token for the gated pyannote model (default: $HF_TOKEN)",
    )
    args = parser.parse_args()
    if not args.hf_token:
        # checked first, the gated diarization model is downloaded last
        parser.error("set HF_TOKEN or --hf-token (needed for the pyannote model)")

    root = Path(args.model_dir).resolve()
    hub = root / "hub"  # the tool passes this as --model_dir
    hub.mkdir(parents=True, exist_ok=True)
    # anything a loader caches without an explicit path also lands in root
    os.environ["HF_HOME"] = str(root)
    os.environ["TORCH_HOME"] = str(root)
    # NLTK only writes below NLTK_DATA (the tool sets the same path)
    os.environ["NLTK_DATA"] = str(root / "nltk_data")
    # otherwise transformers also fetches a safetensors copy of .bin-only
    # models from the Hub's conversion PR, which offline loading never uses
    os.environ["DISABLE_SAFETENSORS_CONVERSION"] = "true"
    # pyannote 4 sends usage metrics to otel.pyannote.ai unless disabled
    os.environ["PYANNOTE_METRICS_ENABLED"] = "false"

    import nltk
    from whisperx.alignment import (
        DEFAULT_ALIGN_MODELS_HF,
        DEFAULT_ALIGN_MODELS_TORCH,
        load_align_model,
    )
    from whisperx.asr import load_model
    from whisperx.diarize import DiarizationPipeline

    failed = []

    def fetch(name, load):
        # loading the model once proves the download is complete
        print(f"Downloading {name} ...", flush=True)
        try:
            load()
            print(f"OK: {name}", flush=True)
        except Exception as e:
            print(f"Failed: {name}: {e}", flush=True)
            failed.append(name)
        gc.collect()

    for model in WHISPER_MODELS:
        fetch(
            model,
            lambda m=model: load_model(
                m, "cpu", compute_type="int8", download_root=str(hub)
            ),
        )
    for lang in sorted({**DEFAULT_ALIGN_MODELS_TORCH, **DEFAULT_ALIGN_MODELS_HF}):
        fetch(
            f"alignment model ({lang})",
            lambda c=lang: load_align_model(c, "cpu", model_dir=str(hub)),
        )
    fetch(
        "diarization model",
        lambda: DiarizationPipeline(
            token=args.hf_token, device="cpu", cache_dir=str(hub)
        ),
    )
    fetch(
        "NLTK punkt_tab",
        lambda: nltk.download(
            "punkt_tab",
            download_dir=str(root / "nltk_data"),
            quiet=True,
            raise_on_error=True,
        ),
    )

    if failed:
        print(f"Failed to download {len(failed)} model(s): {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
