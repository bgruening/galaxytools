"""Download the models of the Galaxy FLUX tool into one directory.

Run it inside the tool's container, so the models are fetched by the same
diffusers version the tool uses:

    docker run --rm --user "$(id -u):$(id -g)" \\
        -v /data/db/flux_models:/models \\
        -v "$PWD/download_models.py:/download_models.py:ro" \\
        -e HF_TOKEN=<token> \\
        quay.io/galaxy/black_forest_labs_flux:2026-0 \\
        python3 /download_models.py --model-dir /models

With rootless Docker, leave out --user. Both models need about 55 GB;
--models downloads only one of them. FLUX 3 Image Klein is gated:
accept their terms on Hugging Face with the account of the token first. The directory is the `path` of the flux rows in the huggingface data
table.
"""

import argparse
import os
import sys
from pathlib import Path

BFL = "black-forest-labs/"
# The value column of the rows in tool-data/huggingface.loc.sample
MODELS = [
    BFL + "FLUX.2-klein-4B",
    BFL + "flux-3-image-klein",
]
UNGATED = {BFL + "FLUX.2-klein-4B"}
# FLUX 3 uses a modular pipeline: no model_index.json
MODULAR = {BFL + "flux-3-image-klein"}
# other models the tool loads with these, which have no row of their own
EXTRA = {
    BFL + "flux-3-image-klein": [BFL + "flux-3-image-klein-pup"],
}


def modular_files(api, model, token):
    # the diffusers layout only; the FLUX 3 repo also has a single-file copy of
    # the transformer, quantized variants and example images, which the tool
    # does not load
    return [
        f
        for f in api.list_repo_files(model, token=token)
        if f.endswith(".json")
        or ("/" in f and not f.startswith(("assets/", "variants/")))
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Download the models of the Galaxy FLUX tool"
    )
    parser.add_argument(
        "--model-dir", required=True, help="Directory to store the models"
    )
    parser.add_argument(
        "--models",
        nargs="+",
        choices=MODELS,
        default=MODELS,
        help="Models to download (default: all)",
    )
    parser.add_argument(
        "--hf-token",
        # an empty HF_TOKEN counts as none; the Hub rejects an empty Bearer header
        default=os.getenv("HF_TOKEN") or None,
        help="Hugging Face token for the gated models (default: $HF_TOKEN)",
    )
    args = parser.parse_args()
    gated = sorted(set(args.models) - UNGATED)
    if gated and not args.hf_token:
        parser.error(f"set HF_TOKEN or --hf-token (needed for {', '.join(gated)})")

    root = Path(args.model_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    # the tool sets HF_HOME to this directory
    os.environ["HF_HOME"] = str(root)
    os.environ["HF_HUB_DISABLE_TELEMETRY"] = "1"

    from diffusers import DiffusionPipeline, ModularPipeline
    from huggingface_hub import HfApi, snapshot_download

    api = HfApi()
    token = args.hf_token

    def fetch(model):
        if model in MODULAR:
            snapshot_download(
                model, allow_patterns=modular_files(api, model, token), token=token
            )
            # resolves every component from the cache, as the tool does
            ModularPipeline.from_pretrained(model, local_files_only=True)
        elif model.startswith(BFL) and not model.endswith("-pup"):
            DiffusionPipeline.download(model, token=token)
            # fails if a file the pipeline needs is missing from the cache
            DiffusionPipeline.download(model, local_files_only=True)
        else:
            # the upsampler is a plain transformers model
            snapshot_download(model, token=token)

    failed = []
    done = set()
    for model in args.models:
        for repo in [model] + EXTRA.get(model, []):
            if repo in done:
                continue
            print(f"Downloading {repo} ...", flush=True)
            try:
                fetch(repo)
                print(f"OK: {repo}", flush=True)
                done.add(repo)
            except Exception as e:
                print(f"Failed: {repo}: {e}", flush=True)
                failed.append(repo)

    if failed:
        print(f"Failed to download {len(failed)} model(s): {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
