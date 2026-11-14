"""Generate or edit an image with a FLUX model from the huggingface data table.

The model is loaded by its Hugging Face ID from the cache in HF_HOME, which
the tool points at the `path` of the data table row (HF_HUB_OFFLINE=1).
"""

import argparse
import sys

import torch
from PIL import Image

BFL = "black-forest-labs/"

# One entry per supported model, keyed by the value column of the data table.
# tasks: what the model can do (the form only offers it for these);
# steps, guidance: the defaults of each model card; references: how many
# reference images it takes for editing.
MODELS = {
    BFL
    + "FLUX.2-klein-4B": dict(
        pipeline="flux2_klein",
        tasks={"generate", "edit", "inpaint"},
        steps=4,
        guidance=1.0,
        references=10,
    ),
    BFL
    + "flux-3-image-klein": dict(
        pipeline="flux3",
        tasks={"generate", "edit"},
        steps=16,
        guidance=None,
        references=10,
    ),
}
IMAGES = {"generate": (0, 0), "inpaint": (1, 1)}


def load_pipeline(model_id, spec, task, upsample):
    if spec["pipeline"] == "flux3":
        from diffusers import ComponentsManager, ModularPipeline

        manager = ComponentsManager()
        manager.enable_auto_cpu_offload(device="cuda")
        blocks = ModularPipeline.from_pretrained(model_id).blocks
        if not upsample:
            # the prompt is used as the full caption, as on the model card
            blocks.sub_blocks.pop("prompt_upsampler")
        pipe = blocks.init_pipeline(model_id, components_manager=manager)
        # diffusers 0.42.1 modular pipelines ignore HF_HUB_OFFLINE for sharded
        # checkpoints and ask the Hub for the shard list
        pipe.load_components(
            dtype={"vae": torch.float32, "default": torch.bfloat16},
            local_files_only=True,
        )
        return pipe
    import diffusers

    name = {
        "flux2_klein": "Flux2KleinPipeline",
    }[spec["pipeline"]]
    if task == "inpaint" and spec["pipeline"] == "flux2_klein":
        name = "Flux2KleinInpaintPipeline"
    pipe = getattr(diffusers, name).from_pretrained(
        model_id, torch_dtype=torch.bfloat16
    )
    pipe.enable_model_cpu_offload()
    return pipe


def default_size(images):
    # the size of the first input image, in multiples of 16 within the form's limits
    if not images:
        return 1024, 1024
    return tuple(max(256, min(2048, side // 16 * 16)) for side in images[0].size)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="Hugging Face model ID")
    parser.add_argument(
        "--task", required=True, choices=["generate", "edit", "inpaint"]
    )
    parser.add_argument("--prompt-file", required=True, help="File with the prompt")
    parser.add_argument(
        "--image", action="append", default=[], help="Input image (repeatable)"
    )
    parser.add_argument("--mask", help="inpaint: white where the image is redrawn")
    parser.add_argument("--width", type=int, help="Default: 1024 or the input image's")
    parser.add_argument("--height", type=int, help="Default: 1024 or the input image's")
    parser.add_argument("--steps", type=int, help="Default: the model's")
    parser.add_argument("--guidance", type=float, help="Default: the model's")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--no-upsample", action="store_true", help="FLUX 3: skip the prompt upsampler"
    )
    parser.add_argument("--output", default="output.png")
    args = parser.parse_args()

    spec = MODELS.get(args.model)
    if spec is None:
        sys.exit(f"Unsupported model: {args.model}. Supported: {', '.join(MODELS)}")
    if args.task not in spec["tasks"]:
        sys.exit(f"{args.model} cannot be used for the task '{args.task}'")
    low, high = IMAGES.get(args.task, (1, spec.get("references", 0)))
    if not low <= len(args.image) <= high:
        sys.exit(
            f"The task '{args.task}' with {args.model} takes {low} to {high} images"
        )
    if args.task == "inpaint" and not args.mask:
        sys.exit("The task 'inpaint' needs a mask")
    if not torch.cuda.is_available():
        sys.exit("FLUX needs a CUDA GPU")

    with open(args.prompt_file) as f:
        prompt = f.read().strip()
    if not prompt:
        sys.exit("The prompt is empty")
    print(f"Model: {args.model}\nTask: {args.task}\nPrompt: {prompt}", flush=True)

    images = [Image.open(path).convert("RGB") for path in args.image]
    width, height = default_size(images)
    kwargs = {
        "prompt": prompt,
        "num_inference_steps": args.steps or spec["steps"],
        "generator": torch.Generator("cpu").manual_seed(args.seed),
        "width": args.width or width,
        "height": args.height or height,
    }
    guidance = args.guidance if args.guidance is not None else spec["guidance"]
    if guidance is not None:
        kwargs["guidance_scale"] = guidance
    if args.task == "inpaint":
        kwargs["image"] = images[0]
        kwargs["mask_image"] = Image.open(args.mask).convert("L")
    elif images:
        kwargs["image"] = images if len(images) > 1 else images[0]

    pipe = load_pipeline(args.model, spec, args.task, not args.no_upsample)
    if spec["pipeline"] == "flux3":
        out = pipe(**kwargs, output=["images", "prompt"])
        if not args.no_upsample:
            print(f"Upsampled prompt: {out['prompt']}", flush=True)
        image = out["images"][0]
    else:
        image = pipe(**kwargs).images[0]
    image.save(args.output)
    print(f"Saved {args.output} ({image.width}x{image.height})")


if __name__ == "__main__":
    main()
