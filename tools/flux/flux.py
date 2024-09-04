import sys
import os

model = sys.argv[1]
prompt = sys.argv[2]
if model == "black-forest-labs/FLUX.1-dev":
    with open(sys.argv[3], "r") as f:
        hf_token = f.read().strip()
    if not hf_token:
        print("HUGGINGFACE HUB TOKEN is not provided in user preferences!")
        sys.exit(1)
    os.environ["HUGGINGFACE_HUB_TOKEN"] = hf_token


import torch
from diffusers import FluxPipeline


pipe = FluxPipeline.from_pretrained(
    model,
    torch_dtype=torch.bfloat16,
)
pipe.enable_sequential_cpu_offload()
pipe.vae.enable_slicing()
pipe.vae.enable_tiling()
pipe.to(torch.float16)

image = pipe(
    prompt,
    num_inference_steps=4,
    generator=torch.Generator("cpu").manual_seed(42),
).images[0]

image.save("output.png")
