import sys

import torch
from diffusers import FluxPipeline

model = sys.argv[1]

prompt_type = sys.argv[2]
if prompt_type == "file":
    with open(sys.argv[3], "r") as f:
        prompt = f.read().strip()
elif prompt_type == "text":
    prompt = sys.argv[3]

if model not in ["black-forest-labs/FLUX.1-dev", "black-forest-labs/FLUX.1-schnell"]:
    print("Invalid model!")
    sys.exit(1)


pipe = FluxPipeline.from_pretrained(model, torch_dtype=torch.bfloat16)
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
