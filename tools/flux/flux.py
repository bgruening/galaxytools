import os
import sys


model = sys.argv[1]
prompt = sys.argv[2]

os.environ["HF_HOME"] = "/scratch/users/arash/hf_home"

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
