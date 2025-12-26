import os
import sys

import torch
from diffusers import FluxPipeline

model_path = sys.argv[1]

prompt_type = sys.argv[2]
if prompt_type == "file":
    with open(sys.argv[3], "r") as f:
        prompt = f.read().strip()
elif prompt_type == "text":
    prompt = sys.argv[3]
print(f"Creating image from prompt: {prompt}")

if "dev" in model_path:
    num_inference_steps = 20
elif "schnell" in model_path:
    num_inference_steps = 4
else:
    print("Invalid model!")
    sys.exit(1)

snapshots = []
for d in os.listdir(os.path.join(model_path, "snapshots")):
    snapshots.append(os.path.join(model_path, "snapshots", d))
latest_snapshot_path = max(snapshots, key=os.path.getmtime)

pipe = FluxPipeline.from_pretrained(latest_snapshot_path, torch_dtype=torch.bfloat16)
pipe.enable_sequential_cpu_offload()
pipe.vae.enable_slicing()
pipe.vae.enable_tiling()
pipe.to(torch.float16)

image = pipe(
    prompt,
    num_inference_steps=num_inference_steps,
    generator=torch.Generator("cpu").manual_seed(42),
).images[0]

image.save("output.png")
