# init_from_local_images.py

import argparse
import json
from fractal_tasks_core.tasks.cellvoyager_to_ome_zarr_init import cellvoyager_to_ome_zarr_init

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--image_dir", required=True)
    parser.add_argument("--allowed_channels", required=True)
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--output_json", required=True)
    parser.add_argument("--num_levels", type=int, default=6)
    parser.add_argument("--coarsening_xy", type=int, default=2)
    parser.add_argument("--image_extension", default="png")
    args = parser.parse_args()

    with open(args.allowed_channels) as f:
        allowed_channels = json.load(f)

    result = cellvoyager_to_ome_zarr_init(
        zarr_urls=[],
        zarr_dir=args.zarr_dir,
        image_dirs=[args.image_dir],
        allowed_channels=allowed_channels,
        num_levels=args.num_levels,
        coarsening_xy=args.coarsening_xy,
        metadata_table_file=None,
        image_extension=args.image_extension
    )

    with open(args.output_json, "w") as f:
        json.dump(result["parallelization_list"], f, indent=2)

