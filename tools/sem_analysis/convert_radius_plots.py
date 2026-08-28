#!/usr/bin/env python3
"""Convert all DiameterJ radius-plot TIFFs in a result directory to PNG."""
import sys
from pathlib import Path

from PIL import Image


def main() -> None:
    result_dir = Path(sys.argv[1])
    for source in result_dir.glob("*_Radius Plot.tif"):
        destination = source.with_suffix(".png")
        with Image.open(source) as image:
            image.save(destination)
        source.unlink()


if __name__ == "__main__":
    main()
