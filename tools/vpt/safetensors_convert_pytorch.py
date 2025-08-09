#!/usr/bin/env python

import argparse
import os
import sys

import torch

from safetensors import safe_open
from safetensors.torch import save_file


def convert_pt_to_safetensors(input_path, output_path):
    """Convert PyTorch model to SafeTensors format."""
    print(f"Loading PyTorch model from: {input_path}")

    try:
        # Load the PyTorch model
        state_dict = torch.load(input_path)

        print(f"Converting {len(state_dict)} tensors to SafeTensors format...")

        # Save as SafeTensors
        save_file(state_dict, output_path)
        print(f"Successfully converted to SafeTensors: {output_path}")

    except Exception as e:
        print(f"Error converting PyTorch to SafeTensors: {e}")
        sys.exit(1)


def convert_safetensors_to_pt(input_path, output_path):
    """Convert SafeTensors to PyTorch format."""
    print(f"Loading SafeTensors model from: {input_path}")

    try:
        state_dict = {}

        # Load SafeTensors
        with safe_open(input_path, framework="pt", device="cpu") as f:
            for key in f.keys():
                state_dict[key] = f.get_tensor(key)

        print(f"Converting {len(state_dict)} tensors to PyTorch format...")

        # Save as PyTorch
        torch.save(state_dict, output_path)
        print(f"Successfully converted to PyTorch: {output_path}")

    except Exception as e:
        print(f"Error converting SafeTensors to PyTorch: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Convert between PyTorch and SafeTensors",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python converter.py --convert pt_safe -i model.pt -o model.safetensors
  python converter.py --convert safe_pt -i model.safetensors -o model.pt
        """
    )

    parser.add_argument(
        "--convert",
        required=True,
        choices=["pt_safe", "safe_pt"],
        help=(
            "Conversion direction: 'pt_safe' (PyTorch to SafeTensors) "
            "or 'safe_pt' (SafeTensors to PyTorch)"
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input model file path"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output model file path"
    )

    args = parser.parse_args()

    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file does not exist: {args.input}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Perform conversion
    if args.convert == "pt_safe":

        convert_pt_to_safetensors(args.input, args.output)

    elif args.convert == "safe_pt":

        convert_safetensors_to_pt(args.input, args.output)

    print("Conversion completed successfully!")


if __name__ == "__main__":
    main()
