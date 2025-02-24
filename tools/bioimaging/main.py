"""
Predict images using AI models from BioImage.IO
"""

import argparse

import imageio
import numpy as np
import torch
import torch.nn.functional as F


def dynamic_resize(image: torch.Tensor, target_shape: tuple):
    """
    Resize an input tensor dynamically to the target shape.

    Parameters:
    - image: Input tensor with shape (C, D1, D2, ..., DN) (any number of spatial dims)
    - target_shape: Tuple specifying the target shape (C', D1', D2', ..., DN')

    Returns:
    - Resized tensor with target shape target_shape.
    """
    # Extract input shape
    input_shape = image.shape
    num_dims = len(input_shape)  # Includes channels and spatial dimensions

    # Ensure target shape matches the number of dimensions
    if len(target_shape) != num_dims:
        raise ValueError(
            f"Target shape {target_shape} must match input dimensions {num_dims}"
        )

    # Extract target channels and spatial sizes
    target_channels = target_shape[0]  # First element is the target channel count
    target_spatial_size = target_shape[1:]  # Remaining elements are spatial dimensions

    # Add batch dim (N=1) for resizing
    image = image.unsqueeze(0)

    # Choose the best interpolation mode based on dimensionality
    if num_dims == 4:
        interp_mode = "trilinear"
    elif num_dims == 3:
        interp_mode = "bilinear"
    elif num_dims == 2:
        interp_mode = "bicubic"
    else:
        interp_mode = "nearest"

    # Resize spatial dimensions dynamically
    image = F.interpolate(
        image, size=target_spatial_size, mode=interp_mode, align_corners=False
    )

    # Adjust channels if necessary
    current_channels = image.shape[1]

    if target_channels > current_channels:
        # Expand channels by repeating existing ones
        expand_factor = target_channels // current_channels
        remainder = target_channels % current_channels
        image = image.repeat(1, expand_factor, *[1] * (num_dims - 1))

        if remainder > 0:
            extra_channels = image[
                :, :remainder, ...
            ]  # Take the first few channels to match target
            image = torch.cat([image, extra_channels], dim=1)

    elif target_channels < current_channels:
        # Reduce channels by averaging adjacent ones
        image = image[:, :target_channels, ...]  # Simply slice to reduce channels
    return image.squeeze(0)  # Remove batch dimension before returning


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-im", "--imaging_model", required=True, help="Input BioImage model"
    )
    arg_parser.add_argument(
        "-ii", "--image_file", required=True, help="Input image file"
    )
    arg_parser.add_argument(
        "-is", "--image_size", required=True, help="Input image file's size"
    )
    arg_parser.add_argument(
        "-ia", "--image_axes", required=True, help="Input image file's axes"
    )

    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["imaging_model"]
    input_image_path = args["image_file"]
    input_size = args["image_size"]

    # load all embedded images in TIF file
    test_data = imageio.v3.imread(input_image_path, index="...")
    test_data = test_data.astype(np.float32)
    test_data = np.squeeze(test_data)

    target_image_dim = input_size.split(",")[::-1]
    target_image_dim = [int(i) for i in target_image_dim if i != "1"]
    target_image_dim = tuple(target_image_dim)

    exp_test_data = torch.tensor(test_data)
    # check if image dimensions are reversed
    reversed_order = list(reversed(range(exp_test_data.dim())))
    exp_test_data_T = exp_test_data.permute(*reversed_order)
    if exp_test_data_T.shape == target_image_dim:
        exp_test_data = exp_test_data_T
    if exp_test_data.shape != target_image_dim:
        for i in range(len(target_image_dim) - exp_test_data.dim()):
            exp_test_data = exp_test_data.unsqueeze(i)
        try:
            exp_test_data = dynamic_resize(exp_test_data, target_image_dim)
        except Exception as e:
            raise RuntimeError(f"Error during resizing: {e}") from e

    current_dimension = len(exp_test_data.shape)
    input_axes = args["image_axes"]
    target_dimension = len(input_axes)
    # expand input image based on the number of target dimensions
    for i in range(target_dimension - current_dimension):
        exp_test_data = torch.unsqueeze(exp_test_data, i)

    # load model
    model = torch.load(model_path)
    model.eval()

    # make prediction
    pred_data = model(exp_test_data)
    pred_data_output = pred_data.detach().numpy()

    # save original image matrix
    np.save("output_predicted_image_matrix.npy", pred_data_output)

    # post process predicted file to correctly save as TIF file
    pred_data = torch.squeeze(pred_data)
    pred_numpy = pred_data.detach().numpy()

    # write predicted TIF image to file
    imageio.v3.imwrite("output_predicted_image.tif", pred_numpy, extension=".tif")
