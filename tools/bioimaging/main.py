"""
Predict images using AI models from BioImage.IO
"""

import argparse

import imageio
import numpy as np
import torch
import torch.nn.functional as F


def find_dim_order(user_in_shape, input_image):
    """
    Find the correct order of input image's
    shape. For a few models, the order of input size
    mentioned in the RDF.yaml file is reversed compared
    to the input image's original size. If it is reversed,
    transpose the image to find correct order of image's
    dimensions.
    """
    image_shape = list(input_image.shape)
    # reverse the input shape provided from RDF.yaml file
    correct_order = user_in_shape.split(",")[::-1]
    # remove 1s from the original dimensions
    correct_order = [int(i) for i in correct_order if i != "1"]
    if (correct_order[0] == image_shape[-1]) and (correct_order != image_shape):
        input_image = torch.tensor(input_image.transpose())
    return input_image, correct_order


def dynamic_resize(image: torch.Tensor, target_shape: tuple):
    """
    Resize an input tensor dynamically to the target shape.
    
    Parameters:
    - image: Input tensor with shape (C, D1, D2, ..., DN) (any number of spatial dims)
    - target_shape: Tuple specifying the target shape (C', D1', D2', ..., DN')
    
    Returns:
    - Resized tensor with shape target_shape.
    """
    # Extract input shape
    input_shape = image.shape
    num_dims = len(input_shape)  # Includes channels and spatial dimensions

    # Ensure target shape matches the number of dimensions
    if len(target_shape) != num_dims:
        raise ValueError(f"Target shape {target_shape} must match input dimensions {num_dims}")

    # Extract target channels and spatial sizes
    target_channels = target_shape[0]  # First element is the target channel count
    target_spatial_size = target_shape[1:]  # Remaining elements are spatial dimensions

    # Step 1: Expand channels dynamically if needed
    if target_channels > input_shape[0]:
        # Expand existing channels to match target_channels
        image = image.expand(target_channels, *input_shape[1:])
    elif target_channels < input_shape[0]:
        # Reduce channels using interpolation
        image = image.unsqueeze(0)  # Add batch dim (1, C, ...)
        image = F.interpolate(image, size=(target_channels, *input_shape[1:]), mode='trilinear' if num_dims == 4 else 'bilinear', align_corners=False)
        image = image.squeeze(0)  # Remove batch dim

    # Step 2: Add batch dim (N=1) for resizing
    image = image.unsqueeze(0)  # Shape: (1, C, D1, D2, ..., DN)

    # Step 3: Resize spatial dimensions dynamically
    image = F.interpolate(image, size=target_spatial_size, mode='trilinear' if num_dims == 4 else 'bilinear', align_corners=False)

    # Step 4: Remove batch dim
    image = image.squeeze(0)  # Shape: (C, D1', D2', ..., DN')

    return image


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-im", "--imaging_model", required=True, help="Input BioImage model")
    arg_parser.add_argument("-ii", "--image_file", required=True, help="Input image file")
    arg_parser.add_argument("-is", "--image_size", required=True, help="Input image file's size")
    arg_parser.add_argument("-ia", "--image_axes", required=True, help="Input image file's axes")

    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["imaging_model"]
    input_image_path = args["image_file"]
    input_size = args["image_size"]
    
    # load all embedded images in TIF file
    test_data = imageio.v3.imread(input_image_path, index="...")
    #test_data = np.squeeze(test_data)
    test_data = test_data.astype(np.float32)
    print("Input image shape:", test_data.shape)
    test_data = np.squeeze(test_data)
    
    #exp_test_data, shape_vals = find_dim_order(input_size, test_data)
    #print(shape_vals)
    target_image_dim = input_size.split(",")[::-1]
    target_image_dim = [int(i) for i in target_image_dim if i != "1"]
    target_image_dim = tuple(target_image_dim)
    print("target_image_dim:", target_image_dim)
    # apply the slices to the reshaped_input
    exp_test_data = torch.tensor(test_data)
    
    reversed_order = list(reversed(range(exp_test_data.dim())))
    print(exp_test_data.permute(*reversed_order).shape)
    exp_test_data_T = exp_test_data.permute(*reversed_order)
    print(exp_test_data_T.shape)
    if exp_test_data_T.shape == target_image_dim:
        exp_test_data = exp_test_data_T
    print("Transposed image dim: ", exp_test_data.shape)
    
    if exp_test_data.shape != target_image_dim:
        for i in range(len(target_image_dim) - exp_test_data.dim()):
            exp_test_data = exp_test_data.unsqueeze(i)   
        exp_test_data = dynamic_resize(exp_test_data, target_image_dim)
        print("Resized image: ", exp_test_data.shape)
    
    current_dimension = len(exp_test_data.shape)
    input_axes = args["image_axes"]
    target_dimension = len(input_axes)
    
    print(current_dimension, target_dimension)
    print("Before expansion: ", exp_test_data.shape)
    
    for i in range(target_dimension - current_dimension):
        exp_test_data = torch.unsqueeze(exp_test_data, i)

    print(exp_test_data.shape)
    print("After expansion: ", exp_test_data.shape)
    
    # load model
    model = torch.load(model_path)
    model.eval()
    
    #for param in model.named_parameters():
    #    print("Num parameters: ", len(param[1].shape))
    
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
