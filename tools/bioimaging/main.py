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
    
    '''if exp_test_data.dim() == 2:
        exp_test_data = exp_test_data.unsqueeze(0)  # Now shape is (1, 1, 512, 512)
        exp_test_data = exp_test_data.unsqueeze(0)
    elif exp_test_data.dim() == 3:
        exp_test_data = exp_test_data.unsqueeze(0)'''
    
    if exp_test_data.shape != target_image_dim:
        for i in range(len(target_image_dim) - exp_test_data.dim()):
            exp_test_data = exp_test_data.unsqueeze(i)
    
        print("Unsqueezed image shape: ", exp_test_data.shape)
        resized_image = F.interpolate(exp_test_data, size=target_image_dim, mode='bilinear', align_corners=False)
        print("Resize image: ", resized_image.shape)
        # Remove the channel dimension if not needed
        exp_test_data = torch.squeeze(resized_image)
        print("Resize and squeezed image: ", exp_test_data.shape)
        #exp_test_data = exp_test_data[slices]
    
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
