"""
Predict images using AI models from BioImage.IO
"""

import argparse

import imageio
import numpy as np
import torch


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

    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["imaging_model"]
    input_image_path = args["image_file"]
    
    # load all embedded images in TIF file
    test_data = imageio.v3.imread(input_image_path, index="...")
    #test_data = np.squeeze(test_data)
    test_data = test_data.astype(np.float32)
    current_dimension = len(test_data.shape)

    # assess the correct dimensions of TIF input image
    input_image_shape = args["image_size"]
    target_dimension = len(input_image_shape)
    print(current_dimension, target_dimension)
    #im_test_data, shape_vals = find_dim_order(input_image_shape, test_data)
    #print(im_test_data.shape, shape_vals)
    # expand input image's dimensions
    exp_test_data = torch.tensor(test_data)
    print(exp_test_data.shape)
    for i in range(target_dimension - current_dimension):
        exp_test_data = torch.unsqueeze(exp_test_data, i)
    print(exp_test_data.shape)
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
