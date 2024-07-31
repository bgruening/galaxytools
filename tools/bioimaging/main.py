"""
Predict images using AI models from BioImage
"""

import argparse

import numpy as np
import PIL
from PIL import Image
import torch
import torchvision.transforms as T

import warnings
from pathlib import Path
from typing import Any, Dict, Literal, Mapping, Optional, Sequence, Tuple, Union

import imageio
from numpy.typing import NDArray
from typing_extensions import assert_never


def find_dim_order(user_in_shape, input_image):    
    image_shape = list(input_image.shape)
    print(user_in_shape, image_shape)
    correct_order = user_in_shape.split(",")[::-1]
    print("correct_order", correct_order)
    correct_order = [int(i) for i in correct_order if i != "1"]
    print("correct order:", correct_order)
    if (correct_order[0] == image_shape[-1]) and (correct_order != image_shape):
        input_image = torch.tensor(input_image.transpose())
    return input_image, correct_order



if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-im", "--imaging_model", required=True, help="Input BioImage model")
    arg_parser.add_argument("-ii", "--image_file", required=True, help="Input image file")
    arg_parser.add_argument("-inp", "--image_file_npy", required=True, help="Input image file as matrix")
    arg_parser.add_argument("-is", "--image_size", required=True, help="Input image file's size")
    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["imaging_model"]
    input_image_path = args["image_file"]
    test_data = imageio.v3.imread(input_image_path, index="...")
    test_data = np.squeeze(test_data)
    test_data = test_data.astype(np.float32)

    test_data_mat = np.load(args["image_file_npy"])
    test_data_mat = torch.Tensor(test_data_mat)
    print("test data mat:", test_data_mat.shape)
    input_image_shape = args["image_size"]
    print(input_image_shape)

    #target_shape = input_image_shape.split(",")
    #target_shape = [i for i in target_shape]
    #print(target_shape)
    # find correct dimensions of input image to be 
    # used by the model
    im_test_data, shape_vals = find_dim_order(input_image_shape, test_data)

    model = torch.load(model_path)
    model.eval()

    target_dimension = 0
    for param in model.named_parameters():
        print(param[1].shape)
        target_dimension = len(param[1].shape)
        break
    current_dimension = len(list(im_test_data.shape))

    slices = tuple(slice(0, s_val) for s_val in shape_vals)
    # Apply the slices to the reshaped_input
    im_test_data = im_test_data[slices]
    exp_test_data = torch.tensor(im_test_data) 

    for i in range(target_dimension - current_dimension):
        exp_test_data = torch.unsqueeze(exp_test_data, i)
    print(exp_test_data.shape)


    #reshaped_input = re_test_data
    #s_vals = target_shape

    #slices = tuple(slice(0, s_val) for s_val in s_vals)
    # Apply the slices to the reshaped_input
    #reshaped_input = reshaped_input[slices]
    #if len(reshaped_input.shape) == 3:
    #    reshaped_input = reshaped_input[:s_vals[0], :s_vals[1], : s_vals[2]]
    #if len(reshaped_input.shape) == 4:
    #    reshaped_input = reshaped_input[:s_vals[0], :s_vals[1], : s_vals[2], : s_vals[3]]
    #if len(reshaped_input.shape) == 5:
    #    reshaped_input = reshaped_input[:s_vals[0], :s_vals[1], : s_vals[2], :s_vals[3], : s_vals[4]]
    #print(reshaped_input.shape)
    #test_data = Image.open(input_image_path)
    #test_data = np.array(test_data)
    #test_data = torch.Tensor(test_data)
    
    #test_data = np.array(test_data)
    #test_data = load_tensor(input_image_path)
    #test_data = torch.Tensor(reshaped_input)
    #test_data = torch.reshape(test_data, input_param.shape)
    #print("Tiff input:", reshaped_input.shape)
    #reshaped_input = torch.Tensor(reshaped_input.astype(np.float32))
    #model = torch.load(model_path)
    #model.eval()
    pred_data = model(exp_test_data)
    pred_data_output = pred_data.detach().numpy()
    # save original image matrix
    np.save("output_predicted_image_matrix.npy", pred_data_output)
    # reshape predicted image matrix to display
    pred_data = torch.squeeze(pred_data)
    pred_numpy = pred_data.detach().numpy()
    pred_numpy = pred_numpy * 255
    pred_numpy = pred_numpy.astype(np.uint8)
    imageio.v3.imwrite("output_predicted_image.tif", pred_numpy, extension=".tif")
