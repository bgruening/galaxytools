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
    #model = torch.load(model_path)
    #model.eval()
    test_data_mat = np.load(args["image_file_npy"])
    test_data_mat = torch.Tensor(test_data_mat)
    print("test data mat:", test_data_mat.shape)
    input_image_shape = args["image_size"]
    print(test_file_size)
    target_shape = input_image_shape.split(",")
    target_shape = [int(i) for i in target_shape]
    print(target_shape)
    target_dimension = 0
    for param in model.named_parameters():
        print(param[1].shape)
        target_dimension = len(param[1].shape)
        break
    current_dimension = len(target_shape)
    for i in range(target_dimension - current_dimension):
        target_shape.append(1)
    re_test_data = torch.tensor(test_data) 
    for i in range(len(target_shape) - len(test_data.shape)):
        re_test_data = torch.unsqueeze(re_test_data, i)
    print(re_test_data.shape)

    reshaped_input = re_test_data


    for i, shp in enumerate(target_shape[::-1]):
        if i == 0:
            reshaped_input = reshaped_input[:shp,:,:,:,:]
        if i == 1:
            reshaped_input = reshaped_input[:,:shp,:,:,:]
        if i == 2:
            reshaped_input = reshaped_input[:,:,:shp,:,:]
        if i == 3:
            reshaped_input = reshaped_input[:,:,:,:shp,:]
        if i == 4:
            reshaped_input = reshaped_input[:,:,:,:,:shp]
    print(reshaped_input.shape)
    #test_data = Image.open(input_image_path)
    #test_data = np.array(test_data)
    #test_data = torch.Tensor(test_data)
    
    #test_data = np.array(test_data)
    #test_data = load_tensor(input_image_path)
    test_data = torch.Tensor(reshaped_input)
    #test_data = torch.reshape(test_data, input_param.shape)
    print("Tiff input:", reshaped_input.shape)
    model = torch.load(model_path)
    model.eval()
    pred_data = model(reshaped_input)
    pred_data_output = pred_data.detach().numpy()
    # save original image matrix
    np.save("output_predicted_image_matrix.npy", pred_data_output)
    # reshape predicted image matrix to display
    pred_data = torch.squeeze(pred_data)
    if len(pred_data.shape) == 3:
        pred_data = pred_data[0, :, :]
    elif len(pred_data.shape) == 4:
        pred_data = pred_data[0, 0, :, :]
    elif len(pred_data.shape) == 5:
        pred_data = pred_data[0, 0, 0, :, :]
    transform = T.ToPILImage()
    pred_img = transform(pred_data)
    pred_img.save("output_predicted_image.png")
