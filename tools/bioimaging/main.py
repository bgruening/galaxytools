"""
Predict images using AI models from BioImage
"""

import argparse
import numpy as np
from PIL import Image

import torch
import torchvision.transforms as T


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-im", "--imaging_model", required=True, help="Input BioImage model")
    arg_parser.add_argument("-ii", "--image_file", required=True, help="Input image file")
    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["imaging_model"]
    input_image_path = args["image_file"]
    model = torch.load(model_path)
    model.eval()
    test_data = np.load(input_image_path)
    test_data = torch.Tensor(test_data)
    pred_data = model(test_data)
    pred_data_output = pred_data.detach().numpy()
    # save original image matrix
    np.savez("output_predicted_image_matrix.npz", pred_data_output)
    # reshape predicted image matrix to display
    pred_data = torch.squeeze(pred_data)
    if len(pred_data.shape) == 3:
        pred_data = pred_data[0, :, :]
    elif len(pred_data.shape) == 4:
        pred_data = pred_data[0, 0, :, :]
    transform = T.ToPILImage()
    pred_img = transform(pred_data)
    pred_img.save("output_predicted_image.png")
