"""
Predict images using AI models from BioImage
"""

import argparse

import cv2
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

#from bioimageio.core.common import Axis, Tensor
from bioimageio.core.tensor import Tensor
from bioimageio.core.axis import Axis
from bioimageio.spec.model import v0_4
from bioimageio.spec.model.v0_4 import InputTensorDescr as InputTensorDescr04
from bioimageio.spec.model.v0_4 import OutputTensorDescr as OutputTensorDescr04
from bioimageio.spec.model.v0_5 import (
    AnyAxis,
    AxisId,
    BatchAxis,
    ChannelAxis,
    Identifier,
    InputTensorDescr,
    OutputTensorDescr,
    SpaceInputAxis,
    convert_axes,
)
from bioimageio.spec.utils import load_array


def interprete_array(
    nd_array: NDArray[Any],
    n_expected_space_axes: Optional[int] = None,
) -> Tensor:

    ndim = nd_array.ndim
    if ndim == 2 and (n_expected_space_axes is None or n_expected_space_axes >= 2):
        current_axes = (
            SpaceInputAxis(id=AxisId("y"), size=nd_array.shape[0]),
            SpaceInputAxis(id=AxisId("x"), size=nd_array.shape[1]),
        )
    elif ndim == 3 and (
        (n_expected_space_axes is None and any(s <= 3 for s in nd_array.shape))
        or n_expected_space_axes == 2
    ):
        current_axes = (
            ChannelAxis(
                channel_names=[
                    Identifier(f"channel{i}") for i in range(nd_array.shape[0])
                ]
            ),
            SpaceInputAxis(id=AxisId("y"), size=nd_array.shape[1]),
            SpaceInputAxis(id=AxisId("x"), size=nd_array.shape[2]),
        )
    elif ndim == 3 and (n_expected_space_axes is None or n_expected_space_axes == 3):
        current_axes = (
            SpaceInputAxis(id=AxisId("z"), size=nd_array.shape[0]),
            SpaceInputAxis(id=AxisId("y"), size=nd_array.shape[1]),
            SpaceInputAxis(id=AxisId("x"), size=nd_array.shape[2]),
        )
    elif ndim == 4:
        current_axes = (
            ChannelAxis(
                channel_names=[
                    Identifier(f"channel{i}") for i in range(nd_array.shape[0])
                ]
            ),
            SpaceInputAxis(id=AxisId("z"), size=nd_array.shape[1]),
            SpaceInputAxis(id=AxisId("y"), size=nd_array.shape[2]),
            SpaceInputAxis(id=AxisId("x"), size=nd_array.shape[3]),
        )
    elif ndim == 5:
        current_axes = (
            BatchAxis(),
            ChannelAxis(
                channel_names=[
                    Identifier(f"channel{i}") for i in range(nd_array.shape[1])
                ]
            ),
            SpaceInputAxis(id=AxisId("z"), size=nd_array.shape[2]),
            SpaceInputAxis(id=AxisId("y"), size=nd_array.shape[3]),
            SpaceInputAxis(id=AxisId("x"), size=nd_array.shape[4]),
        )
    else:
        raise ValueError(
            f"Could not guess an axis mapping for {nd_array.shape} with {n_expected_space_axes} expected space axes"
        )
    print("current_axes", current_axes)

    current_axes_ids = tuple(AxisId(str(a.id)) for a in current_axes)
    print("current_axes_ids", current_axes_ids)

    return Tensor(nd_array, dims=current_axes_ids)

def load_tensor(
    path: Path,
    axes: Optional[Sequence[Axis]] = None,
) -> Tensor:

    #ext = path.suffix
    #if ext == ".npy":
    #    array = load_array(path)
    #else:
    is_volume = True if axes is None else sum(a.type != "channel" for a in axes) > 2
    array = imageio.volread(path) if is_volume else imageio.imread(path)

    if axes is None:
        return interprete_array(array)
    else:
        return Tensor(array, dims=tuple(a.id for a in axes))


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-im", "--imaging_model", required=True, help="Input BioImage model")
    arg_parser.add_argument("-ii", "--image_file", required=True, help="Input image file")
    arg_parser.add_argument("-inp", "--image_file_npy", required=True, help="Input image file as matrix")
    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["imaging_model"]
    input_image_path = args["image_file"]
    #model = torch.load(model_path)
    #model.eval()
    test_data_mat = np.load(args["image_file_npy"])
    test_data_mat = torch.Tensor(test_data_mat)
    print("test data mat:", test_data_mat.shape)
    print()
    #input_param = model.parameters()
    #for para in input_param:
    #    print("First layer shape:", para.size())
    #    break
    #    break
    #print("", input_param.shape)
    #for layer in model.modules():
        #print(layer)
    #test_data = Image.open(input_image_path)
    #test_data = np.array(test_data)
    #test_data = torch.Tensor(test_data)
    test_data = cv2.imread(input_image_path, -1)
    #test_data = np.array(test_data)
    #test_data = load_tensor(input_image_path)
    test_data = torch.Tensor(test_data)
    #test_data = torch.reshape(test_data, input_param.shape)
    print("Tiff input:", test_data.shape)
    model = torch.load(model_path)
    model.eval()
    pred_data = model(test_data_mat)
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
