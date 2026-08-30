import argparse
import json
import os
import warnings

# Make sure, that `MKL_NUM_THREADS` is set to 1, to ensure reproducibility:
# https://forum.image.sc/t/reproducibility-how-we-spent-years-building-containers-and-then-mkl-decided-to-screw-our-results/109599
if str(os.environ['MKL_NUM_THREADS']) != '1':
    warnings.warn('MKL_NUM_THREADS must be set to 1 to ensure reproducibility, and will be adjusted accordingly for now.')
    os.environ['MKL_NUM_THREADS'] = '1'

# Load the remaining packages *after* adjusting `MKL_NUM_THREADS` (this likely necessary for it to take effect)
import matplotlib.pyplot as plt
import numpy as np
import skimage.io
import torch
from cellpose import models, plot

# Apply PyTorch guidelines for reproducibility
torch.backends.cudnn.benchmark = True
torch.backends.cudnn.deterministic = True
torch.manual_seed(0)


def main(inputs, img_path, output_dir):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter
    img_path : str
        File path for the input image
    output_dir : str
        Folder to save the outputs.
    """
    warnings.simplefilter('ignore')
    np.random.seed(42)
    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    gpu = params['use_gpu']
    model_type = params['model_type']
    options = params['options']
    img = skimage.io.imread(img_path)

    print(f"Image shape: {img.shape}")

    model = models.Cellpose(gpu=gpu, model_type=model_type)
    masks, flows, styles, diams = model.eval(img, channels=[0, 0], **options)

    # save masks to tiff
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        skimage.io.imsave(os.path.join(output_dir, 'cp_masks.tif'),
                          masks.astype(np.uint16))

    # make segmentation show #
    if params['show_segmentation']:
        img = skimage.io.imread(img_path)

        maski = masks
        flowi = flows[0]
        fig = plt.figure(figsize=(8, 2))
        # can save images (set save_dir=None if not)
        plot.show_segmentation(fig, img, maski, flowi, channels=[0, 0])
        fig.savefig(os.path.join(output_dir, 'segm_show.png'), dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-p", "--img_path", dest="img_path")
    aparser.add_argument("-O", "--output_dir", dest="output_dir")
    args = aparser.parse_args()

    main(args.inputs, args.img_path, args.output_dir)
