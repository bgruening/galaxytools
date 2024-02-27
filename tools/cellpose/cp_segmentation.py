import argparse
import json
import numpy as np
import skimage.io
import matplotlib.pyplot as plt
import os
import warnings

from cellpose import models, plot, transforms


def main(inputs, img_path, img_format, output_dir):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter
    img_path : str
        File path for the input image
    img_format : str
        One of the ['ome.tiff', 'tiff', 'png', 'jpg']
    output_dir : str
        Folder to save the outputs.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    gpu = params['use_gpu']
    model_type = params['model_type']
    chan = params['chan']
    chan2 = params['chan2']
    chan_first = params['chan_first']
    if chan is None:
        channels = None
    else:
        channels = [int(chan), int(chan2) if chan2 is not None else None]

    options = params['options']

    img = skimage.io.imread(img_path)

    print(f"Image shape: {img.shape}")
    # transpose to Ly x Lx x nchann and reshape based on channels
    if img_format.endswith('tiff'):
        img = np.transpose(img, (1, 2, 0))
        img = transforms.reshape(img, channels=channels, chan_first=chan_first)

    print(f"Image shape: {img.shape}")
    model = models.Cellpose(gpu=gpu, model_type=model_type)
    masks, flows, styles, diams = model.eval(img, channels=channels, **options)

    # save masks to tiff
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        skimage.io.imsave(os.path.join(output_dir, 'cp_masks.tif'),
                          masks.astype(np.uint16))

    # make segmentation show #
    #show_seg = False
    #if show_seg:
    if params['show_segmentation']:
        img = skimage.io.imread(img_path)
        # uniform image
        if img_format.endswith('tiff'):
            img = np.transpose(img, (1, 2, 0))
            img = transforms.reshape(img, channels=channels, chan_first=chan_first)

        maski = masks
        flowi = flows[0]
        fig = plt.figure(figsize=(12, 3))
        # can save images (set save_dir=None if not)
        plot.show_segmentation(fig, img, maski, flowi, channels=channels)
        fig.savefig(os.path.join(output_dir, 'segm_show.png'), dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-p", "--img_path", dest="img_path")
    aparser.add_argument("-f", "--img_format", dest="img_format")
    aparser.add_argument("-O", "--output_dir", dest="output_dir")
    args = aparser.parse_args()

    main(args.inputs, args.img_path, args.img_format, args.output_dir)
