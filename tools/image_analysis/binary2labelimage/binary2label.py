import argparse
import sys
import skimage.io
from skimage.measure import label
import numpy as np
import warnings
from PIL import Image

parser = argparse.ArgumentParser()
parser.add_argument('input_file', type=argparse.FileType('r'), default=sys.stdin, help='input file')
parser.add_argument('out_file', type=argparse.FileType('w'), default=sys.stdin, help='out file (TIFF)')
args = parser.parse_args()

img_in = skimage.io.imread(args.input_file.name) > 0
res = label(img_in).astype(np.int32)

res = Image.fromarray(res)
res.save(args.out_file.name, "tiff")