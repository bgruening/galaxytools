import argparse
import json
import warnings

import pandas as pd
from galaxy_ml.utils import to_categorical


def main(inputs, infile, outfile, num_classes=None):
    """
    Parameter
    ---------
    input : str
        File path to galaxy tool parameter

    infile : str
        File paths of input vector

    outfile : str
        File path to output matrix

    num_classes : str
        Total number of classes. If None, this would be inferred as the (largest number in y) + 1

    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    input_header = params['header0']
    header = 'infer' if input_header else None

    input_vector = pd.read_csv(infile, sep='\t', header=header)

    output_matrix = to_categorical(input_vector, num_classes=num_classes)

    print("Input vector shape: %s" % repr(input_vector.shape))
    print("Output matrix shape: %s" % repr(output_matrix.shape))

    output_matrix.to_csv(outfile, sep='\t', header=input_header, index=False)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-y", "--infile", dest="infile")
    aparser.add_argument("-n", "--num_classes", dest="num_classes", type=int, default=None)
    aparser.add_argument("-o", "--outfile", dest="outfile")
    args = aparser.parse_args()

    main(args.inputs, args.infile, args.outfile, args.num_classes)
