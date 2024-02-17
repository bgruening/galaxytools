import argparse
import json
import warnings

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder


def main(inputs, infile, outfile):
    """
    Parameter
    ---------
    input : str
        File path to galaxy tool parameter

    infile : str
        File paths of input vector

    outfile : str
        File path to output vector

    """
    warnings.simplefilter("ignore")

    with open(inputs, "r") as param_handler:
        params = json.load(param_handler)

    input_header = params["header0"]
    header = "infer" if input_header else None

    input_vector = pd.read_csv(infile, sep="\t", header=header)

    le = LabelEncoder()

    output_vector = le.fit_transform(input_vector)

    np.savetxt(outfile, output_vector, fmt="%d", delimiter="\t")


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-y", "--infile", dest="infile")
    aparser.add_argument("-o", "--outfile", dest="outfile")
    args = aparser.parse_args()

    main(args.inputs, args.infile, args.outfile)
