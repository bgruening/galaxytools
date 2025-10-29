import argparse
import json
import warnings

import numpy as np
import pandas as pd
from galaxy_ml.model_persist import dump_model_to_h5, load_model_from_h5
from sklearn.preprocessing import LabelEncoder


def main(inputs, infile, outfile, encoder_outfile=None, model=None):
    """
    Parameter
    ---------
    input : str
        File path to galaxy tool parameter

    infile : str
        File paths of input vector

    outfile : str
        File path to output vector

    encoder_outfile : str
        File path to encoder hdf5 output

    model : str
        File path to prefitted model
    """

    warnings.simplefilter("ignore")

    with open(inputs, "r") as param_handler:
        params = json.load(param_handler)

    input_header = params["header0"]
    header = "infer" if input_header else None

    input_vector = pd.read_csv(infile, sep="\t", header=header)

    if model:
        le = load_model_from_h5(model)
        output_vector = le.fit_transform(input_vector)

    else:
        le = LabelEncoder()
        output_vector = le.fit_transform(input_vector)

    np.savetxt(outfile, output_vector, fmt="%d", delimiter="\t")

    if encoder_outfile:
        dump_model_to_h5(le, encoder_outfile)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-y", "--infile", dest="infile")
    aparser.add_argument("-o", "--outfile", dest="outfile")
    aparser.add_argument("--encoder_outfile", dest="encoder_outfile")
    aparser.add_argument("--model", dest="model")
    args = aparser.parse_args()

    main(args.inputs, args.infile, args.outfile,
         args.encoder_outfile, args.model)
