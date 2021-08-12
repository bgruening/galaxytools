"""
"""

import argparse
import numpy as np
import h5py
import tensorflow as tf


def read_loaded_file(p_loaded_file):
    dynamic_vars = dict()
    exec(open(p_loaded_file).read(), dynamic_vars)
    return dynamic_vars 
        

if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-op", "--output", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    output_file = args["output"]

    dynamic_vars = read_loaded_file(loaded_file)
    if dynamic_vars is not None:
        weights = dynamic_vars["weights"]
        hf_file = h5py.File(output_file, "w")
        for i in range(len(weights)):
            d_name = "layer_{}".format(str(i))
            hf_file.create_dataset(d_name, data=weights[i])
        hf_file.close()
    
