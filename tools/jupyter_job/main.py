"""
"""

import argparse
import numpy as np
import h5py
import tensorflow as tf


def read_loaded_file(p_loaded_file):
    re_py_code = exec(open(p_loaded_file).read())
    return re_py_code
        

if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-op", "--output", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    output_file = args["output"]

    output = read_loaded_file(loaded_file)
    print(output)
    re_output = np.zeros((4, 4))
    print(re_output)

    hf_file = h5py.File(output_file, "w")
    hf_file.create_dataset("zeros", data=re_output)
    hf_file.close()
