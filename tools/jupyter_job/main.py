"""
"""

import argparse
import time
import numpy as np
import h5py


def read_loaded_file(p_loaded_file):

    with open(p_loaded_file, "r") as fl:
        py_code = fl.read()
        re_py_code = exec(py_code)
        print(re_py_code)
    return re_py_code
        

if __name__ == "__main__":
    start_time = time.time()

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-op", "--output", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    output_file = args["output"]

    output = read_loaded_file(loaded_file)
    re_output = np.zeros((4, 4))
    print(re_output)

    hf_file = h5py.File(output_file, "w")
    hf_file.create_dataset("zeros", data=re_output)
    hf_file.close()

    end_time = time.time()
    #print("Program finished in %s seconds" % str(end_time - start_time))
