"""
"""
import os
import subprocess
import argparse
import tensorflow as tf
import tf2onnx
import yaml


def read_loaded_file(p_loaded_file):
    global_vars = dict()
    input_file = yaml.safe_load(p_loaded_file)
    exec(open(input_file).read(), global_vars)
    return global_vars


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-op", "--output", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    output_file = args["output"]

    # dict of global variables 
    global_vars = read_loaded_file(loaded_file)

    # save model
    if global_vars is not None:
        trained_model = global_vars["model"]
        model_type = str(type(trained_model))
        if "tensorflow" in model_type or "keras" in model_type:
            curr_path = os.path.abspath(os.getcwd())
            tf_new_path = "{}/{}".format(curr_path, "model")
            if not os.path.exists(tf_new_path):
                os.makedirs(tf_new_path)
            # save model as tf model
            tf.saved_model.save(trained_model, tf_new_path)
            # OPSET level defines a level of tensorflow operations supported by ONNX
            python_shell_script = "python -m tf2onnx.convert --saved-model " + tf_new_path + " --output " + output_file + " --opset 7 "
            # convert tf/keras model to ONNX and save it to output file
            subprocess.run(python_shell_script, shell=True, check=True)
