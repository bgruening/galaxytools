"""
"""
import os
import subprocess
import argparse
import tensorflow as tf
import tf2onnx
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType
import yaml


sklearn_ensemble_model = "sklearn.ensemble"
tf_model = "tensorflow"


def read_loaded_file(p_loaded_file):
    global_vars = dict()
    locals_vars = dict()
    input_file = yaml.safe_load(p_loaded_file)
    exec(open(input_file).read(), global_vars, locals_vars)
    return global_vars, locals_vars


def save_sklearn_model(obj, output_file):
    initial_type = [('float_input', FloatTensorType([None, 4]))]
    onx = convert_sklearn(obj, initial_types=initial_type)
    with open(output_file, "wb") as f:
        f.write(onx.SerializeToString())


def save_tf_model(obj, output_file):
    curr_path = os.path.abspath(os.getcwd())
    tf_new_path = "{}/{}".format(curr_path, "model")
    if not os.path.exists(tf_new_path):
        os.makedirs(tf_new_path)
    # save model as tf model
    tf.saved_model.save(obj, tf_new_path)
    # OPSET level defines a level of tensorflow operations supported by ONNX
    python_shell_script = "python -m tf2onnx.convert --saved-model " + tf_new_path + " --output " + output_file + " --opset 7 "
    # convert tf/keras model to ONNX and save it to output file
    subprocess.run(python_shell_script, shell=True, check=True)


def check_vars(var_dict):
    if var_dict is not None:
        for key in var_dict:
            obj = var_dict[key]
            obj_class = str(obj.__class__)
            if tf_model in obj_class:
                save_tf_model(obj, output_file)
            elif sklearn_ensemble_model in obj_class: 
                save_sklearn_model(obj, output_file)


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-op", "--output", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    output_file = args["output"]
    # dict of global variables 
    global_vars, locals_vars = read_loaded_file(loaded_file)
    check_vars(locals_vars)
    check_vars(global_vars)
