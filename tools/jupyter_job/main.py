"""
Train and save machine learning models as ONNX file
"""


import os
import subprocess
import argparse
import tensorflow as tf
import tf2onnx
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType
import yaml


SKLEARN_MODELS = [
    "sklearn.ensemble",
    "sklearn.tree",
    "sklearn.linear_model",
    "sklearn.svm",
    "sklearn.neighbors",
    "sklearn.preprocessing",
    "sklearn.cluster"
]

TF_MODELS = [
    "tensorflow.python.keras.engine.training.Model",
    "tensorflow.python.keras.engine.sequential.Sequential",
    "tensorflow.python.keras.layers"
]


def read_loaded_file(p_loaded_file):
    global_vars = dict()
    locals_vars = dict()
    input_file = yaml.safe_load(p_loaded_file)
    with open(input_file, "r") as f:
        exec(f.read(), global_vars, locals_vars)
        check_vars(locals_vars)
        check_vars(global_vars)


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
            # save tf model
            if len([item for item in TF_MODELS if item in obj_class]) > 0:
                save_tf_model(obj, output_file)
            # save scikit-learn model
            elif len([item for item in SKLEARN_MODELS if item in obj_class]) > 0:
                save_sklearn_model(obj, output_file)


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-op", "--output", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    output_file = args["output"]
    read_loaded_file(loaded_file)
