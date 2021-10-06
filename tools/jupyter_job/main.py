import argparse
import os
import subprocess

import h5py
import yaml
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

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
    "tensorflow.python.keras.layers",
    "keras.engine.sequential.Sequential",
    "keras.engine.training.Model",
    "keras.layers"
]

ARRAYS = [
    "numpy.ndarray",
    "list"
]

DATAFRAME = [
    "pandas.core.frame.DataFrame"
]


def read_loaded_file(p_loaded_file, m_file, a_file):
    global_vars = dict()
    input_file = yaml.safe_load(p_loaded_file)
    code_string = open(input_file, "r").read()
    compiled_code = compile(code_string, input_file, 'exec')
    exec(compiled_code, global_vars)
    check_vars(global_vars, m_file, a_file)


def save_sklearn_model(obj, output_file):
    initial_type = [('float_input', FloatTensorType([None, 4]))]
    onx = convert_sklearn(obj, initial_types=initial_type)
    with open(output_file, "wb") as f:
        f.write(onx.SerializeToString())


def save_tf_model(obj, output_file):
    import tensorflow as tf
    curr_path = os.path.abspath(os.getcwd())
    tf_new_path = "{}/{}".format(curr_path, "model")
    if not os.path.exists(tf_new_path):
        os.makedirs(tf_new_path)
    # save model as tf model
    tf.saved_model.save(obj, tf_new_path)
    # OPSET level defines a level of tensorflow operations supported by ONNX
    python_shell_script = "python -m tf2onnx.convert --saved-model " + tf_new_path + " --output " + output_file + " --opset 9 "
    # convert tf/keras model to ONNX and save it to output file
    subprocess.run(python_shell_script, shell=True, check=True)


def save_array(payload, a_file):
    hf_file = h5py.File(a_file, "w")
    for key in payload:
        try:
            hf_file.create_dataset(key, data=payload[key])
        except Exception as e:
            print(e)
            continue
    hf_file.close()


def save_dataframe(payload, a_file):
    for key in payload:
        payload[key].to_hdf(a_file, key=key)


def check_vars(var_dict, m_file, a_file):
    if var_dict is not None:
        arr_payload = dict()
        dataframe_payload = dict()
        for key in var_dict:
            obj = var_dict[key]
            obj_class = str(obj.__class__)
            # save tf model
            if len([item for item in TF_MODELS if item in obj_class]) > 0:
                save_tf_model(obj, m_file)
            # save scikit-learn model
            elif len([item for item in SKLEARN_MODELS if item in obj_class]) > 0:
                save_sklearn_model(obj, m_file)
            # save arrays and lists
            elif len([item for item in ARRAYS if item in obj_class]) > 0:
                if key not in arr_payload:
                    arr_payload[key] = obj
            elif len([item for item in DATAFRAME if item in obj_class]) > 0:
                if key not in dataframe_payload:
                    dataframe_payload[key] = obj
        save_array(arr_payload, a_file)
        save_dataframe(dataframe_payload, a_file)


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-om", "--output_model", required=True, help="")
    arg_parser.add_argument("-oa", "--output_array", required=True, help="")

    # get argument values
    args = vars(arg_parser.parse_args())
    loaded_file = args["loaded_file"]
    model_output_file = args["output_model"]
    array_output_file = args["output_array"]
    read_loaded_file(loaded_file, model_output_file, array_output_file)
