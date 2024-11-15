import argparse
import os
import subprocess
import warnings
from zipfile import ZipFile

import h5py
import yaml
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType


warnings.filterwarnings("ignore")
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
    "tensorflow.python.keras.engine.functional.Functional",
    "tensorflow.python.keras.layers",
    "keras.engine.functional.Functional",
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

SCALAR_TYPES = [
    "int",
    "float",
    "str"
]


def find_replace_paths(script_file, updated_data_dict):
    for item in updated_data_dict:
        g_path = updated_data_dict[item]
        script_file = script_file.replace(item, g_path)
    return script_file


def update_ml_files_paths(old_file_paths, new_file_paths):
    if old_file_paths == "" or old_file_paths is None or new_file_paths == "" or new_file_paths is None:
        return dict()
    o_files = old_file_paths.split(",")
    n_files = new_file_paths.split(",")
    new_paths_dict = dict()
    for i, o_f in enumerate(o_files):
        new_paths_dict[o_f] = n_files[i]
    return new_paths_dict


def read_loaded_file(new_paths_dict, p_loaded_file, a_file, w_dir, z_file):
    global_vars = dict()
    input_file = yaml.safe_load(p_loaded_file)
    code_string = open(input_file, "r").read()
    re_code_string = find_replace_paths(code_string, new_paths_dict)
    compiled_code = compile(re_code_string, input_file, 'exec')
    exec(compiled_code, global_vars)
    check_vars(w_dir, global_vars, a_file)
    zip_files(w_dir, z_file)


def zip_files(w_dir, z_file):
    with ZipFile(z_file, 'w') as zip_file:
        for f_path in os.listdir(w_dir):
            zip_file.write(f_path)


def create_model_path(curr_path, key):
    onnx_path = curr_path + "/model_outputs"
    if not os.path.exists(onnx_path):
        os.makedirs(onnx_path)
    onnx_model_path = curr_path + "/model_outputs/" + "onnx_model_{}.onnx".format(key)
    return onnx_model_path


def save_sklearn_model(w_dir, key, obj):
    initial_type = [('float_input', FloatTensorType([None, 4]))]
    onx = convert_sklearn(obj, initial_types=initial_type)
    sk_model_path = create_model_path(w_dir, key)
    with open(sk_model_path, "wb") as f:
        f.write(onx.SerializeToString())


def save_tf_model(w_dir, key, obj):
    import tensorflow as tf
    tf_file_key = "tf_model_{}".format(key)
    tf_model_path = "{}/{}".format(w_dir, tf_file_key)
    if not os.path.exists(tf_model_path):
        os.makedirs(tf_model_path)
    # save model as tf model
    tf.saved_model.save(obj, tf_model_path)
    # save model as ONNX
    tf_onnx_model_p = create_model_path(w_dir, key)
    # OPSET level defines a level of tensorflow operations supported by ONNX
    python_shell_script = "python -m tf2onnx.convert --saved-model " + tf_model_path + " --output " + tf_onnx_model_p + " --opset 15 "
    # convert tf/keras model to ONNX and save it to output file
    subprocess.run(python_shell_script, shell=True, check=True)


def save_primitives(payload, a_file):
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


def check_vars(w_dir, var_dict, a_file):
    if var_dict is not None:
        primitive_payload = dict()
        dataframe_payload = dict()
        for key in var_dict:
            obj = var_dict[key]
            obj_class = str(obj.__class__)
            # save tf model
            if len([item for item in TF_MODELS if item in obj_class]) > 0:
                save_tf_model(w_dir, key, obj)
            # save scikit-learn model
            elif len([item for item in SKLEARN_MODELS if item in obj_class]) > 0:
                save_sklearn_model(w_dir, key, obj)
            # save arrays and lists
            elif len([item for item in ARRAYS if item in obj_class]) > 0:
                if key not in primitive_payload:
                    primitive_payload[key] = obj
            elif len([item for item in DATAFRAME if item in obj_class]) > 0:
                if key not in dataframe_payload:
                    dataframe_payload[key] = obj
            elif len([item for item in SCALAR_TYPES if item in obj_class]) > 0:
                if key not in primitive_payload:
                    primitive_payload[key] = obj
        save_primitives(primitive_payload, a_file)
        save_dataframe(dataframe_payload, a_file)


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-mlp", "--ml_paths", required=True, help="")
    arg_parser.add_argument("-ldf", "--loaded_file", required=True, help="")
    arg_parser.add_argument("-wd", "--working_dir", required=True, help="")
    arg_parser.add_argument("-oz", "--output_zip", required=True, help="")
    arg_parser.add_argument("-oa", "--output_array", required=True, help="")
    arg_parser.add_argument("-mlf", "--ml_h5_files", required=True, help="")
    # get argument values
    args = vars(arg_parser.parse_args())
    ml_paths = args["ml_paths"]
    loaded_file = args["loaded_file"]
    array_output_file = args["output_array"]
    zip_output_file = args["output_zip"]
    working_dir = args["working_dir"]
    ml_h5_files = args["ml_h5_files"]
    new_paths_dict = update_ml_files_paths(ml_paths, ml_h5_files)
    read_loaded_file(new_paths_dict, loaded_file, array_output_file, working_dir, zip_output_file)
