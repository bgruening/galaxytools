import os
import numpy as np
import json
import h5py

from keras.models import model_from_json, Sequential
from keras.layers import Dense, GRU, Dropout
from keras.layers.embeddings import Embedding
from keras.layers.core import SpatialDropout1D
from keras.optimizers import RMSprop
from keras import backend as K


def read_file(file_path):
    """
    Read a file
    """
    with open(file_path, "r") as json_file:
        file_content = json.loads(json_file.read())
    return file_content


def write_file(file_path, content):
    """
    Write a file
    """
    remove_file(file_path)
    with open(file_path, "w") as json_file:
        json_file.write(json.dumps(content))


def save_processed_workflows(file_path, unique_paths):
    workflow_paths_unique = ""
    for path in unique_paths:
        workflow_paths_unique += path + "\n"
    with open(file_path, "w") as workflows_file:
        workflows_file.write(workflow_paths_unique)


def load_saved_model(model_config, model_weights):
    """
    Load the saved trained model using the saved network and its weights
    """
    # load the network
    loaded_model = model_from_json(model_config)
    # load the saved weights into the model
    loaded_model.set_weights(model_weights)
    return loaded_model


def format_tool_id(tool_link):
    """
    Extract tool id from tool link
    """
    tool_id_split = tool_link.split("/")
    tool_id = tool_id_split[-2] if len(tool_id_split) > 1 else tool_link
    return tool_id


def get_HDF5(hf, d_key):
    """
    Read h5 file to get train and test data
    """
    return hf.get(d_key).value


def save_HDF5(hf_file, d_key, data, d_type=""):
    """
    Save datasets as h5 file
    """
    if (d_type == 'json'):
        data = json.dumps(data)
    hf_file.create_dataset(d_key, data=data)


def set_trained_model(dump_file, model_values):
    """
    Create an h5 file with the trained weights and associated dicts
    """
    hf_file = h5py.File(dump_file, 'w')
    for key in model_values:
        value = model_values[key]
        if key == 'model_weights':
            for idx, item in enumerate(value):
                w_key = "weight_" + str(idx)
                if w_key in hf_file:
                    hf_file.modify(w_key, item)
                else:
                    hf_file.create_dataset(w_key, data=item)
        else:
            if key in hf_file:
                hf_file.modify(key, json.dumps(value))
            else:
                hf_file.create_dataset(key, data=json.dumps(value))
    hf_file.close()


def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)


def extract_configuration(config_object):
    config_loss = dict()
    for index, item in enumerate(config_object):
        config_loss[index] = list()
        d_config = dict()
        d_config['loss'] = item['result']['loss']
        d_config['params_config'] = item['misc']['vals']
        config_loss[index].append(d_config)
    return config_loss


def get_best_parameters(mdl_dict):
    """
    Get param values (defaults as well)
    """
    lr = float(mdl_dict.get("learning_rate", "0.001"))
    embedding_size = int(mdl_dict.get("embedding_size", "512"))
    dropout = float(mdl_dict.get("dropout", "0.2"))
    recurrent_dropout = float(mdl_dict.get("recurrent_dropout", "0.2"))
    spatial_dropout = float(mdl_dict.get("spatial_dropout", "0.2"))
    units = int(mdl_dict.get("units", "512"))
    batch_size = int(mdl_dict.get("batch_size", "512"))
    activation_recurrent = mdl_dict.get("activation_recurrent", "elu")
    activation_output = mdl_dict.get("activation_output", "sigmoid")

    return {
        "lr": lr,
        "embedding_size": embedding_size,
        "dropout": dropout,
        "recurrent_dropout": recurrent_dropout,
        "spatial_dropout": spatial_dropout,
        "units": units,
        "batch_size": batch_size,
        "activation_recurrent": activation_recurrent,
        "activation_output": activation_output,
    }


def weighted_loss(class_weights):
    """
    Create a weighted loss function. Penalise the misclassification
    of classes more with the higher usage
    """
    weight_values = list(class_weights.values())

    def weighted_binary_crossentropy(y_true, y_pred):
        # add another dimension to compute dot product
        expanded_weights = K.expand_dims(weight_values, axis=-1)
        return K.dot(K.binary_crossentropy(y_true, y_pred), expanded_weights)
    return weighted_binary_crossentropy


def set_recurrent_network(mdl_dict, reverse_dictionary, class_weights):
    """
    Create a RNN network and set its parameters
    """
    dimensions = len(reverse_dictionary) + 1
    model_params = get_best_parameters(mdl_dict)

    # define the architecture of the neural network
    model = Sequential()
    model.add(Embedding(dimensions, model_params["embedding_size"], mask_zero=True))
    model.add(SpatialDropout1D(model_params["spatial_dropout"]))
    model.add(GRU(model_params["units"], dropout=model_params["spatial_dropout"], recurrent_dropout=model_params["recurrent_dropout"], activation=model_params["activation_recurrent"], return_sequences=True))
    model.add(Dropout(model_params["dropout"]))
    model.add(GRU(model_params["units"], dropout=model_params["spatial_dropout"], recurrent_dropout=model_params["recurrent_dropout"], activation=model_params["activation_recurrent"], return_sequences=False))
    model.add(Dropout(model_params["dropout"]))
    model.add(Dense(dimensions, activation=model_params["activation_output"]))
    optimizer = RMSprop(lr=model_params["lr"])
    model.compile(loss=weighted_loss(class_weights), optimizer=optimizer)
    return model, model_params


def compute_precision(model, x, y, reverse_data_dictionary, next_compatible_tools, usage_scores, actual_classes_pos, topk):
    """
    Compute absolute and compatible precision
    """
    absolute_precision = 0.0
    test_sample = np.reshape(x, (1, len(x)))

    # predict next tools for a test path
    prediction = model.predict(test_sample, verbose=0)

    nw_dimension = prediction.shape[1]

    # remove the 0th position as there is no tool at this index
    prediction = np.reshape(prediction, (nw_dimension,))

    prediction_pos = np.argsort(prediction, axis=-1)
    topk_prediction_pos = prediction_pos[-topk:]

    # remove the wrong tool position from the predicted list of tool positions
    topk_prediction_pos = [x for x in topk_prediction_pos if x > 0]

    # read tool names using reverse dictionary
    actual_next_tool_names = [reverse_data_dictionary[int(tool_pos)] for tool_pos in actual_classes_pos]
    top_predicted_next_tool_names = [reverse_data_dictionary[int(tool_pos)] for tool_pos in topk_prediction_pos]

    # compute the class weights of predicted tools
    mean_usg_score = 0
    usg_wt_scores = list()
    for t_id in topk_prediction_pos:
        t_name = reverse_data_dictionary[int(t_id)]
        if t_id in usage_scores and t_name in actual_next_tool_names:
            usg_wt_scores.append(np.log(usage_scores[t_id] + 1.0))
    if len(usg_wt_scores) > 0:
            mean_usg_score = np.sum(usg_wt_scores) / float(topk)
    false_positives = [tool_name for tool_name in top_predicted_next_tool_names if tool_name not in actual_next_tool_names]
    absolute_precision = 1 - (len(false_positives) / float(topk))
    return mean_usg_score, absolute_precision


def verify_model(model, x, y, reverse_data_dictionary, next_compatible_tools, usage_scores, topk_list=[1, 2, 3]):
    """
    Verify the model on test data
    """
    print("Evaluating performance on test data...")
    print("Test data size: %d" % len(y))
    size = y.shape[0]
    precision = np.zeros([len(y), len(topk_list)])
    usage_weights = np.zeros([len(y), len(topk_list)])
    # loop over all the test samples and find prediction precision
    for i in range(size):
        actual_classes_pos = np.where(y[i] > 0)[0]
        for index, abs_topk in enumerate(topk_list):
            abs_mean_usg_score, absolute_precision = compute_precision(model, x[i, :], y, reverse_data_dictionary, next_compatible_tools, usage_scores, actual_classes_pos, abs_topk)
            precision[i][index] = absolute_precision
            usage_weights[i][index] = abs_mean_usg_score
    mean_precision = np.mean(precision, axis=0)
    mean_usage = np.mean(usage_weights, axis=0)
    return mean_precision, mean_usage


def save_model(results, data_dictionary, compatible_next_tools, trained_model_path, class_weights):
    # save files
    trained_model = results["model"]
    best_model_parameters = results["best_parameters"]
    model_config = trained_model.to_json()
    model_weights = trained_model.get_weights()

    model_values = {
        'data_dictionary': data_dictionary,
        'model_config': model_config,
        'best_parameters': best_model_parameters,
        'model_weights': model_weights,
        "compatible_tools": compatible_next_tools,
        "class_weights": class_weights
    }
    set_trained_model(trained_model_path, model_values)
