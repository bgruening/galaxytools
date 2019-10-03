"""
Predict next tools in the Galaxy workflows
using machine learning (recurrent neural network)
"""

import numpy as np
import argparse
import time

# machine learning library
import keras.callbacks as callbacks

import extract_workflow_connections
import prepare_data
import optimise_hyperparameters
import utils


class PredictTool:

    @classmethod
    def __init__(self):
        """ Init method. """

    @classmethod
    def find_train_best_network(self, network_config, reverse_dictionary, train_data, train_labels, test_data, test_labels, n_epochs, class_weights, usage_pred, compatible_next_tools):
        """
        Define recurrent neural network and train sequential data
        """
        print("Start hyperparameter optimisation...")
        hyper_opt = optimise_hyperparameters.HyperparameterOptimisation()
        best_params = hyper_opt.train_model(network_config, reverse_dictionary, train_data, train_labels, class_weights)

        # retrieve the model and train on complete dataset without validation set
        model, best_params = utils.set_recurrent_network(best_params, reverse_dictionary, class_weights)

        # define callbacks
        early_stopping = callbacks.EarlyStopping(monitor='loss', mode='min', min_delta=1e-1, verbose=1, patience=0)
        predict_callback_test = PredictCallback(test_data, test_labels, reverse_dictionary, n_epochs, compatible_next_tools, usage_pred)

        callbacks_list = [predict_callback_test, early_stopping]

        print("Start training on the best model...")
        train_performance = dict()
        if len(test_data) > 0:
            trained_model = model.fit(
                train_data,
                train_labels,
                batch_size=int(best_params["batch_size"]),
                epochs=n_epochs,
                verbose=2,
                callbacks=callbacks_list,
                shuffle="batch",
                validation_data=(test_data, test_labels)
            )
            train_performance["validation_loss"] = np.array(trained_model.history["val_loss"])
            train_performance["precision"] = predict_callback_test.precision
            train_performance["usage_weights"] = predict_callback_test.usage_weights
        else:
            trained_model = model.fit(
                train_data,
                train_labels,
                batch_size=int(best_params["batch_size"]),
                epochs=n_epochs,
                verbose=2,
                callbacks=callbacks_list,
                shuffle="batch"
            )
        train_performance["train_loss"] = np.array(trained_model.history["loss"])
        train_performance["model"] = model
        train_performance["best_parameters"] = best_params
        return train_performance


class PredictCallback(callbacks.Callback):
    def __init__(self, test_data, test_labels, reverse_data_dictionary, n_epochs, next_compatible_tools, usg_scores):
        self.test_data = test_data
        self.test_labels = test_labels
        self.reverse_data_dictionary = reverse_data_dictionary
        self.precision = list()
        self.usage_weights = list()
        self.n_epochs = n_epochs
        self.next_compatible_tools = next_compatible_tools
        self.pred_usage_scores = usg_scores

    def on_epoch_end(self, epoch, logs={}):
        """
        Compute absolute and compatible precision for test data
        """
        if len(self.test_data) > 0:
            precision, usage_weights = utils.verify_model(self.model, self.test_data, self.test_labels, self.reverse_data_dictionary, self.next_compatible_tools, self.pred_usage_scores)
            self.precision.append(precision)
            self.usage_weights.append(usage_weights)
            print("Epoch %d precision: %s" % (epoch + 1, precision))
            print("Epoch %d usage weights: %s" % (epoch + 1, usage_weights))


if __name__ == "__main__":
    start_time = time.time()
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-wf", "--workflow_file", required=True, help="workflows tabular file")
    arg_parser.add_argument("-tu", "--tool_usage_file", required=True, help="tool usage file")
    arg_parser.add_argument("-om", "--output_model", required=True, help="trained model file")
    # data parameters
    arg_parser.add_argument("-cd", "--cutoff_date", required=True, help="earliest date for taking tool usage")
    arg_parser.add_argument("-pl", "--maximum_path_length", required=True, help="maximum length of tool path")
    arg_parser.add_argument("-ep", "--n_epochs", required=True, help="number of iterations to run to create model")
    arg_parser.add_argument("-oe", "--optimize_n_epochs", required=True, help="number of iterations to run to find best model parameters")
    arg_parser.add_argument("-me", "--max_evals", required=True, help="maximum number of configuration evaluations")
    arg_parser.add_argument("-ts", "--test_share", required=True, help="share of data to be used for testing")
    arg_parser.add_argument("-vs", "--validation_share", required=True, help="share of data to be used for validation")
    # neural network parameters
    arg_parser.add_argument("-bs", "--batch_size", required=True, help="size of the tranining batch i.e. the number of samples per batch")
    arg_parser.add_argument("-ut", "--units", required=True, help="number of hidden recurrent units")
    arg_parser.add_argument("-es", "--embedding_size", required=True, help="size of the fixed vector learned for each tool")
    arg_parser.add_argument("-dt", "--dropout", required=True, help="percentage of neurons to be dropped")
    arg_parser.add_argument("-sd", "--spatial_dropout", required=True, help="1d dropout used for embedding layer")
    arg_parser.add_argument("-rd", "--recurrent_dropout", required=True, help="dropout for the recurrent layers")
    arg_parser.add_argument("-lr", "--learning_rate", required=True, help="learning rate")
    arg_parser.add_argument("-ar", "--activation_recurrent", required=True, help="activation function for recurrent layers")
    arg_parser.add_argument("-ao", "--activation_output", required=True, help="activation function for output layers")
    # get argument values
    args = vars(arg_parser.parse_args())
    tool_usage_path = args["tool_usage_file"]
    workflows_path = args["workflow_file"]
    cutoff_date = args["cutoff_date"]
    maximum_path_length = int(args["maximum_path_length"])
    trained_model_path = args["output_model"]
    n_epochs = int(args["n_epochs"])
    optimize_n_epochs = int(args["optimize_n_epochs"])
    max_evals = int(args["max_evals"])
    test_share = float(args["test_share"])
    validation_share = float(args["validation_share"])
    batch_size = args["batch_size"]
    units = args["units"]
    embedding_size = args["embedding_size"]
    dropout = args["dropout"]
    spatial_dropout = args["spatial_dropout"]
    recurrent_dropout = args["recurrent_dropout"]
    learning_rate = args["learning_rate"]
    activation_recurrent = args["activation_recurrent"]
    activation_output = args["activation_output"]

    config = {
        'cutoff_date': cutoff_date,
        'maximum_path_length': maximum_path_length,
        'n_epochs': n_epochs,
        'optimize_n_epochs': optimize_n_epochs,
        'max_evals': max_evals,
        'test_share': test_share,
        'validation_share': validation_share,
        'batch_size': batch_size,
        'units': units,
        'embedding_size': embedding_size,
        'dropout': dropout,
        'spatial_dropout': spatial_dropout,
        'recurrent_dropout': recurrent_dropout,
        'learning_rate': learning_rate,
        'activation_recurrent': activation_recurrent,
        'activation_output': activation_output
    }

    # Extract and process workflows
    connections = extract_workflow_connections.ExtractWorkflowConnections()
    workflow_paths, compatible_next_tools = connections.read_tabular_file(workflows_path)
    # Process the paths from workflows
    print("Dividing data...")
    data = prepare_data.PrepareData(maximum_path_length, test_share)
    train_data, train_labels, test_data, test_labels, data_dictionary, reverse_dictionary, class_weights, usage_pred = data.get_data_labels_matrices(workflow_paths, tool_usage_path, cutoff_date, compatible_next_tools)
    # find the best model and start training
    predict_tool = PredictTool()
    # start training with weighted classes
    print("Training with weighted classes and samples ...")
    results_weighted = predict_tool.find_train_best_network(config, reverse_dictionary, train_data, train_labels, test_data, test_labels, n_epochs, class_weights, usage_pred, compatible_next_tools)
    print()
    print("Best parameters \n")
    print(results_weighted["best_parameters"])
    print()
    utils.save_model(results_weighted, data_dictionary, compatible_next_tools, trained_model_path, class_weights)
    end_time = time.time()
    print()
    print("Program finished in %s seconds" % str(end_time - start_time))
