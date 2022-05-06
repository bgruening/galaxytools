"""
Predict next tools in the Galaxy workflows
using machine learning (recurrent neural network)
"""

import argparse
import time

import extract_workflow_connections
import keras.callbacks as callbacks
import numpy as np
import optimise_hyperparameters
import prepare_data
import utils


class PredictTool:
    def __init__(self, num_cpus):
        """ Init method. """

    def find_train_best_network(
        self,
        network_config,
        reverse_dictionary,
        train_data,
        train_labels,
        test_data,
        test_labels,
        n_epochs,
        class_weights,
        usage_pred,
        standard_connections,
        tool_freq,
        tool_tr_samples,
    ):
        """
        Define recurrent neural network and train sequential data
        """
        # get tools with lowest representation
        lowest_tool_ids = utils.get_lowest_tools(tool_freq)

        print("Start hyperparameter optimisation...")
        hyper_opt = optimise_hyperparameters.HyperparameterOptimisation()
        best_params, best_model = hyper_opt.train_model(
            network_config,
            reverse_dictionary,
            train_data,
            train_labels,
            test_data,
            test_labels,
            tool_tr_samples,
            class_weights,
        )

        # define callbacks
        early_stopping = callbacks.EarlyStopping(
            monitor="loss",
            mode="min",
            verbose=1,
            min_delta=1e-1,
            restore_best_weights=True,
        )
        predict_callback_test = PredictCallback(
            test_data,
            test_labels,
            reverse_dictionary,
            n_epochs,
            usage_pred,
            standard_connections,
            lowest_tool_ids,
        )

        callbacks_list = [predict_callback_test, early_stopping]
        batch_size = int(best_params["batch_size"])

        print("Start training on the best model...")
        train_performance = dict()
        trained_model = best_model.fit_generator(
            utils.balanced_sample_generator(
                train_data,
                train_labels,
                batch_size,
                tool_tr_samples,
                reverse_dictionary,
            ),
            steps_per_epoch=len(train_data) // batch_size,
            epochs=n_epochs,
            callbacks=callbacks_list,
            validation_data=(test_data, test_labels),
            verbose=2,
            shuffle=True,
        )
        train_performance["validation_loss"] = np.array(
            trained_model.history["val_loss"]
        )
        train_performance["precision"] = predict_callback_test.precision
        train_performance["usage_weights"] = predict_callback_test.usage_weights
        train_performance[
            "published_precision"
        ] = predict_callback_test.published_precision
        train_performance[
            "lowest_pub_precision"
        ] = predict_callback_test.lowest_pub_precision
        train_performance[
            "lowest_norm_precision"
        ] = predict_callback_test.lowest_norm_precision
        train_performance["train_loss"] = np.array(trained_model.history["loss"])
        train_performance["model"] = best_model
        train_performance["best_parameters"] = best_params
        return train_performance


class PredictCallback(callbacks.Callback):
    def __init__(
        self,
        test_data,
        test_labels,
        reverse_data_dictionary,
        n_epochs,
        usg_scores,
        standard_connections,
        lowest_tool_ids,
    ):
        self.test_data = test_data
        self.test_labels = test_labels
        self.reverse_data_dictionary = reverse_data_dictionary
        self.precision = list()
        self.usage_weights = list()
        self.published_precision = list()
        self.n_epochs = n_epochs
        self.pred_usage_scores = usg_scores
        self.standard_connections = standard_connections
        self.lowest_tool_ids = lowest_tool_ids
        self.lowest_pub_precision = list()
        self.lowest_norm_precision = list()

    def on_epoch_end(self, epoch, logs={}):
        """
        Compute absolute and compatible precision for test data
        """
        if len(self.test_data) > 0:
            (
                usage_weights,
                precision,
                precision_pub,
                low_pub_prec,
                low_norm_prec,
                low_num,
            ) = utils.verify_model(
                self.model,
                self.test_data,
                self.test_labels,
                self.reverse_data_dictionary,
                self.pred_usage_scores,
                self.standard_connections,
                self.lowest_tool_ids,
            )
            self.precision.append(precision)
            self.usage_weights.append(usage_weights)
            self.published_precision.append(precision_pub)
            self.lowest_pub_precision.append(low_pub_prec)
            self.lowest_norm_precision.append(low_norm_prec)
            print("Epoch %d usage weights: %s" % (epoch + 1, usage_weights))
            print("Epoch %d normal precision: %s" % (epoch + 1, precision))
            print("Epoch %d published precision: %s" % (epoch + 1, precision_pub))
            print("Epoch %d lowest published precision: %s" % (epoch + 1, low_pub_prec))
            print("Epoch %d lowest normal precision: %s" % (epoch + 1, low_norm_prec))
            print(
                "Epoch %d number of test samples with lowest tool ids: %s"
                % (epoch + 1, low_num)
            )


if __name__ == "__main__":
    start_time = time.time()

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-wf", "--workflow_file", required=True, help="workflows tabular file"
    )
    arg_parser.add_argument(
        "-tu", "--tool_usage_file", required=True, help="tool usage file"
    )
    arg_parser.add_argument(
        "-om", "--output_model", required=True, help="trained model file"
    )
    # data parameters
    arg_parser.add_argument(
        "-cd",
        "--cutoff_date",
        required=True,
        help="earliest date for taking tool usage",
    )
    arg_parser.add_argument(
        "-pl",
        "--maximum_path_length",
        required=True,
        help="maximum length of tool path",
    )
    arg_parser.add_argument(
        "-ep",
        "--n_epochs",
        required=True,
        help="number of iterations to run to create model",
    )
    arg_parser.add_argument(
        "-oe",
        "--optimize_n_epochs",
        required=True,
        help="number of iterations to run to find best model parameters",
    )
    arg_parser.add_argument(
        "-me",
        "--max_evals",
        required=True,
        help="maximum number of configuration evaluations",
    )
    arg_parser.add_argument(
        "-ts",
        "--test_share",
        required=True,
        help="share of data to be used for testing",
    )
    # neural network parameters
    arg_parser.add_argument(
        "-bs",
        "--batch_size",
        required=True,
        help="size of the tranining batch i.e. the number of samples per batch",
    )
    arg_parser.add_argument(
        "-ut", "--units", required=True, help="number of hidden recurrent units"
    )
    arg_parser.add_argument(
        "-es",
        "--embedding_size",
        required=True,
        help="size of the fixed vector learned for each tool",
    )
    arg_parser.add_argument(
        "-dt", "--dropout", required=True, help="percentage of neurons to be dropped"
    )
    arg_parser.add_argument(
        "-sd",
        "--spatial_dropout",
        required=True,
        help="1d dropout used for embedding layer",
    )
    arg_parser.add_argument(
        "-rd",
        "--recurrent_dropout",
        required=True,
        help="dropout for the recurrent layers",
    )
    arg_parser.add_argument(
        "-lr", "--learning_rate", required=True, help="learning rate"
    )

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
    batch_size = args["batch_size"]
    units = args["units"]
    embedding_size = args["embedding_size"]
    dropout = args["dropout"]
    spatial_dropout = args["spatial_dropout"]
    recurrent_dropout = args["recurrent_dropout"]
    learning_rate = args["learning_rate"]
    num_cpus = 16

    config = {
        "cutoff_date": cutoff_date,
        "maximum_path_length": maximum_path_length,
        "n_epochs": n_epochs,
        "optimize_n_epochs": optimize_n_epochs,
        "max_evals": max_evals,
        "test_share": test_share,
        "batch_size": batch_size,
        "units": units,
        "embedding_size": embedding_size,
        "dropout": dropout,
        "spatial_dropout": spatial_dropout,
        "recurrent_dropout": recurrent_dropout,
        "learning_rate": learning_rate,
    }

    # Extract and process workflows
    connections = extract_workflow_connections.ExtractWorkflowConnections()
    (
        workflow_paths,
        compatible_next_tools,
        standard_connections,
    ) = connections.read_tabular_file(workflows_path)
    # Process the paths from workflows
    print("Dividing data...")
    data = prepare_data.PrepareData(maximum_path_length, test_share)
    (
        train_data,
        train_labels,
        test_data,
        test_labels,
        data_dictionary,
        reverse_dictionary,
        class_weights,
        usage_pred,
        train_tool_freq,
        tool_tr_samples,
    ) = data.get_data_labels_matrices(
        workflow_paths,
        tool_usage_path,
        cutoff_date,
        compatible_next_tools,
        standard_connections,
    )
    # find the best model and start training
    predict_tool = PredictTool(num_cpus)
    # start training with weighted classes
    print("Training with weighted classes and samples ...")
    results_weighted = predict_tool.find_train_best_network(
        config,
        reverse_dictionary,
        train_data,
        train_labels,
        test_data,
        test_labels,
        n_epochs,
        class_weights,
        usage_pred,
        standard_connections,
        train_tool_freq,
        tool_tr_samples,
    )
    utils.save_model(
        results_weighted,
        data_dictionary,
        compatible_next_tools,
        trained_model_path,
        class_weights,
        standard_connections,
    )
    end_time = time.time()
    print("Program finished in %s seconds" % str(end_time - start_time))
