"""
Predict next tools in the Galaxy workflows
using deep learning learning (Transformers)
"""
import argparse
import time

import extract_workflow_connections
import prepare_data
import train_transformer

if __name__ == "__main__":
    start_time = time.time()

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-wf", "--workflow_file", required=True, help="workflows tabular file")
    arg_parser.add_argument("-tu", "--tool_usage_file", required=True, help="tool usage file")
    # data parameters
    arg_parser.add_argument("-cd", "--cutoff_date", required=True, help="earliest date for taking tool usage")
    arg_parser.add_argument("-pl", "--maximum_path_length", required=True, help="maximum length of tool path")
    arg_parser.add_argument("-om", "--output_model", required=True, help="trained model path")
    # neural network parameters
    arg_parser.add_argument("-ti", "--n_train_iter", required=True, help="Number of training iterations run to create model")
    arg_parser.add_argument("-nhd", "--n_heads", required=True, help="Number of head in transformer's multi-head attention")
    arg_parser.add_argument("-ed", "--n_embed_dim", required=True, help="Embedding dimension")
    arg_parser.add_argument("-fd", "--n_feed_forward_dim", required=True, help="Feed forward network dimension")
    arg_parser.add_argument("-dt", "--dropout", required=True, help="Percentage of neurons to be dropped")
    arg_parser.add_argument("-lr", "--learning_rate", required=True, help="Learning rate")
    arg_parser.add_argument("-ts", "--te_share", required=True, help="Share of data to be used for testing")
    arg_parser.add_argument("-trbs", "--tr_batch_size", required=True, help="Train batch size")
    arg_parser.add_argument("-trlg", "--tr_logging_step", required=True, help="Train logging frequency")
    arg_parser.add_argument("-telg", "--te_logging_step", required=True, help="Test logging frequency")
    arg_parser.add_argument("-tebs", "--te_batch_size", required=True, help="Test batch size")

    # get argument values
    args = vars(arg_parser.parse_args())
    tool_usage_path = args["tool_usage_file"]
    workflows_path = args["workflow_file"]
    cutoff_date = args["cutoff_date"]
    maximum_path_length = int(args["maximum_path_length"])

    n_train_iter = int(args["n_train_iter"])
    te_share = float(args["te_share"])
    tr_batch_size = int(args["tr_batch_size"])
    te_batch_size = int(args["te_batch_size"])

    n_heads = int(args["n_heads"])
    feed_forward_dim = int(args["n_feed_forward_dim"])
    embedding_dim = int(args["n_embed_dim"])
    dropout = float(args["dropout"])
    learning_rate = float(args["learning_rate"])
    te_logging_step = int(args["te_logging_step"])
    tr_logging_step = int(args["tr_logging_step"])
    trained_model_path = args["output_model"]

    config = {
        'cutoff_date': cutoff_date,
        'maximum_path_length': maximum_path_length,
        'n_train_iter': n_train_iter,
        'n_heads': n_heads,
        'feed_forward_dim': feed_forward_dim,
        'embedding_dim': embedding_dim,
        'dropout': dropout,
        'learning_rate': learning_rate,
        'te_share': te_share,
        'te_logging_step': te_logging_step,
        'tr_logging_step': tr_logging_step,
        'tr_batch_size': tr_batch_size,
        'te_batch_size': te_batch_size,
        'trained_model_path': trained_model_path
    }
    print("Preprocessing workflows...")
    # Extract and process workflows
    connections = extract_workflow_connections.ExtractWorkflowConnections()
    # Process raw workflow file
    wf_dataframe, usage_df = connections.process_raw_files(workflows_path, tool_usage_path, config)
    workflow_paths, pub_conn = connections.read_tabular_file(wf_dataframe, config)
    # Process the paths from workflows
    print("Dividing data...")
    data = prepare_data.PrepareData(maximum_path_length, te_share)
    train_data, train_labels, test_data, test_labels, f_dict, r_dict, c_wts, c_tools, tr_tool_freq = data.get_data_labels_matrices(workflow_paths, usage_df, cutoff_date, pub_conn)
    print(train_data.shape, train_labels.shape, test_data.shape, test_labels.shape)
    train_transformer.create_enc_transformer(train_data, train_labels, test_data, test_labels, f_dict, r_dict, c_wts, c_tools, pub_conn, tr_tool_freq, config)
    end_time = time.time()
    print("Program finished in %s seconds" % str(end_time - start_time))
