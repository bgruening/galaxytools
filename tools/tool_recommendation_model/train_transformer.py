import os
import numpy as np

import tensorflow as tf
from tensorflow.keras.layers import Dense, Dropout, GlobalAveragePooling1D, Input
from tensorflow.keras.models import Model

import utils
import transformer_network

base_path = "log/"
model_path = base_path + "saved_model/"


def create_model(vocab_size, config):
    embed_dim = config["embedding_dim"]
    ff_dim = config["feed_forward_dim"]
    max_len = config["maximum_path_length"]
    dropout = config["dropout"]

    inputs = Input(shape=(max_len,))
    embedding_layer = transformer_network.TokenAndPositionEmbedding(max_len, vocab_size, embed_dim)
    x = embedding_layer(inputs)
    transformer_block = transformer_network.TransformerBlock(embed_dim, config["n_heads"], ff_dim)
    x, weights = transformer_block(x)
    x = GlobalAveragePooling1D()(x)
    x = Dropout(dropout)(x)
    x = Dense(ff_dim, activation="relu")(x)
    x = Dropout(dropout)(x)
    outputs = Dense(vocab_size, activation="sigmoid")(x)
    return Model(inputs=inputs, outputs=[outputs, weights])


def create_enc_transformer(train_data, train_labels, test_data, test_labels, f_dict, r_dict, c_wts, c_tools, pub_conn, tr_t_freq, config):
    print("Train transformer...")
    vocab_size = len(f_dict) + 1
    maxlen = config["maximum_path_length"]

    enc_optimizer = tf.keras.optimizers.Adam(learning_rate=config["learning_rate"])

    model = create_model(vocab_size, config)

    u_tr_y_labels, u_tr_y_labels_dict = utils.get_u_tr_labels(train_labels)
    u_te_y_labels, u_te_y_labels_dict = utils.get_u_tr_labels(test_labels)

    trained_on_labels = [int(item) for item in list(u_tr_y_labels_dict.keys())]

    epo_tr_batch_loss = list()
    epo_tr_batch_acc = list()
    epo_tr_batch_categorical_loss = list()
    epo_te_batch_loss = list()
    epo_te_batch_acc = list()
    epo_te_batch_categorical_loss = list()
    all_sel_tool_ids = list()
    epo_te_precision = list()
    epo_low_te_precision = list()

    te_lowest_t_ids = utils.get_low_freq_te_samples(test_data, test_labels, tr_t_freq)
    tr_log_step = config["tr_logging_step"]
    te_log_step = config["te_logging_step"]
    n_train_steps = config["n_train_iter"]
    te_batch_size = config["te_batch_size"]
    tr_batch_size = config["tr_batch_size"]
    is_transformer = config["is_transformer"]

    sel_tools = list()
    for batch in range(n_train_steps):
        print("Total train data size: ", train_data.shape, train_labels.shape)
        x_train, y_train, sel_tools = utils.sample_balanced_tr_y(train_data, train_labels, u_tr_y_labels_dict, tr_batch_size, tr_t_freq, sel_tools)
        print("Batch train data size: ", x_train.shape, y_train.shape)
        all_sel_tool_ids.extend(sel_tools)
        with tf.GradientTape() as model_tape:
            prediction, att_weights = model(x_train, training=True)
            tr_loss, tr_cat_loss = utils.compute_loss(y_train, prediction)
            tr_acc = tf.reduce_mean(utils.compute_acc(y_train, prediction))
        trainable_vars = model.trainable_variables
        model_gradients = model_tape.gradient(tr_loss, trainable_vars)
        enc_optimizer.apply_gradients(zip(model_gradients, trainable_vars))
        epo_tr_batch_loss.append(tr_loss.numpy())
        epo_tr_batch_acc.append(tr_acc.numpy())
        epo_tr_batch_categorical_loss.append(tr_cat_loss.numpy())
        print("Step {}/{}, training binary loss: {}, categorical_loss: {}, training accuracy: {}".format(batch+1, n_train_steps, tr_loss.numpy(), tr_cat_loss.numpy(), tr_acc.numpy()))
        if (batch+1) % te_log_step == 0:
            print("Predicting on test data...")
            te_loss, te_acc, test_cat_loss, te_prec, low_te_prec = utils.validate_model(test_data, test_labels, te_batch_size, model, f_dict, r_dict, u_te_y_labels_dict, trained_on_labels, te_lowest_t_ids, is_transformer)
            epo_te_batch_loss.append(te_loss)
            epo_te_batch_acc.append(te_acc)
            epo_te_batch_categorical_loss.append(test_cat_loss)
            epo_te_precision.append(te_prec)
            epo_low_te_precision.append(low_te_prec)
        print()
        if (batch+1) % tr_log_step == 0:
            print("Saving model at training step {}/{}".format(batch + 1, n_train_steps))
            tf_path = model_path + "{}/".format(batch+1)
            tf_model_save = model_path + "{}/tf_model/".format(batch+1)
            tf_model_save_h5 = model_path + "{}/tf_model_h5/".format(batch+1)
            if not os.path.isdir(tf_path):
                os.mkdir(tf_path)
                os.mkdir(tf_model_save)
                os.mkdir(tf_model_save_h5)
            tf.saved_model.save(model, tf_model_save)
            utils.save_model_file(tf_model_save_h5, model, r_dict, c_wts, c_tools, pub_conn)
    new_dict = dict()
    for k in u_tr_y_labels_dict:
        new_dict[str(k)] = ",".join([str(item) for item in u_tr_y_labels_dict[k]])

    utils.write_file(base_path + "data/epo_tr_batch_loss.txt", ",".join([str(item) for item in epo_tr_batch_loss]))
    utils.write_file(base_path + "data/epo_tr_batch_acc.txt", ",".join([str(item) for item in epo_tr_batch_acc]))
    utils.write_file(base_path + "data/epo_te_batch_loss.txt", ",".join([str(item) for item in epo_te_batch_loss]))
    utils.write_file(base_path + "data/epo_te_batch_acc.txt", ",".join([str(item) for item in epo_te_batch_acc]))
    utils.write_file(base_path + "data/epo_tr_batch_categorical_loss.txt", ",".join([str(item) for item in epo_tr_batch_categorical_loss]))
    utils.write_file(base_path + "data/epo_te_batch_categorical_loss.txt", ",".join([str(item) for item in epo_te_batch_categorical_loss]))
    utils.write_file(base_path + "data/epo_te_precision.txt", ",".join([str(item) for item in epo_te_precision]))
    utils.write_file(base_path + "data/all_sel_tool_ids.txt", ",".join([str(item) for item in all_sel_tool_ids]))
    utils.write_file(base_path + "data/epo_low_te_precision.txt", ",".join([str(item) for item in epo_low_te_precision]))
    utils.write_file(base_path + "data/u_tr_y_labels_dict.txt", new_dict)
    utils.write_file(base_path + "data/te_lowest_t_ids.txt", ",".join([str(item) for item in te_lowest_t_ids]))
