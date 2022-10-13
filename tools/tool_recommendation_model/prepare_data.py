"""
Prepare the workflow paths to be used by downstream
machine learning algorithm. The paths are divided
into the test and training sets
"""

import collections
import numpy as np
import random

from sklearn.model_selection import train_test_split

import predict_tool_usage


class PrepareData:

    def __init__(self, max_seq_length, test_data_share):
        """ Init method. """
        self.max_tool_sequence_len = max_seq_length
        self.test_share = test_data_share

    def process_workflow_paths(self, workflow_paths):
        """
        Get all the tools and complete set of individual paths for each workflow
        """
        tokens = list()
        raw_paths = workflow_paths
        raw_paths = [x.replace("\n", '') for x in raw_paths]
        for item in raw_paths:
            split_items = item.split(",")
            for token in split_items:
                if token != "":
                    tokens.append(token)
        tokens = list(set(tokens))
        tokens = np.array(tokens)
        tokens = np.reshape(tokens, [-1, ])
        return tokens, raw_paths

    def create_new_dict(self, new_data_dict):
        """
        Create new data dictionary
        """
        reverse_dict = dict((v, k) for k, v in new_data_dict.items())
        return new_data_dict, reverse_dict

    def assemble_dictionary(self, new_data_dict, old_data_dictionary={}):
        """
        Create/update tools indices in the forward and backward dictionary
        """
        new_data_dict, reverse_dict = self.create_new_dict(new_data_dict)
        return new_data_dict, reverse_dict

    def create_data_dictionary(self, words, old_data_dictionary={}):
        """
        Create two dictionaries having tools names and their indexes
        """
        count = collections.Counter(words).most_common()
        dictionary = dict()
        for index, (word, _) in enumerate(count):
            word = word.lstrip()
            word = word.rstrip()
            dictionary[word] = len(dictionary) + 1
        dictionary, reverse_dictionary = self.assemble_dictionary(dictionary, old_data_dictionary)
        return dictionary, reverse_dictionary

    def decompose_paths(self, paths, dictionary):
        """
        Decompose the paths to variable length sub-paths keeping the first tool fixed
        """
        max_len = 0
        sub_paths_pos = list()
        for index, item in enumerate(paths):
            tools = item.split(",")
            len_tools = len(tools)
            if len_tools > max_len:
                max_len = len_tools
            if len_tools < self.max_tool_sequence_len:
                sequence = tools[0: len_tools]
                tools_pos = [str(dictionary[str(tool_item)]) for tool_item in sequence]
                if len(tools_pos) > 1:
                    sub_paths_pos.append(",".join(tools_pos))
        sub_paths_pos = list(set(sub_paths_pos))
        print("Max length of tools: ", max_len)
        return sub_paths_pos

    def prepare_input_one_target_paths(self, dictionary, reverse_dictionary, paths):
        input_target_paths = dict()
        compatible_tools = dict()
        d_size = 0
        for i, item in enumerate(paths):
            input_tools = item.split(",")
            tool_seq = input_tools
            i_tools = ",".join(tool_seq[0:-1])
            last_i_tool = i_tools.split(",")[-1]
            if last_i_tool not in compatible_tools:
                compatible_tools[last_i_tool] = list()
            t_tools = tool_seq[-1]
            if t_tools not in compatible_tools[last_i_tool]:
                compatible_tools[last_i_tool].append(t_tools)
            if i_tools not in input_target_paths:
                input_target_paths[i_tools] = list()
            if t_tools not in input_target_paths[i_tools]:
                input_target_paths[i_tools].append(t_tools)
            if i_tools not in input_target_paths:
                input_target_paths[i_tools] = list()
            if t_tools not in input_target_paths[i_tools]:
                input_target_paths[i_tools].append(t_tools)
        for item in input_target_paths:
            d_size += len(input_target_paths[item])
        print("Dataset size:", d_size)
        return input_target_paths, compatible_tools, d_size

    def prepare_input_target_paths(self, dictionary, reverse_dictionary, paths):
        input_target_paths = dict()
        compatible_tools = dict()
        d_size = 0
        for i, item in enumerate(paths):
            input_tools = item.split(",")
            ctr = 0
            for ctr in range(len(input_tools) - 1):
                # uncomment this for one token target idea
                tool_seq = input_tools[0: ctr + 2]
                i_tools = ",".join(tool_seq[0:-1])
                last_i_tool = i_tools.split(",")[-1]
                if last_i_tool not in compatible_tools:
                    compatible_tools[last_i_tool] = list()
                t_tools = tool_seq[-1]
                if t_tools not in compatible_tools[last_i_tool]:
                    compatible_tools[last_i_tool].append(t_tools)
                if i_tools not in input_target_paths:
                    input_target_paths[i_tools] = list()
                if t_tools not in input_target_paths[i_tools]:
                    input_target_paths[i_tools].append(t_tools)
                if i_tools not in input_target_paths:
                    input_target_paths[i_tools] = list()
                if t_tools not in input_target_paths[i_tools]:
                    input_target_paths[i_tools].append(t_tools)
        for item in input_target_paths:
            d_size += len(input_target_paths[item])
        print("Dataset size:", d_size)
        return input_target_paths, compatible_tools, d_size

    def pad_paths_one_tool_target(self, multi_paths, compatible_tools, d_size, rev_dict, dictionary):
        d_size = len(multi_paths)
        input_mat = np.zeros([d_size, self.max_tool_sequence_len])
        target_mat = np.zeros([d_size, len(dictionary) + 1])
        train_counter = 0
        for input_seq, target_seq_tools in list(multi_paths.items()):
            input_seq_tools = input_seq.split(",")
            last_i_tool = input_seq_tools[-1]
            for id_pos, pos in enumerate(input_seq_tools):
                input_mat[train_counter][id_pos] = int(pos)
            if last_i_tool in compatible_tools:
                compatible_targets = compatible_tools[last_i_tool]
            for k, t_label in enumerate(target_seq_tools):
                target_mat[train_counter][int(t_label)] = 1
            for c_tool in compatible_targets:
                target_mat[train_counter][int(c_tool)] = 1
            train_counter += 1
        print("Final data size: ", input_mat.shape, target_mat.shape)
        train_data, test_data, train_labels, test_labels = train_test_split(input_mat, target_mat, test_size=self.test_share, random_state=42)
        return train_data, train_labels, test_data, test_labels

    def split_test_train_data(self, multilabels_paths):
        """
        Split into test and train data randomly for each run
        """
        train_dict = dict()
        test_dict = dict()
        all_paths = multilabels_paths.keys()
        random.shuffle(list(all_paths))
        split_number = int(self.test_share * len(all_paths))
        for index, path in enumerate(list(all_paths)):
            if index < split_number:
                test_dict[path] = multilabels_paths[path]
            else:
                train_dict[path] = multilabels_paths[path]
        return train_dict, test_dict

    def get_predicted_usage(self, data_dictionary, predicted_usage):
        """
        Get predicted usage for tools
        """
        usage = dict()
        epsilon = 0.0
        # index 0 does not belong to any tool
        usage[0] = epsilon
        for k, v in data_dictionary.items():
            try:
                usg = predicted_usage[k]
                if usg < epsilon:
                    usg = epsilon
                usage[v] = usg
            except Exception:
                usage[v] = epsilon
                continue
        return usage

    def assign_class_weights(self, n_classes, predicted_usage):
        """
        Compute class weights using usage
        """
        class_weights = dict()
        class_weights[str(0)] = 0.0
        for key in range(1, n_classes + 1):
            u_score = predicted_usage[key]
            if u_score < 1.0:
                u_score += 1.0
            class_weights[key] = np.round(np.log(u_score), 6)
        return class_weights

    def get_train_tool_labels_freq(self, train_paths, reverse_dictionary):
        """
        Get the frequency of last tool of each tool sequence
        to estimate the frequency of tool sequences
        """
        last_tool_freq = dict()
        freq_dict_names = dict()
        for path in train_paths:
            tools_pos = np.where(path > 0)[0]
            path_pos = tools_pos
            path_pos = [str(int(item)) for item in path_pos]

            for tool_pos in path_pos:
                if tool_pos not in last_tool_freq:
                    last_tool_freq[tool_pos] = 0
                    freq_dict_names[reverse_dictionary[int(tool_pos)]] = 0
                last_tool_freq[tool_pos] += 1
                freq_dict_names[reverse_dictionary[int(tool_pos)]] += 1
        sorted_dict = dict(sorted(last_tool_freq.items(), key=lambda kv: kv[1], reverse=True))
        return sorted_dict

    def get_train_last_tool_freq(self, train_paths, reverse_dictionary):
        """
        Get the frequency of last tool of each tool sequence
        to estimate the frequency of tool sequences
        """
        last_tool_freq = dict()
        freq_dict_names = dict()
        for path in train_paths:
            tools_pos = np.where(path > 0)[0]
            path_pos = path[tools_pos]
            path_pos = [str(int(item)) for item in path_pos]
            last_tool = path_pos[-1]
            if last_tool not in last_tool_freq:
                last_tool_freq[last_tool] = 0
                freq_dict_names[reverse_dictionary[int(last_tool)]] = 0
            last_tool_freq[last_tool] += 1
            freq_dict_names[reverse_dictionary[int(last_tool)]] += 1
        sorted_dict = dict(sorted(last_tool_freq.items(), key=lambda kv: kv[1], reverse=True))
        return sorted_dict

    def get_toolid_samples(self, train_data, l_tool_freq):
        l_tool_tr_samples = dict()
        for tool_id in l_tool_freq:
            for index, tr_sample in enumerate(train_data):
                last_tool_id = str(int(tr_sample[-1]))
                if last_tool_id == tool_id:
                    if last_tool_id not in l_tool_tr_samples:
                        l_tool_tr_samples[last_tool_id] = list()
                    l_tool_tr_samples[last_tool_id].append(index)
        return l_tool_tr_samples

    def get_data_labels_matrices(self, workflow_paths, usage_df, cutoff_date, standard_connections, old_data_dictionary={}):
        """
        Convert the training and test paths into corresponding numpy matrices
        """
        processed_data, raw_paths = self.process_workflow_paths(workflow_paths)
        dictionary, rev_dict = self.create_data_dictionary(processed_data, old_data_dictionary)

        num_classes = len(dictionary)

        print("Raw paths: %d" % len(raw_paths))
        random.shuffle(raw_paths)

        print("Decomposing paths...")
        all_unique_paths = self.decompose_paths(raw_paths, dictionary)
        random.shuffle(all_unique_paths)

        print("Creating dictionaries...")
        multilabels_paths, compatible_tools, d_size = self.prepare_input_target_paths(dictionary, rev_dict, all_unique_paths)

        print("Complete data: %d" % d_size)

        print("Padding train and test data...")
        # pad training and test data with trailing zeros
        train_data, train_labels, test_data, test_labels = self.pad_paths_one_tool_target(multilabels_paths, compatible_tools, d_size, rev_dict, dictionary)

        print("Train data: ", train_data.shape)
        print("Test data: ", test_data.shape)

        print("Estimating sample frequency...")
        tr_tool_freq = self.get_train_tool_labels_freq(train_labels, rev_dict)

        # Predict tools usage
        print("Predicting tools' usage...")
        usage_pred = predict_tool_usage.ToolPopularity()
        usage = usage_pred.extract_tool_usage(usage_df, cutoff_date, dictionary)
        tool_usage_prediction = usage_pred.get_pupularity_prediction(usage)
        t_pred_usage = self.get_predicted_usage(dictionary, tool_usage_prediction)
        # get class weights using the predicted usage for each tool
        class_weights = self.assign_class_weights(num_classes, t_pred_usage)
        return train_data, train_labels, test_data, test_labels, dictionary, rev_dict, class_weights, compatible_tools, tr_tool_freq
