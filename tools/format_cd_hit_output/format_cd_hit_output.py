#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


def extract_mapping_info(input_mapping_filepath):
    mapping_info = {}
    categories = set([])

    with open(input_mapping_filepath, 'r') as mapping_file:
        for line in mapping_file.readlines():
            split_line = line[:-1].split('\t')
            mapping_info.setdefault(split_line[0], split_line[1])
            categories.add(split_line[1])

    return mapping_info, categories


def init_category_distribution(categories=None):
    cluster_categ_distri = {}
    if categories is not None:
        for category in categories:
            cluster_categ_distri[category] = 0
    return cluster_categ_distri


def flush_cluster_info(cluster_name, cluster_ref_seq, ref_seq_cluster,
                       cluster_categ_distri, categories,
                       output_category_distribution_file, cluster_seq_number):
    if cluster_name != '':
        if categories is not None:
            string = cluster_name
            string += '\t' + str(cluster_seq_number)
            for category in categories:
                string += '\t'
                string += str(cluster_categ_distri[category])
            string += '\n'
            output_category_distribution_file.write(string)

        if cluster_ref_seq == '':
            string = "No reference sequence found for "
            string += cluster_name
            raise ValueError(string)

        ref_seq_cluster.setdefault(cluster_ref_seq, cluster_name)


def extract_cluster_info(args, mapping_info=None, categories=None):
    ref_seq_cluster = {}

    if args.output_category_distribution is not None:
        if mapping_info is None or categories is None:
            string = "A file with category distribution is expected but "
            string += "no mapping information are available"
            raise ValueError(string)
        output_cat_distri_file = open(args.output_category_distribution, 'w')
        output_cat_distri_file.write('Cluster\tSequence_number')
        for category in categories:
            output_cat_distri_file.write('\t' + category)

        output_cat_distri_file.write('\n')
    else:
        output_cat_distri_file = None

    with open(args.input_cluster_info, 'r') as cluster_info_file:
        cluster_name = ''
        cluster_categ_distri = init_category_distribution(categories)
        cluster_ref_seq = ''
        cluster_seq_number = 0
        for line in cluster_info_file.readlines():
            if line[0] == '>':
                flush_cluster_info(
                    cluster_name,
                    cluster_ref_seq,
                    ref_seq_cluster,
                    cluster_categ_distri,
                    categories,
                    output_cat_distri_file,
                    cluster_seq_number)
                cluster_name = line[1:-1]
                cluster_name = cluster_name.replace(' ', '_')
                cluster_categ_distri = init_category_distribution(categories)
                cluster_ref_seq = ''
                cluster_seq_number = 0
            else:
                seq_info = line[:-1].split('\t')[1].split(' ')
                seq_name = seq_info[1][1:-3]
                cluster_seq_number += 1

                if categories is not None:
                    seq_count = 1
                    if args.number_sum is not None:
                        if seq_name.find('size') != -1:
                            substring = seq_name[seq_name.find('size'):-1]
                            seq_count = int(substring.split('=')[1])
                    if seq_name not in mapping_info:
                        string = seq_name + " not found in mapping"
                        raise ValueError(string)
                    category = mapping_info[seq_name]
                    cluster_categ_distri[category] += seq_count

                if seq_info[-1] == '*':
                    if cluster_ref_seq != '':
                        string = "A reference sequence (" + cluster_ref_seq
                        string += ") already found for cluster " + cluster_name
                        string += " (" + seq_name + ")"
                        raise ValueError(string)
                    cluster_ref_seq = seq_name

        flush_cluster_info(
            cluster_name,
            cluster_ref_seq,
            ref_seq_cluster,
            cluster_categ_distri,
            categories,
            output_cat_distri_file,
            cluster_seq_number)

    if args.output_category_distribution is not None:
        output_cat_distri_file.close()

    return ref_seq_cluster


def rename_representative_sequences(args, ref_seq_cluster):
    with open(args.input_representative_sequences, 'r') as input_sequences:
        with open(args.output_representative_sequences, 'w') as output_seq:
            for line in input_sequences.readlines():
                if line[0] == '>':
                    seq_name = line[1:-1]
                    if seq_name not in ref_seq_cluster:
                        string = seq_name + " not found as reference sequence"
                        raise ValueError(string)
                    string = '>' + ref_seq_cluster[seq_name] + '\n'
                    output_seq.write(string)
                else:
                    output_seq.write(line)


def format_cd_hit_outputs(args):
    if args.input_mapping is not None:
        mapping_info, categories = extract_mapping_info(args.input_mapping)
    else:
        mapping_info = None
        categories = None

    ref_seq_cluster = extract_cluster_info(args, mapping_info, categories)

    if args.input_representative_sequences is not None:
        rename_representative_sequences(args, ref_seq_cluster)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_cluster_info', required=True)
    parser.add_argument('--input_representative_sequences')
    parser.add_argument('--output_representative_sequences')
    parser.add_argument('--input_mapping')
    parser.add_argument('--output_category_distribution')
    parser.add_argument('--number_sum')
    args = parser.parse_args()

    format_cd_hit_outputs(args)
