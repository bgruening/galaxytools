#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


def extract_abundances(fp, nb_charact_to_extract):
    abundances = {}
    more_abund_charact = []
    abund_sum = 0
    with open(fp, 'r') as abundance_f:
        for line in abundance_f.readlines()[1:]:
            split_line = line[:-1].split('\t')
            charact_id = split_line[0]
            abund = float(split_line[1])
            abundances[charact_id] = 100*abund
            abund_sum += abundances[charact_id]

            if len(more_abund_charact) < nb_charact_to_extract:
                more_abund_charact.append(charact_id)
            else:
                best_pos = None
                for i in range(len(more_abund_charact)-1, -1, -1):
                    if abundances[more_abund_charact[i]] < abund:
                        best_pos = i
                    else:
                        break
                if best_pos is not None:
                    tmp_more_abund_charact = more_abund_charact
                    more_abund_charact = tmp_more_abund_charact[:best_pos]
                    more_abund_charact += [charact_id]
                    more_abund_charact += tmp_more_abund_charact[best_pos:-1]
    return abundances, more_abund_charact


def format_characteristic_name(all_name):
    if all_name.find(':') != -1:
        charact_id = all_name.split(':')[0]
        char_name = all_name.split(':')[1][1:]
    else:
        charact_id = all_name
        char_name = ''

    char_name = char_name.replace('/', ' ')
    char_name = char_name.replace('-', ' ')
    char_name = char_name.replace("'", '')
    if char_name.find('(') != -1 and char_name.find(')') != -1:
        open_bracket = char_name.find('(')
        close_bracket = char_name.find(')')+1
        char_name = char_name[:open_bracket] + char_name[close_bracket:]
    return charact_id, char_name


def write_more_abundant_charat(abundances, more_abund_charact, output_fp):
    with open(output_fp, 'w') as output_f:
        output_f.write('id\tname\t%s\n' % '\t'.join(abundances.keys()))

        for mac in more_abund_charact:
            charact_id, charact_name = format_characteristic_name(mac)
            output_f.write('%s\t%s' % (charact_id, charact_name))
            for sample in abundances:
                abund = abundances[sample].get(mac, 0)
                output_f.write('\t%s' % (abund))
            output_f.write('\n')


def extract_similar_characteristics(abund, sim_output_fp, output_files):
    abund_keys = list(abund)
    sim_characteristics = set(abund[abund_keys[0]].keys())
    for sample in abund_keys[1:]:
        sim_characteristics.intersection_update(abund[sample].keys())
    print('Similar between all samples: %s' % len(sim_characteristics))

    with open(sim_output_fp, 'w') as sim_output_f:
        sim_output_f.write('id\tname\t%s\n' % '\t'.join(abund_keys))
        for charact in list(sim_characteristics):
            charact_id, charact_name = format_characteristic_name(charact)
            sim_output_f.write('%s\t%s' % (charact_id, charact_name))
            for sample in abund_keys:
                sim_output_f.write('\t%s' % abund[sample][charact])
            sim_output_f.write('\n')

    print('Specific to samples:')
    diff_char = {}
    for i in range(len(abund_keys)):
        sample = abund_keys[i]
        print(' %s' % sample )
        print('    All: %s' % len(abund[sample].keys()))
        diff_char[sample] = set(abund[sample].keys())
        diff_char[sample].difference_update(sim_characteristics)
        perc = 100*len(diff_char[sample])/(1.*len(abund[sample].keys()))
        print('    Number of specific characteristics: %s' % len(diff_char[sample]))
        print('    Percentage of specific characteristics: %s' % perc)

        relative_abundance = 0
        with open(output_files[i], 'w') as output_f:
            output_f.write('id\tname\tabundances\n')
            for charact in list(diff_char[sample]):
                charact_id, charact_name = format_characteristic_name(charact)
                output_f.write('%s\t%s' % (charact_id, charact_name))
                output_f.write('%s\n' % abund[sample][charact])
                relative_abundance += abund[sample][charact]
        print('    Relative abundance of specific characteristics: %s' % relative_abundance)

    return sim_characteristics


def compare_humann2_output(args):
    abund = {}
    more_abund_charact = []

    for i in range(len(args.sample_name)):
        abund[args.sample_name[i]], mac = extract_abundances(
            args.charact_input_fp[i],
            args.most_abundant_characteristics_to_extract)
        more_abund_charact += mac

    write_more_abundant_charat(
        abund,
        list(set(more_abund_charact)),
        args.more_abundant_output_fp)
    extract_similar_characteristics(
        abund,
        args.similar_output_fp,
        args.specific_output_fp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', required=True, action='append')
    parser.add_argument('--charact_input_fp', required=True, action='append')
    parser.add_argument(
        '--most_abundant_characteristics_to_extract',
        required=True,
        type=int)
    parser.add_argument('--more_abundant_output_fp', required=True)
    parser.add_argument('--similar_output_fp', required=True)
    parser.add_argument(
        '--specific_output_fp',
        required=True,
        action='append')
    args = parser.parse_args()

    if len(args.sample_name) != len(args.charact_input_fp):
        string = "Same number of values (in same order) are expected for "
        string += "--sample_name and --charact_input_fp"
        raise ValueError(string)
    if len(args.sample_name) != len(args.specific_output_fp):
        string = "Same number of values (in same order) are expected for "
        string += "--sample_name and --specific_output_fp"
        raise ValueError(string)

    compare_humann2_output(args)
