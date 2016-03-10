#!/usr/bin/env python

description = """
Explicit Decomposition with Neighborhood (EDeN) utility program.
Protease modelling driver.
"""

epilog = """
Author: Fabrizio Costa
Copyright: 2015
License: GPL
Maintainer: Fabrizio Costa
Email: costa@informatik.uni-freiburg.de
Status: Production

Cite:  Costa, Fabrizio, and Kurt De Grave, 'Fast neighborhood subgraph pairwise
distance kernel', Proceedings of the 26th International Conference on Machine
Learning. 2010. """

import os
import logging

from eden.graph import Vectorizer
from eden.model_base import ModelInitializerBase, main_script
from eden.converter.fasta import fasta_to_sequence
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier
from eden.modifier.seq import mark_modifier
from eden.converter.fasta import sequence_to_eden


class ModelInitializer(ModelInitializerBase):

    def load_data(self, args):
        seqs = fasta_to_sequence(args.input_file)
        return seqs

    def load_positive_data(self, args):
        return self.load_data(args)

    def load_negative_data(self, args):
        seqs = self.load_data(args)
        return seq_to_seq(seqs,
                          modifier=shuffle_modifier,
                          times=args.negative_ratio,
                          order=args.shuffle_order)

    def pre_processor_init(self, args):
        def pre_processor(seqs, **args):
            seqs = seq_to_seq(seqs, modifier=mark_modifier, position=0.5, mark='%')
            seqs = seq_to_seq(seqs, modifier=mark_modifier, position=0.0, mark='@')
            seqs = seq_to_seq(seqs, modifier=mark_modifier, position=1.0, mark='*')
            graphs = sequence_to_eden(seqs)
            return graphs

        pre_processor_parameters = {}
        return pre_processor, pre_processor_parameters

    def vectorizer_init(self, args):
        vectorizer = Vectorizer()
        vectorizer_parameters = {'complexity': [2, 3, 4, 5, 6]}
        return vectorizer, vectorizer_parameters

    def add_arguments(self, parser):
        parser.add_argument('--version', action='version', version='0.1')
        return parser

    def add_arguments_fit(self, parser):
        parser.add_argument("-i", "--input-file",
                            dest="input_file",
                            help="Path to FASTA file containing input sequences.",
                            required=True)
        parser.add_argument("--negative-ratio",
                            dest="negative_ratio",
                            type=int,
                            help="Relative size ration for the randomly permuted negative instances w.r.t.\
                            the positive instances.",
                            default=2)
        parser.add_argument("--shuffle-order",
                            dest="shuffle_order",
                            type=int,
                            help="Order of the k-mer for the random shuffling procedure.",
                            default=2)
        return parser

    def add_arguments_estimate(self, parser):
        return self.add_arguments_fit(parser)

if __name__ == "__main__":
    model_initializer = ModelInitializer()
    main_script(model_initializer=model_initializer,
                description=description,
                epilog=epilog,
                prog_name=os.path.basename(__file__),
                logger=logging.getLogger())
