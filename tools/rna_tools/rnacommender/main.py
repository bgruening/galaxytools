#!/usr/bin/env python
"""Recommendation."""

import argparse
import sys
from rbpfeatures import RBPVectorizer
from data import PredictDataset
from recommend import Predictor

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', metavar='fasta', type=str,
                        help="""Fasta file containing the RBP \
                        sequences.""")

    args = parser.parse_args()

    v = RBPVectorizer(fasta=args.fasta)
    rbp_fea = v.vectorize()

    if rbp_fea is not None:
        # Define and load dataset
        D = PredictDataset(
            fp=rbp_fea, fr="AURA_Human_data/RNA_features/HT_utrs.h5")
        dataset = D.load()

        model = "AURA_Human_data/model/trained_model.pkl"

        # Define the Trainer and train the model
        P = Predictor(predict_dataset=dataset,
                      trained_model=model,
                      serendipity_dic=model + '_',
                      output="output.txt")
        P.predict()
    else:
        sys.stdout.write("""
        ########################################
        WARNING: The queried protein has no domain similarity with the proteins in the training dataset. It cannot be predicted.
        ########################################
        """)
