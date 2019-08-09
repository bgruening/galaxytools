import argparse
import json
import numpy as np
import pandas as pd
import warnings

from scipy.io import mmread
from sklearn.pipeline import Pipeline

from galaxy_ml.utils import (load_model, read_columns,
                             get_module, try_get_attr)


N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))


def main(inputs, infile_estimator, outfile_predict,
         infile_weights=None, infile1=None,
         fasta_path=None, ref_seq=None,
         vcf_path=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_estimator : strgit
        File path to trained estimator input

    outfile_predict : str
        File path to save the prediction results, tabular

    infile_weights : str
        File path to weights input

    infile1 : str
        File path to dataset containing features

    fasta_path : str
        File path to dataset containing fasta file

    ref_seq : str
        File path to dataset containing the reference genome sequence.

    vcf_path : str
        File path to dataset containing variants info.
    """
    warnings.filterwarnings('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    # load model
    with open(infile_estimator, 'rb') as est_handler:
        estimator = load_model(est_handler)

    main_est = estimator
    if isinstance(estimator, Pipeline):
        main_est = estimator.steps[-1][-1]
    if hasattr(main_est, 'config') and hasattr(main_est, 'load_weights'):
        if not infile_weights or infile_weights == 'None':
            raise ValueError("The selected model skeleton asks for weights, "
                             "but dataset for weights wan not selected!")
        main_est.load_weights(infile_weights)

    # handle data input
    input_type = params['input_options']['selected_input']
    # tabular input
    if input_type == 'tabular':
        header = 'infer' if params['input_options']['header1'] else None
        column_option = (params['input_options']
                               ['column_selector_options_1']
                               ['selected_column_selector_option'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = params['input_options']['column_selector_options_1']['col1']
        else:
            c = None

        df = pd.read_csv(infile1, sep='\t', header=header, parse_dates=True)

        X = read_columns(df, c=c, c_option=column_option).astype(float)

        if params['method'] == 'predict':
            preds = estimator.predict(X)
        else:
            preds = estimator.predict_proba(X)

    # sparse input
    elif input_type == 'sparse':
        X = mmread(open(infile1, 'r'))
        if params['method'] == 'predict':
            preds = estimator.predict(X)
        else:
            preds = estimator.predict_proba(X)

    # fasta input
    elif input_type == 'seq_fasta':
        if not hasattr(estimator, 'data_batch_generator'):
            raise ValueError(
                "To do prediction on sequences in fasta input, "
                "the estimator must be a `KerasGBatchClassifier`"
                "equipped with data_batch_generator!")
        pyfaidx = get_module('pyfaidx')
        sequences = pyfaidx.Fasta(fasta_path)
        n_seqs = len(sequences.keys())
        X = np.arange(n_seqs)[:, np.newaxis]
        seq_length = estimator.data_batch_generator.seq_length
        batch_size = getattr(estimator, 'batch_size', 32)
        steps = (n_seqs + batch_size - 1) // batch_size

        seq_type = params['input_options']['seq_type']
        klass = try_get_attr(
            'galaxy_ml.preprocessors', seq_type)

        pred_data_generator = klass(
            fasta_path, seq_length=seq_length)

        if params['method'] == 'predict':
            preds = estimator.predict(
                X, data_generator=pred_data_generator, steps=steps)
        else:
            preds = estimator.predict_proba(
                X, data_generator=pred_data_generator, steps=steps)

    # vcf input
    elif input_type == 'variant_effect':
        klass = try_get_attr('galaxy_ml.preprocessors',
                             'GenomicVariantBatchGenerator')

        options = params['input_options']
        options.pop('selected_input')
        if options['blacklist_regions'] == 'none':
            options['blacklist_regions'] = None

        pred_data_generator = klass(
            ref_genome_path=ref_seq, vcf_path=vcf_path, **options)

        pred_data_generator.fit()

        preds = estimator.model_.predict_generator(
            pred_data_generator.flow(batch_size=32),
            workers=N_JOBS,
            use_multiprocessing=True)

        if preds.min() < 0. or preds.max() > 1.:
            warnings.warn('Network returning invalid probability values. '
                          'The last layer might not normalize predictions '
                          'into probabilities '
                          '(like softmax or sigmoid would).')

        if params['method'] == 'predict_proba' and preds.shape[1] == 1:
            # first column is probability of class 0 and second is of class 1
            preds = np.hstack([1 - preds, preds])

        elif params['method'] == 'predict':
            if preds.shape[-1] > 1:
                # if the last activation is `softmax`, the sum of all
                # probibilities will 1, the classification is considered as
                # multi-class problem, otherwise, we take it as multi-label.
                act = getattr(estimator.model_.layers[-1], 'activation', None)
                if act and act.__name__ == 'softmax':
                    classes = preds.argmax(axis=-1)
                else:
                    preds = (preds > 0.5).astype('int32')
            else:
                classes = (preds > 0.5).astype('int32')

            preds = estimator.classes_[classes]
    # end input

    # output
    if input_type == 'variant_effect':   # TODO: save in batchs
        rval = pd.DataFrame(preds)
        meta = pd.DataFrame(
            pred_data_generator.variants,
            columns=['chrom', 'pos', 'name', 'ref', 'alt', 'strand'])

        rval = pd.concat([meta, rval], axis=1)

    elif len(preds.shape) == 1:
        rval = pd.DataFrame(preds, columns=['Predicted'])
    else:
        rval = pd.DataFrame(preds)

    rval.to_csv(outfile_predict, sep='\t',
                header=True, index=False)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--infile_estimator", dest="infile_estimator")
    aparser.add_argument("-w", "--infile_weights", dest="infile_weights")
    aparser.add_argument("-X", "--infile1", dest="infile1")
    aparser.add_argument("-O", "--outfile_predict", dest="outfile_predict")
    aparser.add_argument("-f", "--fasta_path", dest="fasta_path")
    aparser.add_argument("-r", "--ref_seq", dest="ref_seq")
    aparser.add_argument("-v", "--vcf_path", dest="vcf_path")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.outfile_predict,
         infile_weights=args.infile_weights, infile1=args.infile1,
         fasta_path=args.fasta_path, ref_seq=args.ref_seq,
         vcf_path=args.vcf_path)
