import argparse
import json
import pandas as pd
import warnings

from galaxy_ml.model_validations import train_test_split
from galaxy_ml.model_validations import OrderedKFold


def main(inputs, infile_array, outfile_train, outfile_test,
         infile_labels=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_array : str
        File paths of input arrays separated by comma

    infile_labels : str
        File path to dataset containing labels

    outfile_train : str
        File path to dataset containing train split

    outfile_test : str
        File path to dataset containing test split
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    input_header = params['header0']
    header = 'infer' if input_header else None
    array = pd.read_csv(infile_array, sep='\t', header=header,
                        parse_dates=True)

    options = params['options']
    shuffle_selection = options.pop('shuffle_selection')
    options['shuffle'] = shuffle_selection['shuffle']
    if infile_labels:
        header = 'infer' if shuffle_selection['header1'] else None
        col_index = shuffle_selection['col'][0] - 1
        df = pd.read_csv(infile_labels, sep='\t', header=header,
                         parse_dates=True)
        labels = df.iloc[:, col_index].values
        options['labels'] = labels

    if shuffle_selection['shuffle'] == 'ordered_target':
        test_size = options['test_size']
        if test_size < 1.0:
            if test_size > 0.5:
                raise ValueError("Ordered Target Split only supports "
                                 "test proportion 0 - 0.5!")
            n_splits = round(1 / test_size)
        else:
            n_samples = array.shape[0]
            n_splits = round(n_samples / test_size)

        splitter = OrderedKFold(n_splits=n_splits, shuffle=True,
                                random_state=options['random_state'])
        train_index, test_index = next(splitter.split(array.values, labels))
        train, test = array.iloc[train_index, :], array.iloc[test_index, :]
    else:
        train, test = train_test_split(array, **options)

    print("Input shape: %s" % repr(array.shape))
    print("Train shape: %s" % repr(train.shape))
    print("Test shape: %s" % repr(test.shape))
    train.to_csv(outfile_train, sep='\t', header=input_header, index=False)
    test.to_csv(outfile_test, sep='\t', header=input_header, index=False)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-X", "--infile_array", dest="infile_array")
    aparser.add_argument("-g", "--infile_labels", dest="infile_labels")
    aparser.add_argument("-o", "--outfile_train", dest="outfile_train")
    aparser.add_argument("-t", "--outfile_test", dest="outfile_test")
    args = aparser.parse_args()

    main(args.inputs, args.infile_array, args.outfile_train,
         args.outfile_test, args.infile_labels)
