import argparse
import json
import warnings

import pandas as pd
from mlxtend.frequent_patterns import association_rules, fpgrowth
from mlxtend.preprocessing import TransactionEncoder


def main(
    inputs,
    infile,
    outfile,
    min_support=0.5,
    min_confidence=0.5,
    min_lift=1.0,
    min_conviction=1.0,
    max_length=None,
):
    """
    Parameter
    ---------
    input : str
        File path to galaxy tool parameter

    infile : str
        File paths of input vector

    outfile : str
        File path to output matrix

    min_support: float
        Minimum support

    min_confidence: float
        Minimum confidence

    min_lift: float
        Minimum lift

    min_conviction: float
        Minimum conviction

    max_length: int
        Maximum length

    """
    warnings.simplefilter("ignore")

    with open(inputs, "r") as param_handler:
        params = json.load(param_handler)

    input_header = params["header0"]
    header = "infer" if input_header else None

    with open(infile) as fp:
        lines = fp.read().splitlines()

    if header is not None:
        lines = lines[1:]

    dataset = []
    for line in lines:
        line_items = line.split("\t")
        dataset.append(line_items)

    # TransactionEncoder learns the unique labels in the dataset and transforms the
    # input dataset (a Python list of lists) into a one-hot encoded NumPy boolean array
    te = TransactionEncoder()
    te_ary = te.fit_transform(dataset)

    # Turn the encoded NumPy array into a DataFrame
    df = pd.DataFrame(te_ary, columns=te.columns_)

    # Extract frequent itemsets for association rule mining
    # use_colnames: Use DataFrames' column names in the returned DataFrame instead of column indices
    frequent_itemsets = fpgrowth(
        df, min_support=min_support, use_colnames=True, max_len=max_length
    )

    # Get association rules, with confidence larger than min_confidence
    rules = association_rules(
        frequent_itemsets, metric="confidence", min_threshold=min_confidence
    )

    # Filter association rules, keeping rules with lift and conviction larger than min_liftand and min_conviction
    rules = rules[(rules["lift"] >= min_lift) & (rules["conviction"] >= min_conviction)]

    # Convert columns from frozenset to list (more readable)
    rules["antecedents"] = rules["antecedents"].apply(list)
    rules["consequents"] = rules["consequents"].apply(list)

    # The next 3 steps are intended to fix the order of the association
    # rules generated, so tests that rely on diff'ing a desired output
    # with an expected output can pass

    # 1) Sort entry in every row/column for columns 'antecedents' and 'consequents'
    rules["antecedents"] = rules["antecedents"].apply(lambda row: sorted(row))
    rules["consequents"] = rules["consequents"].apply(lambda row: sorted(row))

    # 2) Create two temporary string columns to sort on
    rules["ant_str"] = rules["antecedents"].apply(lambda row: " ".join(row))
    rules["con_str"] = rules["consequents"].apply(lambda row: " ".join(row))

    # 3) Sort results so they are re-producable
    rules.sort_values(by=["ant_str", "con_str"], inplace=True)
    del rules["ant_str"]
    del rules["con_str"]
    rules.reset_index(drop=True, inplace=True)

    # Write association rules and metrics to file
    rules.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-y", "--infile", dest="infile", required=True)
    aparser.add_argument("-o", "--outfile", dest="outfile", required=True)
    aparser.add_argument("-s", "--support", dest="support", default=0.5)
    aparser.add_argument("-c", "--confidence", dest="confidence", default=0.5)
    aparser.add_argument("-l", "--lift", dest="lift", default=1.0)
    aparser.add_argument("-v", "--conviction", dest="conviction", default=1.0)
    aparser.add_argument("-t", "--length", dest="length", default=5)
    args = aparser.parse_args()

    main(
        args.inputs,
        args.infile,
        args.outfile,
        min_support=float(args.support),
        min_confidence=float(args.confidence),
        min_lift=float(args.lift),
        min_conviction=float(args.conviction),
        max_length=int(args.length),
    )
