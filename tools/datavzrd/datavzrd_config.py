#!/usr/bin/env python3
"""Helper script to adapt datavzrd report configuration files for Galaxy.

Subcommands:
    rewrite: rewrite the dataset paths of a config file to point to the
             Galaxy input files. Datasets are matched to the input files
             in the order in which they appear in the config file.
    set:     set one or more top-level key-value pairs in a config file.
"""

import argparse
import sys

import yaml


def load_config(path):
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def dump_config(config, path):
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(
            config,
            fh,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )


def command_rewrite(args):
    config = load_config(args.config)
    inputs = [value for value in args.inputs.split(",") if value]
    datasets = config.get("datasets") or {}
    if len(datasets) != len(inputs):
        sys.exit(
            "Error: the config file defines %d dataset(s) (%s), but %d input "
            "file(s) were provided. Please provide exactly one input file "
            "per dataset, in the order in which the datasets appear in the "
            "config file." % (len(datasets), ", ".join(datasets), len(inputs))
        )
    for (name, dataset), path in zip(datasets.items(), inputs):
        dataset["path"] = path
    dump_config(config, args.output)


def command_set(args):
    config = load_config(args.config)
    for assignment in args.assignments:
        if "=" not in assignment:
            sys.exit("Error: invalid assignment '%s', expected key=value." % assignment)
        key, value = assignment.split("=", 1)
        config[key] = yaml.safe_load(value)
    dump_config(config, args.output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_rewrite = subparsers.add_parser(
        "rewrite", help="rewrite dataset paths to point to Galaxy input files"
    )
    parser_rewrite.add_argument("config", help="path of the datavzrd config file")
    parser_rewrite.add_argument(
        "--inputs", required=True, help="comma-separated list of input files"
    )
    parser_rewrite.add_argument(
        "--output", required=True, help="output path of the rewritten config file"
    )
    parser_rewrite.set_defaults(func=command_rewrite)

    parser_set = subparsers.add_parser(
        "set", help="set top-level key-value pairs in a config file"
    )
    parser_set.add_argument("config", help="path of the datavzrd config file")
    parser_set.add_argument(
        "--output", required=True, help="output path of the updated config file"
    )
    parser_set.add_argument(
        "assignments",
        nargs="+",
        help="key=value pairs; values are parsed as YAML scalars",
    )
    parser_set.set_defaults(func=command_set)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
