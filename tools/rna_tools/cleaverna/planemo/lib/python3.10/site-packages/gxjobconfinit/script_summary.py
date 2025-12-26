#!/usr/bin/env python

# how to use this from venv
# PYTHONPATH=lib python gxjobconfinit/script_summary.py
import argparse

from gxjobconfinit.generate import summary_markdown

DESCRIPTION = """Generate summary of galaxy-job-config-init example outputs.

These will be generated in markdown format along with the CLI that would generate
them. It would be pretty impossible to test all possible combination but this script
can be used to visually inspect possible generated configurations for manual review.
"""


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.parse_args()

    markdown_output = summary_markdown()
    print(markdown_output)


if __name__ == "__main__":
    main()
