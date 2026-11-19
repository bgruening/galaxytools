#!/usr/bin/env python
"""Write the ELN rules and instructions of pynxtools-spm, for the ELN tools.

They come from the installed package (``pynxtools_spm/galaxy/``), or from an override file, e.g. a
fixed version given by the NOMAD Action. An override of the rules must be YAML with `definitions`.
"""

import argparse
import logging
import shutil
import sys

import yaml

logger = logging.getLogger("eln_files")


def package_files():
    """(rules, instructions) paths of the installed pynxtools-spm; stops if it has none."""
    try:
        from pynxtools_spm.galaxy import ELN_INSTRUCTIONS_FILE, ELN_RULES_FILE
    except ImportError as error:
        sys.exit(f"The installed pynxtools-spm has no ELN files (pynxtools_spm.galaxy): {error}")
    return ELN_RULES_FILE, ELN_INSTRUCTIONS_FILE


def check_rules(path):
    """Stop unless `path` is YAML with a `definitions` mapping."""
    try:
        with open(path, encoding="utf-8") as handle:
            rules = yaml.safe_load(handle)
    except (OSError, UnicodeDecodeError, yaml.YAMLError) as error:
        sys.exit(f"ELN rules not readable: {error}")
    if not isinstance(rules, dict) or not isinstance(rules.get("definitions"), dict):
        sys.exit("ELN rules have no 'definitions' mapping")


def main():
    parser = argparse.ArgumentParser(description="Write the ELN rules and instructions of pynxtools-spm.")
    parser.add_argument("--rules-override", help="Rules YAML to use instead of the package's.")
    parser.add_argument("--instructions-override", help="Instructions to use instead of the package's.")
    parser.add_argument("--rules", required=True, help="Output rules YAML.")
    parser.add_argument("--instructions", required=True, help="Output instructions.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr)])

    rules, instructions = args.rules_override, args.instructions_override
    if not (rules and instructions):
        package_rules, package_instructions = package_files()
        rules, instructions = rules or package_rules, instructions or package_instructions
    check_rules(rules)
    shutil.copyfile(rules, args.rules)
    shutil.copyfile(instructions, args.instructions)
    logger.info("rules: %s", "override" if args.rules_override else rules)
    logger.info("instructions: %s", "override" if args.instructions_override else instructions)


if __name__ == "__main__":
    main()
