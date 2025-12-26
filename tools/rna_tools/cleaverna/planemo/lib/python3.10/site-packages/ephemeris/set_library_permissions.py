#!/usr/bin/env python
"""Tool to set permissions for all datasets of a given Galaxy Data Library"""

import argparse
import logging as log
import sys

from bioblend import galaxy
from rich.progress import Progress

from .common_parser import (
    get_common_args,
    HideUnderscoresHelpFormatter,
)

# Print iterations progress


def get_datasets(gi, library_id) -> list[str]:
    objects = gi.libraries.show_dataset(library_id=library_id, dataset_id="")
    datasets = []
    for index in range(len(objects)):
        if objects[index]["type"] == "file":
            datasets.append(objects[index]["id"])
    if datasets == []:
        sys.exit("No datasets in library!")
    else:
        return datasets


def set_permissions(gi, library_id, role_ids, auto):
    log.info("Your library_id is %s", library_id)
    log.info("Your roles are: %s", " ".join(role_ids))
    datasets = get_datasets(gi, library_id)
    total = len(datasets)
    est = total * 3 / 60
    # Give User time to abort
    log.info(
        "\nSuccess! %d datasets found. Processing can take up to %0.02f min\n",
        total,
        est,
    )
    if auto:
        for current in range(total):
            log.debug("Processing dataset %d of %d, ID=%s", current, total, datasets[current])
            gi.libraries.set_dataset_permissions(
                dataset_id=datasets[current],
                access_in=role_ids,
                modify_in=role_ids,
                manage_in=role_ids,
            )
    else:
        if input("Do you want to continue? (y/n) ") == "y":
            with Progress() as progress:
                task = progress.add_task("[green]Processing datasets...", total=total)
                for current in range(total):
                    log.debug(
                        "Processing dataset %d of %d, ID=%s",
                        current,
                        total,
                        datasets[current],
                    )
                    gi.libraries.set_dataset_permissions(
                        dataset_id=datasets[current],
                        access_in=role_ids,
                        modify_in=role_ids,
                        manage_in=role_ids,
                    )
                    progress.update(task, advance=1)
        else:
            log.info("Operation cancelled by user. No changes were applied.\n")


def _parser():
    """Constructs the parser object"""
    parent = get_common_args()
    parser = argparse.ArgumentParser(
        parents=[parent],
        formatter_class=HideUnderscoresHelpFormatter,
        description="Populate the Galaxy data library with data.",
    )
    parser.add_argument("library", help="Specify the data library ID")
    parser.add_argument("--roles", nargs="+", help="Specify a list of comma separated role IDs")
    parser.add_argument(
        "-y",
        "--yes",
        default=False,
        action="store_true",
        help="Set the -y flag for auto-accept and skip manual approvement",
    )
    parser.add_argument(
        "-s",
        "--silent",
        default=False,
        action="store_true",
        help="sets loglevel to ERROR",
    )
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.user and args.password:
        gi = galaxy.GalaxyInstance(url=args.galaxy, email=args.user, password=args.password)
    elif args.api_key:
        gi = galaxy.GalaxyInstance(url=args.galaxy, key=args.api_key)
    else:
        sys.exit("Please specify either a valid Galaxy username/password or an API key.")

    if args.verbose:
        log.basicConfig(level=log.DEBUG)
    elif args.silent:
        log.basicConfig(level=log.ERROR)
    else:
        log.basicConfig(level=log.INFO)

    if args.roles and args.library:
        args.roles = [r.strip() for r in args.roles.split(",")]
    else:
        sys.exit("Specify library ID (--library myLibraryID) and (list of) role(s) (--roles roleId1,roleId2)")
    set_permissions(gi, library_id=args.library, role_ids=args.roles, auto=args.yes)
    log.info(
        "\nThis script uses bioblend to update ALL permissions of ALL datasets in a"
        "specified library to the given roles. Be careful and cancel if unsure\n"
    )


if __name__ == "__main__":
    main()
