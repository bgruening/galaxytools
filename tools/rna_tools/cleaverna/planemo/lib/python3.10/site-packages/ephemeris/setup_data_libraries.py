#!/usr/bin/env python
"""Tool to setup data libraries on a galaxy instance"""
import argparse
import logging as log
import sys
import time

import yaml
from bioblend import galaxy

from .common_parser import (
    DEFAULT_JOB_SLEEP,
    get_common_args,
    HideUnderscoresHelpFormatter,
)


def create_legacy(gi, desc):
    destination = desc["destination"]
    if destination["type"] != "library":
        raise Exception("Only libraries may be created with pre-18.05 Galaxies using this script.")
    library_name = destination.get("name")
    library_description = destination.get("description")
    library_synopsis = destination.get("synopsis")

    # Check to see if the library already exists. If it does, do not recreate it. If it doesn't, create it.
    lib_id = None
    print("Library name: " + str(library_name))
    rmt_lib_list = gi.libraries.get_libraries(name=library_name, deleted=False)
    # Now we need to check if the library has been deleted since deleted=False still returns the deleted libraries!
    not_deleted_rmt_lib_list = []
    folder_id = None

    if rmt_lib_list:
        for x in rmt_lib_list:
            if not x["deleted"]:
                not_deleted_rmt_lib_list.append(x)
    if not_deleted_rmt_lib_list:
        lib_id = not_deleted_rmt_lib_list[0]["id"]
        print("Library already exists! id: " + str(lib_id))
        folder_id = gi.libraries.show_library(lib_id)["root_folder_id"]
    else:
        lib = gi.libraries.create_library(library_name, library_description, library_synopsis)
        lib_id = lib["id"]
        folder_id = lib["root_folder_id"]

    def populate_items(base_folder_id, has_items):
        if "items" in has_items:
            name = has_items.get("name")
            description = has_items.get("description")
            folder_id = base_folder_id
            if name:
                # Check to see if the folder already exists, if it doesn't create it.
                rmt_folder_list = []
                folder = gi.libraries.get_folders(lib_id)
                new_folder_name = "/" + name
                if folder and not folder[0]["name"] == "/":
                    new_folder_name = folder[0]["name"] + "/" + name
                rmt_folder_list = gi.libraries.get_folders(lib_id, name=new_folder_name)
                if rmt_folder_list:
                    folder_id = rmt_folder_list[0]["id"]
                else:
                    folder = gi.libraries.create_folder(lib_id, name, description, base_folder_id=base_folder_id)
                    folder_id = folder[0]["id"]
            for item in has_items["items"]:
                populate_items(folder_id, item)
        else:
            src = has_items["src"]
            if src != "url":
                raise Exception("For pre-18.05 Galaxies only support URLs src items are supported.")
            rmt_library_files = gi.folders.show_folder(base_folder_id, contents=True)["folder_contents"]
            file_names = []
            for item in rmt_library_files:
                if item["type"] == "file":
                    file_names.append(item["name"])
            if has_items["url"] not in file_names:
                try:
                    gi.libraries.upload_file_from_url(
                        lib_id,
                        has_items["url"],
                        folder_id=base_folder_id,
                        file_type=has_items["ext"],
                    )
                except Exception:
                    log.exception(
                        "Could not upload %s to %s/%s",
                        has_items["url"],
                        lib_id,
                        base_folder_id,
                    )
        return None

    populate_items(folder_id, desc)
    return []


def create_batch_api(gi, desc):
    hc = galaxy.histories.HistoryClient(gi)
    tc = galaxy.tools.ToolClient(gi)

    history = hc.create_history()
    url = f"{gi.url}/tools/fetch"
    payload = {"targets": [desc], "history_id": history["id"]}
    yield tc._post(payload=payload, url=url)


def setup_data_libraries(gi, data, training=False, legacy=False):
    """
    Load files into a Galaxy data library.
    By default all test-data tools from all installed tools
    will be linked into a data library.
    """

    log.info("Importing data libraries.")
    jc = galaxy.jobs.JobsClient(gi)
    config = galaxy.config.ConfigClient(gi)
    version = config.get_version()

    if legacy:
        create_func = create_legacy
    else:
        version_major = version.get("version_major", "16.01")
        create_func = create_batch_api if version_major >= "18.05" else create_legacy

    library_def = yaml.safe_load(data)

    def normalize_items(has_items):
        # Synchronize Galaxy batch format with older training material style.
        if "files" in has_items:
            items = has_items.pop("files")
            has_items["items"] = items

        items = has_items.get("items", [])
        for item in items:
            normalize_items(item)
            src = item.get("src")
            url = item.get("url")
            if src is None and url:
                item["src"] = "url"
            if "file_type" in item:
                ext = item.pop("file_type")
                item["ext"] = ext

    # Normalize library definitions to allow older ephemeris style and native Galaxy batch
    # upload formats.
    if "libraries" in library_def:
        # File contains multiple definitions.
        library_def["items"] = library_def.pop("libraries")

    if "destination" not in library_def:
        library_def["destination"] = {"type": "library"}
    destination = library_def["destination"]

    if training:
        destination["name"] = destination.get("name", "Training Data")
        destination["description"] = destination.get("description", "Data pulled from online archives.")
    else:
        destination["name"] = destination.get("name", "New Data Library")
        destination["description"] = destination.get("description", "")

    normalize_items(library_def)

    if library_def:
        jobs = list(create_func(gi, library_def))

        job_ids = []
        if legacy:
            for job in jc.get_jobs():
                # Fetch all upload job IDs, ignoring complete ones.
                if job["tool_id"] == "upload1" and job["state"] not in ("ok", "error"):
                    job_ids.append(job["id"])

            # Just have to check that all upload1 jobs are termianl.
        else:
            # Otherwise get back an actual list of jobs
            for job in jobs:
                if "jobs" in job:
                    for subjob in job["jobs"]:
                        job_ids.append(subjob["id"])

        while True:
            job_states = [jc.get_state(job) in ("ok", "error", "deleted") for job in job_ids]
            log.debug(
                f'Job states: {",".join([f"{job_id}={job_state}" for (job_id, job_state) in zip(job_ids, job_states)])}'
            )

            if all(job_states):
                break
            time.sleep(DEFAULT_JOB_SLEEP)

        log.info("Finished importing test data.")


def _parser():
    """Constructs the parser object"""
    parent = get_common_args()
    parser = argparse.ArgumentParser(
        parents=[parent],
        formatter_class=HideUnderscoresHelpFormatter,
        description="Populate the Galaxy data library with data.",
    )
    parser.add_argument("-i", "--infile", required=True, type=argparse.FileType("r"))
    parser.add_argument(
        "--training",
        default=False,
        action="store_true",
        help="Set defaults that make sense for training data.",
    )
    parser.add_argument(
        "--legacy",
        default=False,
        action="store_true",
        help="Use legacy APIs even for newer Galaxies that should have a batch upload API enabled.",
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

    setup_data_libraries(gi, args.infile, training=args.training, legacy=args.legacy)


if __name__ == "__main__":
    main()
