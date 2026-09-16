#!/usr/bin/env python

import argparse
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
from urllib.request import urlopen

import wget


def rmtree_tolerant(path, retries=5, delay=2):
    """Remove ``path`` recursively, tolerating the transient "Directory not
    empty" (OSError errno 39) race that ``shutil.rmtree`` hits on parallel/network
    filesystems (e.g. Lustre), where ``readdir`` can return stale entries
    mid-walk. Retries a few times, then gives up quietly: ``path`` is scratch
    space and the fetched DB has already been copied elsewhere, so a leftover
    temp dir must never fail the job."""
    for attempt in range(retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            if attempt == retries - 1:
                shutil.rmtree(path, ignore_errors=True)
                return
            time.sleep(delay)


version_mapping = {
    "3.1.0": "https://zenodo.org/records/7778108/files/db_mOTU_v3.1.0.tar.gz",
    "3.0.1": "https://zenodo.org/records/5140350/files/db_mOTU_v3.0.1.tar.gz",
    "3.0.0": "https://zenodo.org/records/5012106/files/db_mOTU_v3.0.0.tar.gz",
}


def zenodo_publication_date(file_url):
    """Return the publication date (YYYY-MM-DD) of the Zenodo record hosting
    ``file_url`` (``.../records/<record_id>/files/<name>``). Used to derive a
    deterministic default DB identifier from the record itself rather than from
    the wall-clock time of the build."""
    record_id = file_url.split("/records/")[1].split("/")[0]
    with urlopen(
        f"https://zenodo.org/api/records/{record_id}"
    ) as fh:  # noqa: S310 (fixed https host)
        record = json.load(fh)
    return record["metadata"]["publication_date"]


def download_untar_store(url, tmp_path, dest_path):
    """
    Download a tar.gz file containing one folder,
    extract that folder and move the content inside dest_path
    """

    extract_path = os.path.join(tmp_path, "extract")

    os.makedirs(tmp_path, exist_ok=True)

    # download data
    filename = wget.download(url, out=tmp_path)
    tarfile_path = os.path.join(tmp_path, filename)
    tar = tarfile.open(tarfile_path)
    tar.extractall(extract_path)

    if len(list(os.listdir(extract_path))) > 1:
        print("More then one folder in zipped file, aborting !")
    else:
        for folder in os.listdir(extract_path):
            folder_path = os.path.join(extract_path, folder)

            print(f"Move data to {dest_path}")
            # folder_path and dest_path are both under workdir (same
            # filesystem), so this is an atomic rename: no multi-GB duplicate
            # copy, and no recursive delete of the freshly-written DB tree -
            # which is what triggers the Lustre readdir/rmtree race below.
            shutil.move(folder_path, dest_path)
            print("Done !")

    rmtree_tolerant(tmp_path)


def main():
    # Parse Command Line
    parser = argparse.ArgumentParser(description="Create data manager JSON.")
    parser.add_argument("--out", dest="output", action="store", help="JSON filename")
    parser.add_argument(
        "--version", dest="version", action="store", help="Version of the DB"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="option to test the script with an lighted database",
    )
    parser.add_argument(
        "--value",
        dest="value",
        action="store",
        default=None,
        help=(
            "Override the data-table value / on-disk directory name. Defaults to "
            "db_from_<Zenodo publication date>. Set it to reproduce or match an "
            "existing entry from another server."
        ),
    )

    args = parser.parse_args()

    # the output file of a DM is a json containing args that can be used by the DM
    # most tools mainly use these args to find the extra_files_path for the DM, which can be used
    # to store the DB data
    with open(args.output) as fh:
        params = json.load(fh)

    workdir = params["output_data"][0]["extra_files_path"]
    os.mkdir(workdir)

    url = version_mapping[args.version]
    # Default the identifier to the Zenodo record's publication date so it is
    # deterministic/reproducible across builds; an explicit --value overrides it
    # (e.g. to match an existing entry already published on another server).
    if args.value:
        db_value = args.value
    else:
        db_value = "db_from_{0}".format(zenodo_publication_date(url))
    db_path = os.path.join(workdir, db_value)
    tmp_path = os.path.join(workdir, "tmp")

    # create DB
    if args.test:  # the test only checks that the pharokka download script is available

        # check if link is there
        command_args = ["wget", "--spider", url]
        proc = subprocess.Popen(args=command_args, shell=False)
        return_code = proc.wait()
        if return_code:
            print("Error downloading motus database.", file=sys.stderr)
            sys.exit(return_code)

        # copy the test DB
        # TODO ones available: https://github.com/motu-tool/mOTUs/issues/121
        test_db_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "motus_test_DB_non_functional"
        )
        command_args = ["cp", "-r", test_db_path, db_path]
        proc = subprocess.Popen(args=command_args, shell=False)
        return_code = proc.wait()
        if return_code:
            print("Error copying motus database.", file=sys.stderr)
            sys.exit(return_code)

    else:

        # download data
        download_untar_store(url, tmp_path, db_path)

    # Update Data Manager JSON and write to file
    data_manager_entry = {
        "data_tables": {
            "motus_db_versioned": {
                "value": db_value,
                "version": args.version,
                "name": f"mOTUs DB version {args.version} ({db_value})",
                "path": db_path,
            }
        }
    }

    with open(os.path.join(args.output), "w+") as fh:
        json.dump(data_manager_entry, fh, sort_keys=True)


if __name__ == "__main__":
    main()
