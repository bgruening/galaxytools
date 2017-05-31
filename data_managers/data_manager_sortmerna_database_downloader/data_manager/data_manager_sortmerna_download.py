#!/usr/bin/env python
# Data manager for reference data for the SortMeRNA Galaxy tools

import argparse
import json
import os
import tarfile
import requests
import subprocess


# Utility functions for interacting with Galaxy JSON
def read_input_json(jsonfile):
    """Read the JSON supplied from the data manager tool

    Returns a tuple (param_dict,extra_files_path)

    'param_dict' is an arbitrary dictionary of parameters
    input into the tool; 'extra_files_path' is the path
    to a directory where output files must be put for the
    receiving data manager to pick them up.

    NB the directory pointed to by 'extra_files_path'
    doesn't exist initially, it is the job of the script
    to create it if necessary.

    """
    params = json.loads(open(jsonfile).read())
    return (params['param_dict'],
            params['output_data'][0]['extra_files_path'])


# Utility functions for creating data table dictionaries
#
# Example usage:
# >>> d = create_data_tables_dict()
# >>> add_data_table(d,'my_data')
# >>> add_data_table_entry(dict(dbkey='hg19',value='human'))
# >>> add_data_table_entry(dict(dbkey='mm9',value='mouse'))
# >>> print str(json.dumps(d))
def create_data_tables_dict():
    """Return a dictionary for storing data table information

    Returns a dictionary that can be used with 'add_data_table'
    and 'add_data_table_entry' to store information about a
    data table. It can be converted to JSON to be sent back to
    the data manager.

    """
    d = {}
    d['data_tables'] = {}
    return d


def add_data_table(d, table):
    """Add a data table to the data tables dictionary

    Creates a placeholder for a data table called 'table'.

    """
    d['data_tables'][table] = []


def add_data_table_entry(d, table, entry):
    """Add an entry to a data table

    Appends an entry to the data table 'table'. 'entry'
    should be a dictionary where the keys are the names of
    columns in the data table.

    Raises an exception if the named data table doesn't
    exist.

    """
    try:
        d['data_tables'][table].append(entry)
    except KeyError:
        raise Exception("add_data_table_entry: no table '%s'" % table)


def download_archive(version):
    """

    """
    filepath = "%s.tar.gz" % (version)
    url = "https://github.com/biocore/sortmerna/archive/%s.tar.gz" % (version)
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(filepath, "wb") as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)
    return filepath


def find_archive_content_path(archive_content_path):
    """
    """
    content = os.listdir(archive_content_path)
    archive_content = []
    for x in content:
        if not x.startswith(".") and not x.startswith("_"):
            archive_content.append(x)
    if len(archive_content) == 1:
        archive_content_path = os.path.join(
            archive_content_path,
            archive_content[0])
    return archive_content_path


def extract_archive(filepath):
    """
    """
    archive_content_path = "tmp"
    tar = tarfile.open(filepath)
    tar.extractall(path=archive_content_path)
    tar.close()
    archive_content_path = find_archive_content_path(archive_content_path)
    return archive_content_path


def move_index_files(archive_content_path, target_dir, data_tables, version):
    """
    """
    file_dir = os.path.join(archive_content_path, "rRNA_databases")
    for filename in os.listdir(file_dir):
        if not filename.endswith("fasta"):
            continue
        input_filepath = os.path.join(file_dir, filename)
        output_filepath = os.path.join(target_dir, filename)
        # Move file
        os.rename(input_filepath, output_filepath)
        # Index the file with indexdb_rna
        command = "indexdb_rna --ref %s,%s" % (
            output_filepath,
            os.path.splitext(output_filepath)[0])
        process = subprocess.call(command, shell=True )
        # Add entry in the data table
        db_name = os.path.splitext(filename)[0]
        add_data_table_entry(
            data_tables,
            "rRNA_databases",
            dict(
                dbkey=db_name,
                value=version,
                name=db_name,
                path=output_filepath))


def download_db(data_tables, version, target_dir):
    """Download SortMeRNA database

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      version: version of the database
      table_name: name of the table
      target_dir: directory to put copy or link to the data file
    """
    print("Download archive")
    filepath = download_archive(version)

    print("Extract archive %s" % filepath)
    archive_content_path = extract_archive(filepath)

    print("Moving fasta file from %s and index them" % archive_content_path)
    move_index_files(
        archive_content_path,
        target_dir,
        data_tables,
        version)


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = argparse.ArgumentParser(
        description='Download QIIME reference database')
    parser.add_argument('--version', help="Database version")
    parser.add_argument('--jsonfile', help="Output JSON file")
    args = parser.parse_args()

    jsonfile = args.jsonfile

    # Read the input JSON
    params, target_dir = read_input_json(jsonfile)

    # Make the target directory
    print("Making %s" % target_dir)
    os.mkdir(target_dir)
    os.mkdir(os.path.join(target_dir, "rRNA_databases"))
    target_dir = os.path.join(target_dir, "rRNA_databases")

    # Set up data tables dictionary
    data_tables = create_data_tables_dict()
    add_data_table(data_tables, "rRNA_databases")

    # Fetch data from specified data sources
    download_db(
        data_tables,
        args.version,
        target_dir)

    # Write output JSON
    print("Outputting JSON")
    print(str(json.dumps(data_tables)))
    with open(jsonfile, 'w') as out:
        json.dump(data_tables, out)
    print("Done.")
