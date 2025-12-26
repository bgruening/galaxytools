# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import errno
import json
import os
import sys

from spython.logger import bot

################################################################################
## FOLDER OPERATIONS ###########################################################
################################################################################


def mkdir_p(path):
    """mkdir_p attempts to get the same functionality as mkdir -p
    :param path: the path to create.
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            bot.error("Error creating path %s, exiting." % path)
            sys.exit(1)


################################################################################
## FILE OPERATIONS #############################################################
################################################################################


def write_file(filename, content, mode="w"):
    """write_file will open a file, "filename" and write content, "content"
    and properly close the file
    """
    with open(filename, mode) as filey:
        filey.writelines(content)
    return filename


def write_json(json_obj, filename, mode="w", print_pretty=True):
    """
    write_json will (optionally,pretty print) a json object to file
    :param json_obj: the dict to print to json
    :param filename: the output file to write to
    :param pretty_print: if True, will use nicer formatting
    """
    with open(filename, mode) as filey:
        if print_pretty:
            filey.writelines(json.dumps(json_obj, indent=4, separators=(",", ": ")))
        else:
            filey.writelines(json.dumps(json_obj))
    return filename


def read_file(filename, mode="r", readlines=True):
    """
    write_file will open a file, "filename" and write content, "content"
    and properly close the file
    """
    with open(filename, mode) as filey:
        if readlines:
            content = filey.readlines()
        else:
            content = filey.read()
    return content


def read_json(filename, mode="r"):
    """
    read_json reads in a json file and returns
    the data structure as dict.
    """
    with open(filename, mode) as filey:
        data = json.load(filey)
    return data
