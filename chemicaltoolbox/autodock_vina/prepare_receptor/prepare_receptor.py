#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Galaxy Tools Project
# Tool implementation for prepare ligands
# @authors: Léo Biscassi - leo.biscassi@gmail.com
#           Rodrigo Faccioli - rodrigo.faccioli@gmail.com
import argparse
import os
import shutil
from subprocess import Popen, PIPE


def load_parameters():
    parser = argparse.ArgumentParser(description="script for prepare receptor")
    parser.add_argument("-r", action="store", dest="receptor", required=True)
    parser.add_argument("-p", action="store",
                        dest="path_output", required=True)
    parser.add_argument("-f", action="store",
                        dest="file_output", required=True)

    arguments = parser.parse_args()

    return arguments


def load_settings():
    pythonsh = os.environ["PYTHONSH"]
    script_receptor4 = os.environ["SCRIPT_RECEPTOR4"]

    return pythonsh, script_receptor4


def check_autodock_programs(pythonsh, script_receptor4):
    if not os.path.isfile(pythonsh):
        raise EnvironmentError("*** ERROR: pythonsh not found! ***")

    if not os.path.isfile(script_receptor4):
        raise EnvironmentError("*** ERROR: script_receptor4 not found! ***")


def check_pdb_path(path_pdb):
    if not os.path.exists(path_pdb):
        os.makedirs(path_pdb)


def convert_dat_to_pdb(path_dat_file, path_pdb):
    if not os.path.isfile(path_dat_file):
        raise EnvironmentError("*** ERROR: dataset not found! ***")

    shutil.copy(path_dat_file, path_pdb)

    for data_file in os.listdir(path_pdb):
        if data_file.endswith(".dat"):
            new_file = data_file.split(".")[0] + ".pdb"
            file_name = os.path.join(path_pdb, data_file)
            new_file_name = os.path.join(path_pdb, new_file)
            os.rename(file_name, new_file_name)


def get_pdb_file(path_pdb):
    if not os.path.exists(path_pdb):
        raise EnvironmentError("*** ERROR: pdb path not found! ***")

    for file_name in os.listdir(path_pdb):
        if file_name.endswith(".pdb"):
            path_pdb_file = os.path.join(path_pdb, file_name)

    return path_pdb_file


def convert_pdb_to_pdbqt(pdb_file, file_output, pythonsh, script_receptor4):
    process = Popen([pythonsh, script_receptor4, '-r', pdb_file, '-o', file_output, '-v', '-U', 'nphs_lps', '-A', 'hydrogens'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()


def main():
    arguments = load_parameters()
    path_pdb = os.path.join(arguments.path_output, "pdb")
    pythonsh, script_receptor4 = load_settings()
    check_autodock_programs(pythonsh, script_receptor4)
    check_pdb_path(path_pdb)
    convert_dat_to_pdb(arguments.receptor, path_pdb)
    pdb_file = get_pdb_file(path_pdb)
    convert_pdb_to_pdbqt(pdb_file, arguments.file_output,
                         pythonsh, script_receptor4)

main()
