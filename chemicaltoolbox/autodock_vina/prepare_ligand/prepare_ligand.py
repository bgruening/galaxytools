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
    parser = argparse.ArgumentParser(description="script for prepare ligands")
    parser.add_argument("-l", action="store", dest="ligand", required=True)
    parser.add_argument("-p", action="store",
                        dest="path_output", required=True)
    parser.add_argument("-f", action="store",
                        dest="file_output", required=True)

    arguments = parser.parse_args()

    return arguments


def load_settings():
    pythonsh = os.environ["PYTHONSH"]
    script_ligand4 = os.environ["SCRIPT_LIGAND4"]

    return pythonsh, script_ligand4


def check_autodock_programs(pythonsh, script_ligand4):
    if not os.path.isfile(pythonsh):
        raise EnvironmentError("*** ERROR: pythonsh not found! ***")

    if not os.path.isfile(script_ligand4):
        raise EnvironmentError("*** ERROR: script_ligand4 not found! ***")


def check_mol2_path(path_mol2):
    if not os.path.exists(path_mol2):
        os.makedirs(path_mol2)


def convert_dat_for_mol2(path_dat_file, path_mol2):
    if not os.path.isfile(path_dat_file):
        raise EnvironmentError("*** ERROR: dataset not found! ***")

    shutil.copy(path_dat_file, path_mol2)

    for data_file in os.listdir(path_mol2):
        if data_file.endswith(".dat"):
            new_file = data_file.split(".")[0] + ".mol2"
            file_name = os.path.join(path_mol2, data_file)
            new_file_name = os.path.join(path_mol2, new_file)
            os.rename(file_name, new_file_name)


def get_mol2_file(path_mol2):
    if not os.path.exists(path_mol2):
        raise EnvironmentError("*** ERROR: mol2 path not found! ***")

    for file_name in os.listdir(path_mol2):
        if file_name.endswith(".mol2"):
            path_mol2_file = os.path.join(path_mol2, file_name)

    return path_mol2_file


def convert_mol2_to_pdbqt(mol2_file, file_output, pythonsh, script_ligand4):
    process = Popen([pythonsh, script_ligand4, '-l', mol2_file, '-v', '-o', file_output, '-U', 'nphs_lps', '-A', 'hydrogens'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()


def main():
    arguments = load_parameters()
    path_mol2 = os.path.join(arguments.path_output, "mol2")
    pythonsh, script_ligand4 = load_settings()
    check_autodock_programs(pythonsh, script_ligand4)
    check_mol2_path(path_mol2)
    convert_dat_for_mol2(arguments.ligand, path_mol2)
    mol2_file = get_mol2_file(path_mol2)
    convert_mol2_to_pdbqt(mol2_file, arguments.file_output,
                          pythonsh, script_ligand4)

main()