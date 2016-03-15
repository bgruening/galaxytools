#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Galaxy Tools Project
# Tool implementation for virtual screening
# @authors: Léo Biscassi - leo.biscassi@gmail.com
#           Rodrigo Faccioli - rodrigo.faccioli@gmail.com

import argparse
import os
import shutil
from subprocess import Popen, PIPE


def load_parameters():
    parser = argparse.ArgumentParser(description="script for run docking")
    # Parameter for ligand
    parser.add_argument("-l", action="store", dest="ligand", required=True)
    # Parameter for receptor
    parser.add_argument("-r", action="store", dest="receptor", required=True)
    # Parameter for box of active site
    parser.add_argument("-b", action="store", dest="box", required=True)
    # Parameter for path output
    parser.add_argument("-p", action="store",
                        dest="path_output", required=True)
    # Parameter for structure file
    parser.add_argument("-f", action="store",
                        dest="file_output1", required=True)
    # Parameter for log file
    parser.add_argument("-g", action="store",
                        dest="file_output2", required=True)

    arguments = parser.parse_args()

    return arguments


def load_settings():
    vina = os.environ["VINA"]
    return vina


def check_autodock_vina(autodock_vina):
    if not os.path.isfile(autodock_vina):
        raise EnvironmentError("*** ERROR: Autodock Vina not found! ***")


def check_pdbqt_path(path_pdbqt):
    if not os.path.exists(path_pdbqt):
        os.makedirs(path_pdbqt)


def convert_dat_for_pdbqt(path_dat_file, path_pdbqt):
    if not os.path.isfile(path_dat_file):
        raise EnvironmentError("*** ERROR: dataset not found! ***")

    shutil.copy(path_dat_file, path_pdbqt)

    for data_file in os.listdir(path_pdbqt):
        if data_file.endswith(".dat"):
            new_file = data_file.split(".")[0] + ".pdbqt"
            file_name = os.path.join(path_pdbqt, data_file)
            new_file_name = os.path.join(path_pdbqt, new_file)
            os.rename(file_name, new_file_name)


def get_pdbqt_file(path_pdbqt):
    if not os.path.exists(path_pdbqt):
        raise EnvironmentError("*** ERROR: path pdbqt not found! ***")

    for file_name in os.listdir(path_pdbqt):
        if file_name.endswith(".pdbqt"):
            path_pdbqt_file = os.path.join(path_pdbqt, file_name)

    return path_pdbqt_file


def run_docking(autodock_vina, config_box, ligand, receptor, structure_file, log_file):
    process = Popen([autodock_vina, '--config', config_box, '--receptor', receptor, '--ligand', ligand, '--out', structure_file, '--log', log_file], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()


def main():
    arguments = load_parameters()
    path_ligand_pdbqt = os.path.join(arguments.path_output, "ligand_pdbqt")
    path_receptor_pdbqt = os.path.join(arguments.path_output, "receptor_pdbqt")
    autodock_vina = load_settings()
    check_autodock_vina(autodock_vina)
    check_pdbqt_path(path_ligand_pdbqt)
    check_pdbqt_path(path_receptor_pdbqt)
    convert_dat_for_pdbqt(arguments.ligand, path_ligand_pdbqt)
    convert_dat_for_pdbqt(arguments.receptor, path_receptor_pdbqt)
    ligand_pdbqt = get_pdbqt_file(path_ligand_pdbqt)
    receptor_pdbqt = get_pdbqt_file(path_receptor_pdbqt)
    run_docking(
        autodock_vina,
        arguments.box,
        ligand_pdbqt,
        receptor_pdbqt,
        arguments.file_output1,
        arguments.file_output2)

main()
