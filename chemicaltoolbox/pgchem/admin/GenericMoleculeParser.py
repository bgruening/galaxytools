#!/usr/bin/env python

__author__ = "Bjoern Gruening"
__version__ = "0.1"
__date__ = "2012"
__license__ = "GLP3+"

import argparse
import os
import sys
from multiprocessing import Pool

import cheminfolib
import GenericMoleculeDB
import pybel
from GenericMoleculeDB import build_database, create_indices
from sqlalchemy.orm import sessionmaker


class GenericMoleculeParser:
    def __init__(self, filepath, schema, conn_string, mol_type="sdf"):
        self.engine = GenericMoleculeDB.init(conn_string, schema)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.filepath = filepath
        self.molecule_type = mol_type
        self.metadata_synonym = ""

    def __del__(self):
        self.session.close()
        self.engine.dispose()

    def parse_into_database(self):
        for mol in pybel.readfile(self.molecule_type, self.filepath):
            DBMolecule = GenericMoleculeDB.Molecule(mol)
            synonyms = {}
            if isinstance(self.metadata_synonym, list):
                for key_name in self.metadata_synonym:
                    syn = mol.data.get(key_name, False)
                    if syn:
                        synonyms.update({syn.strip(): key_name})
            else:
                if self.metadata_synonym:
                    synonyms.update({self.metadata_synonym.strip(): ""})

            if mol.title.strip():
                synonyms.update({mol.title.strip(): ""})

            syn_obj_list = list()
            [
                syn_obj_list.append(GenericMoleculeDB.Synonym(key, value))
                for key, value in synonyms.items()
            ]

            DBMolecule.synonyms = (
                syn_obj_list  # map(GenericMoleculeParser.Synonym, synonyms.items())
            )

            try:
                self.session.add(DBMolecule)
                self.session.commit()
            except Exception:
                self.session.rollback()
                sys.stderr.write(mol.write("smi"))

        self.session.add(GenericMoleculeDB.File(self.filepath))
        try:
            self.session.commit()
        except Exception:
            self.session.rollback()
            sys.stderr.write(mol.write("smi"))
        self.session.close()
        self.engine.dispose()
        return True


def start_parser(args):
    (path, schema, conn_string, molecule_type) = args
    cheminfolib.pybel_stop_logging()
    gmp = GenericMoleculeParser(path, schema, conn_string, molecule_type)
    gmp.parse_into_database()
    return path


def run(compound_path, procs, schema, conn_string, filetype):
    clean_up_files = False
    if os.path.isdir(compound_path):
        paths = []
        for root, dirs, files in os.walk(compound_path):
            for compound_file in files:
                path = os.path.join(root, compound_file)
                paths.append((path, schema, conn_string, filetype))
        paths.sort()
    else:
        if filetype in ["smi", "inchi"]:
            clean_up_files = True
            compound_count = cheminfolib.CountLines(compound_path)
            sys.stdout.write(
                "Splitting inputfile (%s) with %s molecules in files with %s molecules.\n"
                % (compound_path, compound_count, int(compound_count / procs))
            )
            paths = cheminfolib.split_smi_library(
                compound_path, structures_in_one_file=int(compound_count / procs)
            )
            paths = [(path, schema, conn_string, filetype) for path in paths]
        if filetype == "sdf":
            paths = [compound_path]

    pool = Pool(processes=procs)
    sys.stdout.write("Process initialized with %s processes.\n" % procs)
    result = pool.map_async(start_parser, paths)
    result.get()
    if clean_up_files:
        for path in paths:
            os.remove(path[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Molecule Parser and pgchem insertion script."
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    parser.add_argument(
        "--processes",
        "-p",
        dest="processes",
        type=int,
        default=4,
        help="# of processes to use for file parsing",
    )
    parser.add_argument("--filetype", "-f", dest="filetype", default="smi")
    parser.add_argument(
        "--molecule-file",
        "-m",
        dest="compound_file",
        required=True,
        help="Compound file to load in the database.",
    )
    options = parser.parse_args()

    (schema, conn_string) = build_database()
    run(options.compound_file, options.processes, schema, conn_string, options.filetype)
    sys.stdout.write("Creating molecule indices.")
    create_indices(conn_string, schema)
