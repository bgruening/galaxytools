#!/usr/bin/env python
"""
    Input:  Molecular pattern of molecule.
    Output: Multiple substructurally related molecules file.
"""
import argparse
import sys

import cheminfolib
import psycopg2.extras
import pybel

cheminfolib.pybel_stop_logging()


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--istyle",
        required=True,
        choices=["pattern", "infile"],
        help="input is file containing multiple molecules as SMILES string or a unique SMILES string pattern?",
    )
    parser.add_argument(
        "-i", "--input", required=True, help="input file name or SMILES patters"
    )
    parser.add_argument("-o", "--output", default="output.sdf", help="output file name")
    parser.add_argument("-oformat", default="smi", help="output file format")
    parser.add_argument(
        "-max", type=int, default=10000, help="Maximum number of retrieved molecules"
    )
    parser.add_argument(
        "-dbname", required=True, help="Specify the name of the db to connect to"
    )
    parser.add_argument(
        "-dbuser",
        required=True,
        help="Specify the user name for the connection to the db",
    )
    parser.add_argument(
        "-dbhost", required=True, help="Specify the host for the db connection"
    )
    parser.add_argument(
        "-dbpasswd",
        required=True,
        help="Specify the password for the connection to the db",
    )
    parser.add_argument(
        "-lib", required=True, help="Name of the libname where the compounds are stored"
    )
    parser.add_argument("-fetch", help="Specify query to perform")
    parser.add_argument(
        "--header",
        type=bool,
        help="Do you want to include the header as the first line of the output table?",
    )
    return parser.parse_args()


def create_query(args, mol):
    """
    Create a substructure search query using pgchem.
    Replace \ by \\, otherwise PostGreSQL changes the canonical SMILES.
    """
    return """SELECT *
                FROM %(libname)s.tbl_molecule molec
                LEFT JOIN (SELECT fk_molecule_id, MAX(synonym) as synonym FROM %(libname)s.tbl_molecule_synonym GROUP BY fk_molecule_id) syn
                ON molec.id = syn.fk_molecule_id
                WHERE E'%(mol)s'::molecule <= mol
                LIMIT %(max_results)i
            """ % {
        "libname": args.lib,
        "mol": mol.write("can").split()[0].strip().replace("\\", "\\\\"),
        "max_results": args.max,
    }


def __main__():
    """
    Perform substructure searches on the selected db based on Fingerprint descriptors.
    """
    args = parse_command_line(sys.argv)
    db_conn = cheminfolib.db_connect(args)

    if args.istyle == "pattern":
        mol = pybel.readstring("smi", args.input)
        cur = db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cur.execute(create_query(args, mol))
        rows = cur.fetchall()
    else:
        old_rows = []
        for mol in pybel.readfile("smi", args.input):
            cur = db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cur.execute(create_query(args, mol))
            old_rows.append(cur.fetchall())

        ### TODO: if i have a testing env, try to not append the cur.fetchall ... + it and save the one loop here
        rows = []
        for first_row in old_rows:
            for row in first_row:
                rows.append(row)

    cheminfolib.print_output(args, rows)


if __name__ == "__main__":
    __main__()
