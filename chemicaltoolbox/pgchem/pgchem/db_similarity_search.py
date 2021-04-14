#!/usr/bin/env python
"""
Input:  File with a single molecule or SMILES string
Output: Similar molecules.
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
        type=str,
        choices=["pattern", "infile"],
        help="input is file containing multiple molecules as SMILES string or a unique SMILES string pattern?",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="input file name or SMILES patters",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output.sdf", help="output file name"
    )
    parser.add_argument("-oformat", required=True, type=str, help="output file format")
    parser.add_argument(
        "-tanimoto", type=float, default=0.85, help="Minimum Tanimoto threshold"
    )
    parser.add_argument(
        "-dbname",
        type=str,
        required=True,
        help="Specify the name of the db to connect to",
    )
    parser.add_argument(
        "-dbuser",
        type=str,
        required=True,
        help="Specify the user name for the connection to the db",
    )
    parser.add_argument(
        "-dbhost",
        type=str,
        required=True,
        help="Specify the host for the db connection",
    )
    parser.add_argument(
        "-dbpasswd",
        type=str,
        required=True,
        help="Specify the password for the connection to the db",
    )
    parser.add_argument(
        "-lib",
        type=str,
        required=True,
        help="Name of the libname where the compounds are stored",
    )
    parser.add_argument("-fetch", type=str, help="Specify query to perform")
    parser.add_argument(
        "--header",
        type=bool,
        help="Do you want to include the header as the first line of the output table?",
    )
    return parser.parse_args()


def create_query(args, mol):
    # the replace statement replaces \ by \\, otherwise PostGreSQL messes the string...
    return """SELECT *, mol @ E'%(mol)s'::molecule as tani
                FROM %(libname)s.tbl_molecule molec
                LEFT JOIN (SELECT fk_molecule_id, MAX(synonym) as synonym FROM %(libname)s.tbl_molecule_synonym GROUP BY fk_molecule_id) syn
                ON molec.id = syn.fk_molecule_id
                WHERE set_bits
                BETWEEN floor(nbits_set(fp2string(E'%(mol)s'::molecule))*%(Tanimoto)f)::integer
                AND ceil(nbits_set(fp2string(E'%(mol)s'::molecule))/%(Tanimoto)f)::integer
                AND mol @ E'%(mol)s'::molecule >= %(Tanimoto)f
                ORDER BY tani desc
            """ % {
        "libname": args.lib,
        "mol": mol.write("can").split()[0].strip().replace("\\", "\\\\"),
        "Tanimoto": args.tanimoto,
    }


def __main__():
    """
    Perform similarity searches on the selected db based on Tanimoto similarity of Fingerprint descriptors.
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

        rows = []
        # this crap deletes one level from the list...
        for first_row in old_rows:
            for row in first_row:
                rows.append(row)

    if args.oformat in ["table", "sdf"]:
        # Add the Tanimoto coefficient as output for the tabular and sdf selections in similarity searches
        args.fetch += ", %s" % "tani"

    cheminfolib.print_output(args, rows)


if __name__ == "__main__":
    __main__()
