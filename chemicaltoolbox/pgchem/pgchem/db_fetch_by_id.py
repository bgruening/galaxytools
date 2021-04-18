#!/usr/bin/env python
"""
    Input:  File with tab separated columns.
            One column contains the id's of a collection of molecules
            that should to be retrieved.
    Output: Retrieved compounds
"""
import argparse
import sys

import cheminfolib
import psycopg2.extras

cheminfolib.pybel_stop_logging()


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="input file name containing one or multiple unique molecule identifiers",
    )
    parser.add_argument(
        "-c",
        "--column",
        required=True,
        type=int,
        help="#column containing the id codes",
    )
    parser.add_argument("-oformat", default="sdf", help="output file format")
    parser.add_argument("-o", "--output", help="output file name")
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


def __main__():
    """
    Query the database and retrieve the requested compounds. A unique identifier for the molecules must be provided in a tab file
    """
    args = parse_command_line(sys.argv)
    identifiers = [
        line.split("\t")[args.column - 1].strip() for line in open(args.input, "r")
    ]

    db_conn = cheminfolib.db_connect(args)

    query = """SELECT *
                FROM %(libname)s.tbl_molecule molec
                LEFT JOIN %(libname)s.tbl_molecule_synonym syn
                ON molec.id = syn.fk_molecule_id
                WHERE syn.synonym IN ('%(id)s')
            """ % {
        "libname": args.lib,
        "id": ("', '").join(identifiers),
    }

    cur = db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(query)
    rows = cur.fetchall()

    cheminfolib.print_output(args, rows)


if __name__ == "__main__":
    __main__()
