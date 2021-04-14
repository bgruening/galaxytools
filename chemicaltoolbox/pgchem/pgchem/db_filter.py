#!/usr/bin/env python
"""
Input:  User-defined set of filters and the database authentifications.
Output: set of compounds that pass all the selected filters.
"""
import sys, os
import argparse
import psycopg2.extras
import cheminfolib
import json

cheminfolib.pybel_stop_logging()


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
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
    parser.add_argument("--filters", help="Specify the filters to be applied")
    return parser.parse_args()


def __main__():
    """
    Filter a db using user-selected physico-chemical properties
    """

    args = parse_command_line(sys.argv)
    db_conn = cheminfolib.db_connect(args)
    filters = json.loads((args.filters).replace(" ", "").replace(",}", "}"))

    query = """SELECT can_smiles, synonym
                FROM %(libname)s.tbl_molecule molec
                LEFT JOIN (SELECT fk_molecule_id, MAX(synonym) as synonym FROM %(libname)s.tbl_molecule_synonym GROUP BY fk_molecule_id) syn
                ON molec.id = syn.fk_molecule_id
                WHERE %(filter_query)s
             """ % {
        "libname": args.lib,
        "filter_query": " AND ".join(
            [
                ("%s BETWEEN %s AND %s" % (key, elem[0], elem[1]))
                for key, elem in filters.items()
            ]
        ),
    }

    cur = db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(query)
    rows = cur.fetchall()

    # Only fetch smiles from the database since SDF retrieval takes ages
    args.oformat = "can_smiles"

    cheminfolib.print_output(args, rows)


if __name__ == "__main__":
    __main__()
