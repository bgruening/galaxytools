#!/usr/bin/env python
"""
Description:
    Convert a gff file into a new database or add it to a existing database.

Usage:
    gff_to_db.py -g <gff_file> -i <sqlite_file> 

Reference:
    http://pythonhosted.org/gffutils/contents.html

"""
import os
import gffutils
import argparse
import sys


def main( gffFile, inputDB, newDB ):
    # creating a new database or adding to an existing one
    if inputDB:
        addToDB( gffFile,inputDB )
    else:
        createDB( gffFile,newDB )


def addToDB( gffFile, inputDB ):
    # opening the database
    try:
        db = gffutils.FeatureDB(dbfn=inputDB)
    except:
        print("can't open input database")
    try:
        db.update(gffFile, make_backup=False, id_spec=':source:', merge_strategy='create_unique')
    except:
        print("can't add files to database; source key is already in use; please change the source description with the \'gff_source_editing\' tool")


def createDB(gffFile, newDB):
    # creating a new database
    gffutils.create_db(gffFile, dbfn=newDB, id_spec=':source:', merge_strategy='create_unique')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'gff_to_db.py', description='Create a database from a gff file or add to an excisting database.', prefix_chars='-+', epilog="")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--gff', '-g', dest='gffFile', required=True, help='input GFF file')
    parser.add_argument('--input_database', '-i', dest='inputDB', help='write to existing database')
    parser.add_argument('--new_database', '-n', dest='newDB', help='create new database')

    options = parser.parse_args()
    if not (options.inputDB or options.newDB):
        parser.error("One of --input_database or --new_database must be given")

    main(options.gffFile, options.inputDB, options.newDB)

