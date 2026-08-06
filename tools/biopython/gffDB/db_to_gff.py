#!/usr/bin/env python

"""
Description:
    Get information from a gffutils database and create a new gff file.

Usage:
    db_to_gff.py -g <gff_file> -d <database_file>

Reference:
    http://pythonhosted.org/gffutils/contents.html

"""
import sys
import gffutils
import argparse
import sqlite3


## main function
def main(gffFile, dataBase, batchSearch, all, features, source, seqid, start, end, strand, within):
    # opening the database
    try:
        db = gffutils.FeatureDB(dbfn=dataBase)
    except sqlite3.DatabaseError:
        print("Error: Can't open database")
        sys.exit()

    # getting all data
    if all:
        data = db.all_features()
        printData(data, gffFile, source, 'w')
    # getting data from batch file
    elif batchSearch:
        batchQuery(db, gffFile, batchSearch, source)
    # getting single request
    else:
        singleQuery(db, gffFile, features, source, seqid, start, end, strand, within)


## check items
def batchQuery(db, gffFile, batchSearch, source):
    with open(batchSearch,'r') as batch:
        first = True
        lineNum = 1
        for line in batch:
            lineList = line.split('\t')
            # check if the batch file is consistent
            if len(lineList) < 6:
                print("Error: Batch file line " + str(lineNum) + " has less than 6 entries")
                sys.exit()
            seqid, features, start, end, strand, within = lineList
            # remove new line character
            within = within.strip()
            # check if the items are consistent
            b_seqid, b_feature, b_start, b_end, b_strand, b_within = checkItems(seqid, features, start, end, strand, within)
            # use best search strategy for the requested task
            if b_seqid or b_start or b_end:
                data = db.region(seqid=b_seqid, start=b_start, end=b_end, strand=b_strand, featuretype=b_feature, completely_within=b_within)
            elif b_feature:
                data = db.features_of_type(b_feature, strand=b_strand)
            else:
                data = db.all_features(strand=b_strand)
            # print data
            if first:
                printData(data, gffFile, source, 'w')
                first = False
            else:
                printData(data, gffFile, source, 'a')
            lineNum += 1


## check items
def singleQuery(db, gffFile, features, source, seqid, start, end, strand, within):
    # check if the items are consistent
    c_seqid, c_feature, c_start, c_end, c_strand, c_within = checkItems(seqid, features, start, end, strand, within)
    featureList = None
    if c_feature:
        featureList = c_feature.split(',')
        if len(featureList) > 1:
            c_feature = featureList
    # use best search strategy for the requested task
    if c_seqid or c_start or c_end:
        data = db.region(seqid=c_seqid, start=c_start, end=c_end, strand=c_strand, featuretype=c_feature, completely_within=c_within)
    elif c_feature:
        data = db.features_of_type(c_feature, strand=c_strand)
    else:
        data = db.all_features(strand=c_strand)
    printData(data, gffFile, source, 'w')


## check items
def checkItems(seqid, feature, start, end, strand, within):
    # check seqid
    if not seqid or str(seqid) == '.':
        seqid = None
    # check features
    if not feature or str(feature) == '.':
        feature = None
    # check start
    if not start or str(start) == '.':
        start = None
    else:
        try:
            n = int(start)
        except ValueError:
            print("Error: Start \'" + start + "\' is not convertible to integer")
            sys.exit()
    # check end
    if not end or str(end) == '.':
        end = None
    else:
        try:
            n = int(end)
        except ValueError:
            print("Error: End \'" + end + "\' is not convertible to integer")
            sys.exit()
    # check strand
    if not strand or str(strand) == '.':
        strand = None
    elif not (str(strand) == '+' or str(strand) == '-' or str(strand) == 'x'):
        print("Error: Strand \'" + strand + "\' is not a valid strand item")
        sys.exit()
    # check within
    if within == True or str(within) == 'w':
        within = True
    else:
        within = False
    return seqid, feature, start, end, strand, within


def printData(data, gffFile, source, mode):

    with open( gffFile, mode ) as handle:
        for entry in data:
            if source and entry.source != source:
                try:
                    next(data)
                except StopIteration:
                    pass
            else:
                handle.write(str(entry) + "\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'db_to_gff.py',
        description = 'Get information from a gffutils database and create a new gff file.', prefix_chars='-+', epilog="")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--gff', '-g', dest='gffFile', required=True, help='GFF file name')
    parser.add_argument('--database', '-d', dest='dataBase', required=True, help='gffutils database')
    parser.add_argument('--batch', '-b', dest='batchSearch',
        help='tabular file with search data; format: seqid - feature - start - end - strand - within; values can be empty or use \'.\' as a place holder; always provide 5 tabs; use \'x\' for unstranded features; use \'w\' for within')
    parser.add_argument('--all', '-a', dest='all', action='store_true', help='returns all db entries')
    parser.add_argument('--features', '-f', dest='features', help='returns all entries with the requested features; can be a comma separated list')
    parser.add_argument('--source', '-so', dest='source', help='returns features with the requested source')
    parser.add_argument('--seqid', '-sq', dest='seqid', help='returns features with the requested seqid')
    parser.add_argument('--start', '-s', dest='start', help='only returns features that start with this region or after this region')
    parser.add_argument('--end', '-e', dest='end', help='only returns features that end with this region or before this region')
    parser.add_argument('--strand', '-r', dest='strand', choices=['+', '-', 'x'], help='returns only features in strand direction; \'x\' returns unstranded features')
    parser.add_argument('--within', '-w', dest='within', action='store_true',
        help='forces the feature to be completely within the provided --start and/or --end region')
    
    options = parser.parse_args()
    main(options.gffFile, options.dataBase, options.batchSearch, options.all, options.features, options.source, options.seqid, options.start, options.end, options.strand, options.within)

