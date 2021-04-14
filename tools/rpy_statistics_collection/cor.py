#!/usr/bin/env python
# Greg Von Kuster
"""
Calculate correlations between numeric columns in a tab delim file.
usage: %prog infile output.txt columns method
"""

import sys

# from rpy import *
import rpy2.robjects as robjects

r = robjects.r


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    method = sys.argv[4]
    assert method in ("pearson", "kendall", "spearman")

    try:
        columns = map(int, sys.argv[3].split(","))
    except:
        stop_err(
            "Problem determining columns, perhaps your query does not contain a column of numerical data."
        )

    matrix = []
    skipped_lines = 0
    first_invalid_line = 0
    invalid_value = ""
    invalid_column = 0

    for i, line in enumerate(file(sys.argv[1])):
        valid = True
        line = line.rstrip("\n\r")

        if line and not line.startswith("#"):
            # Extract values and convert to floats
            row = []
            for column in columns:
                column -= 1
                fields = line.split("\t")
                if len(fields) <= column:
                    valid = False
                else:
                    val = fields[column]
                    if val.lower() == "na":
                        row.append(float("nan"))
                    else:
                        try:
                            row.append(float(fields[column]))
                        except:
                            valid = False
                            skipped_lines += 1
                            if not first_invalid_line:
                                first_invalid_line = i + 1
                                invalid_value = fields[column]
                                invalid_column = column + 1
        else:
            valid = False
            skipped_lines += 1
            if not first_invalid_line:
                first_invalid_line = i + 1

        if valid:
            # matrix.append( row )
            matrix += row

    if skipped_lines < i:
        try:
            out = open(sys.argv[2], "w")
        except:
            stop_err("Unable to open output file")

        # Run correlation
        # print >> sys.stderr, "matrix: %s" % matrix
        # print >> sys.stderr, "array: %s" % array( matrix )
        try:
            # value = r.cor( array( matrix ), use="pairwise.complete.obs", method=method )
            fv = robjects.FloatVector(matrix)
            m = r["matrix"](fv, ncol=len(columns), byrow=True)
            rslt_mat = r.cor(m, use="pairwise.complete.obs", method=method)
            value = []
            for ri in range(1, rslt_mat.nrow + 1):
                row = []
                for ci in range(1, rslt_mat.ncol + 1):
                    row.append(rslt_mat.rx(ri, ci)[0])
                value.append(row)
        except Exception, exc:
            out.close()
            stop_err("%s" % str(exc))
        for row in value:
            print >> out, "\t".join(map(str, row))
        out.close()

    if skipped_lines > 0:
        msg = "..Skipped %d lines starting with line #%d. " % (
            skipped_lines,
            first_invalid_line,
        )
        if invalid_value and invalid_column > 0:
            msg += "Value '%s' in column %d is not numeric." % (
                invalid_value,
                invalid_column,
            )
        print msg


if __name__ == "__main__":
    main()
