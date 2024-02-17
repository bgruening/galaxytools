#!/usr/bin/env python

# from galaxy import eggs

import sys

import rpy2.rlike.container as rlc
import rpy2.robjects as robjects

# from rpy import *


r = robjects.r

# export PYTHONPATH=~/galaxy/lib/
# running command python partialR_square.py reg_inp.tab 4 1,2,3 partialR_result.tabular


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def sscombs(s):
    if len(s) == 1:
        return [s]
    else:
        ssc = sscombs(s[1:])
        return [s[0]] + [s[0] + comb for comb in ssc] + ssc


infile = sys.argv[1]
y_col = int(sys.argv[2]) - 1
x_cols = sys.argv[3].split(",")
outfile = sys.argv[4]

print("Predictor columns: %s; Response column: %d" % (x_cols, y_col + 1))
fout = open(outfile, "w")

for i, line in enumerate(file(infile)):  # noqa F821
    line = line.rstrip("\r\n")
    if len(line) > 0 and not line.startswith("#"):
        elems = line.split("\t")
        break
    if i == 30:
        break  # Hopefully we'll never get here...

if len(elems) < 1:
    stop_err(
        "The data in your input dataset is either missing or not formatted properly."
    )

y_vals = []
x_vals = []
x_vector = []
for k, col in enumerate(x_cols):
    x_cols[k] = int(col) - 1
    x_vals.append([])
    """
    try:
        float( elems[x_cols[k]] )
    except Exception:
        try:
            msg = "This operation cannot be performed on non-numeric column %d containing value '%s'." %( col, elems[x_cols[k]] )
        except Exception:
            msg = "This operation cannot be performed on non-numeric data."
        stop_err( msg )
    """
NA = "NA"
for ind, line in enumerate(file(infile)):  # noqa F821
    if line and not line.startswith("#"):
        try:
            fields = line.split("\t")
            try:
                yval = float(fields[y_col])
            except Exception:
                yval = r("NA")
                # print >>sys.stderr, "ey = %s" %ey
            y_vals.append(yval)
            for k, col in enumerate(x_cols):
                try:
                    xval = float(fields[col])
                except Exception:
                    xval = r("NA")
                    # print >>sys.stderr, "ex = %s" %ex
                x_vals[k].append(xval)
                x_vector.append(xval)
        except Exception:
            pass

# x_vals1 = numpy.asarray(x_vals).transpose()
# dat= r.list(x=array(x_vals1), y=y_vals)

# set_default_mode(NO_CONVERSION)
# try:
#    full = r.lm(r("y ~ x"), data= r.na_exclude(dat))    #full model includes all the predictor variables specified by the user
# except Exception as rex:
#    stop_err("Error performing linear regression on the input data.\nEither the response column or one of the predictor columns contain no numeric values.")
# set_default_mode(BASIC_CONVERSION)

fv = robjects.FloatVector(x_vector)
m = r["matrix"](fv, ncol=len(x_cols), byrow=True)
# ensure order for generating formula
od = rlc.OrdDict([("y", robjects.FloatVector(y_vals)), ("x", m)])
dat = robjects.DataFrame(od)
# convert dat.names: ["y","x.1","x.2"] to formula string: 'y ~ x.1 + x.2'
formula = " + ".join(dat.names).replace("+", "~", 1)
try:
    full = r.lm(formula, data=r["na.exclude"](dat))
except Exception:
    stop_err(
        "Error performing linear regression on the input data.\nEither the response column or one of the predictor columns contain only non-numeric or invalid values."
    )


summary = r.summary(full)
# fullr2 = summary.get('r.squared','NA')
fullr2 = summary.rx2("r.squared")[0]

if fullr2 == "NA":
    stop_err("Error in linear regression")

if len(x_vals) < 10:
    s = ""
    for ch in range(len(x_vals)):
        s += str(ch)
else:
    stop_err("This tool only works with less than 10 predictors.")

print("#Model\tR-sq\tpartial_R_Terms\tpartial_R_Value", file=fout)
all_combos = sorted(sscombs(s), key=len)
all_combos.reverse()
for j, cols in enumerate(all_combos):
    # if len(cols) == len(s):    #Same as the full model above
    #    continue
    if len(cols) == 1:
        # x_vals1 = x_vals[int(cols)]
        x_v = x_vals[int(cols)]
    else:
        x_v = []
        for col in cols:
            # x_v.append(x_vals[int(col)])
            x_v.extend(x_vals[int(col)])
        # x_vals1 = numpy.asarray(x_v).transpose()
    # dat= r.list(x=array(x_vals1), y=y_vals)
    # set_default_mode(NO_CONVERSION)
    # red = r.lm(r("y ~ x"), data= dat)    #Reduced model
    # set_default_mode(BASIC_CONVERSION)
    fv = robjects.FloatVector(x_v)
    m = r["matrix"](fv, ncol=len(cols), byrow=False)
    # ensure order for generating formula
    od = rlc.OrdDict([("y", robjects.FloatVector(y_vals)), ("x", m)])
    dat = robjects.DataFrame(od)
    # convert dat.names: ["y","x.1","x.2"] to formula string: 'y ~ x.1 + x.2'
    formula = " + ".join(dat.names).replace("+", "~", 1)
    try:
        red = r.lm(formula, data=r["na.exclude"](dat))
    except Exception:
        stop_err(
            "Error performing linear regression on the input data.\nEither the response column or one of the predictor columns contain only non-numeric or invalid values."
        )

    summary = r.summary(red)
    # redr2 = summary.get('r.squared','NA')
    redr2 = summary.rx2("r.squared")[0]

    try:
        partial_R = (float(fullr2) - float(redr2)) / (1 - float(redr2))
    except Exception:
        partial_R = "NA"
    col_str = ""
    for col in cols:
        col_str = col_str + str(int(x_cols[int(col)]) + 1) + " "
    col_str.strip()
    partial_R_col_str = ""
    for col in s:
        if col not in cols:
            partial_R_col_str = partial_R_col_str + str(int(x_cols[int(col)]) + 1) + " "
    partial_R_col_str.strip()
    if len(cols) == len(s):  # full model
        partial_R_col_str = "-"
        partial_R = "-"
    try:
        redr2 = "%.4f" % (float(redr2))
    except Exception:
        pass
    try:
        partial_R = "%.4f" % (float(partial_R))
    except Exception:
        pass
    print("%s\t%s\t%s\t%s" % (col_str, redr2, partial_R_col_str, partial_R), file=fout)
