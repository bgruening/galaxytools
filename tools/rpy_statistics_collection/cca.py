#!/usr/bin/env python

import string
import sys

import numpy
import rpy2.rlike.container as rlc

# from rpy import *
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

r = robjects.r
grdevices = importr("grDevices")


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


infile = sys.argv[1]
x_cols = sys.argv[2].split(",")
y_cols = sys.argv[3].split(",")

x_scale = x_center = "FALSE"
if sys.argv[4] == "both":
    x_scale = x_center = "TRUE"
elif sys.argv[4] == "center":
    x_center = "TRUE"
elif sys.argv[4] == "scale":
    x_scale = "TRUE"

y_scale = y_center = "FALSE"
if sys.argv[5] == "both":
    y_scale = y_center = "TRUE"
elif sys.argv[5] == "center":
    y_center = "TRUE"
elif sys.argv[5] == "scale":
    y_scale = "TRUE"

std_scores = "FALSE"
if sys.argv[6] == "yes":
    std_scores = "TRUE"

outfile = sys.argv[7]
outfile2 = sys.argv[8]

fout = open(outfile, "w")
elems = []
for i, line in enumerate(file(infile)):
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

x_vals = []

for k, col in enumerate(x_cols):
    x_cols[k] = int(col) - 1
    # x_vals.append([])

y_vals = []

for k, col in enumerate(y_cols):
    y_cols[k] = int(col) - 1
    # y_vals.append([])

skipped = 0
for ind, line in enumerate(file(infile)):
    if line and not line.startswith("#"):
        try:
            fields = line.strip().split("\t")
            valid_line = True
            for col in x_cols + y_cols:
                try:
                    assert float(fields[col])
                except Exception:
                    skipped += 1
                    valid_line = False
                    break
            if valid_line:
                for k, col in enumerate(x_cols):
                    try:
                        xval = float(fields[col])
                    except Exception:
                        xval = NaN  #
                    # x_vals[k].append(xval)
                    x_vals.append(xval)
                for k, col in enumerate(y_cols):
                    try:
                        yval = float(fields[col])
                    except Exception:
                        yval = NaN  #
                    # y_vals[k].append(yval)
                    y_vals.append(yval)
        except Exception:
            skipped += 1

# x_vals1 = numpy.asarray(x_vals).transpose()
# y_vals1 = numpy.asarray(y_vals).transpose()

# x_dat= r.list(array(x_vals1))
# y_dat= r.list(array(y_vals1))

x_dat = r["matrix"](robjects.FloatVector(x_vals), ncol=len(x_cols), byrow=True)
y_dat = r["matrix"](robjects.FloatVector(y_vals), ncol=len(y_cols), byrow=True)

try:
    r.suppressWarnings(r.library("yacca"))
except Exception:
    stop_err("Missing R library yacca.")

# set_default_mode(NO_CONVERSION)
try:
    xcolnames = ["c%d" % (el + 1) for el in x_cols]
    ycolnames = ["c%d" % (el + 1) for el in y_cols]
    # cc = r.cca(x=x_dat, y=y_dat, xlab=xcolnames, ylab=ycolnames, xcenter=r(x_center), ycenter=r(y_center), xscale=r(x_scale), yscale=r(y_scale), standardize_scores=r(std_scores))
    cc = r.cca(
        x=x_dat,
        y=y_dat,
        xlab=xcolnames,
        ylab=ycolnames,
        xcenter=r(x_center),
        ycenter=r(y_center),
        xscale=r(x_scale),
        yscale=r(y_scale),
        **{"standardize.scores": r(std_scores)}
    )
    # ftest = r.F_test_cca(cc)
    ftest = r["F.test.cca"](cc)
except RException, rex:
    stop_err("Encountered error while performing CCA on the input data: %s" % (rex))

# set_default_mode(BASIC_CONVERSION)
summary = r.summary(cc)

# ncomps = len(summary['corr'])
ncomps = len(summary.rx2("corr"))
# comps = summary['corr'].keys()
# comps = summary.rx2('corr').names
comps = (",".join(summary.rx2("corr").names)).split(",")
# corr = summary['corr'].values()
corr = summary.rx2("corr")
# xlab = summary['xlab']
xlab = summary.rx2("xlab")
# ylab = summary['ylab']
ylab = summary.rx2("ylab")

for i in range(ncomps):
    corr[comps.index("CV %s" % (i + 1))] = summary.rx2("corr")[i]
    # corr[comps.index('CV %s' %(i+1))] = summary['corr'].values()[i]

# ftest=ftest.as_py()
print("#Component\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
print("#Correlation\t%s" % ("\t".join(["%.4g" % el for el in corr])), file=fout)
# print >>fout, "#F-statistic\t%s" %("\t".join(["%.4g" % el for el in ftest['statistic']]))
print("#F-statistic\t%s" % (
    "\t".join(["%.4g" % el for el in ftest.rx2("statistic")])
), file=fout)
# print >>fout, "#p-value\t%s" %("\t".join(["%.4g" % el for el in ftest['p.value']]))
print("#p-value\t%s" % (
    "\t".join(["%.4g" % el for el in ftest.rx2("p.value")])
), file=fout)


print("#X-Coefficients\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# for i,val in enumerate(summary['xcoef']):
#    print >>fout, "%s\t%s" %(xlab[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2("xcoef")
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (xlab[i][0], "\t".join(vals)), file=fout)

print("#Y-Coefficients\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# for i,val in enumerate(summary['ycoef']):
#    print >>fout, "%s\t%s" %(ylab[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2("ycoef")
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (ylab[i][0], "\t".join(vals)), file=fout)

print("#X-Loadings\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# for i,val in enumerate(summary['xstructcorr']):
#    print >>fout, "%s\t%s" %(xlab[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2("xstructcorr")
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (xlab[i][0], "\t".join(vals)), file=fout)

print("#Y-Loadings\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# for i,val in enumerate(summary['ystructcorr']):
#    print >>fout, "%s\t%s" %(ylab[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2("ystructcorr")
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (ylab[i][0], "\t".join(vals)), file=fout)

print("#X-CrossLoadings\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# for i,val in enumerate(summary['xcrosscorr']):
#    print >>fout, "%s\t%s" %(xlab[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2("xcrosscorr")
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (xlab[i][0], "\t".join(vals)), file=fout)

print("#Y-CrossLoadings\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# for i,val in enumerate(summary['ycrosscorr']):
#    print >>fout, "%s\t%s" %(ylab[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2("ycrosscorr")
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (ylab[i][0], "\t".join(vals)), file=fout)

r.pdf(outfile2, 8, 8)
# r.plot(cc)
for i in range(ncomps):
    r["helio.plot"](
        cc, cv=i + 1, main=r.paste("Explained Variance for CV", i + 1), type="variance"
    )
# r.dev_off()
grdevices.dev_off()
