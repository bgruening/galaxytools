#!/usr/bin/env python

import sys

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
method = sys.argv[3]
outfile = sys.argv[4]
outfile2 = sys.argv[5]

if method == "svd":
    scale = center = "FALSE"
    if sys.argv[6] == "both":
        scale = center = "TRUE"
    elif sys.argv[6] == "center":
        center = "TRUE"
    elif sys.argv[6] == "scale":
        scale = "TRUE"

fout = open(outfile, "w")
elems = []
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

x_vals = []

for k, col in enumerate(x_cols):
    x_cols[k] = int(col) - 1
    # x_vals.append([])

NA = "NA"
skipped = 0
for ind, line in enumerate(file(infile)):  # noqa F821
    if line and not line.startswith("#"):
        try:
            fields = line.strip().split("\t")
            valid_line = True
            for k, col in enumerate(x_cols):
                try:
                    xval = float(fields[col])
                except Exception:
                    skipped += 1
                    valid_line = False
                    break
            if valid_line:
                for k, col in enumerate(x_cols):
                    xval = float(fields[col])
                    # x_vals[k].append(xval)
                    x_vals.append(xval)
        except Exception:
            skipped += 1

# x_vals1 = numpy.asarray(x_vals).transpose()
# dat= r.list(array(x_vals1))
dat = r["matrix"](robjects.FloatVector(x_vals), ncol=len(x_cols), byrow=True)

# set_default_mode(NO_CONVERSION)
try:
    if method == "cor":
        # pc = r.princomp(r.na_exclude(dat), cor = r("TRUE"))
        pc = r.princomp(r["na.exclude"](dat), cor=r("TRUE"))
    elif method == "cov":
        # pc = r.princomp(r.na_exclude(dat), cor = r("FALSE"))
        pc = r.princomp(r["na.exclude"](dat), cor=r("FALSE"))
    elif method == "svd":
        # pc = r.prcomp(r.na_exclude(dat), center = r(center), scale = r(scale))
        pc = r.prcomp(r["na.exclude"](dat), center=r(center), scale=r(scale))
# except Exception as rex:
except Exception as rex:  # need to find rpy2 RException
    stop_err("Encountered error while performing PCA on the input data: %s" % (rex))

# set_default_mode(BASIC_CONVERSION)
summary = r.summary(pc, loadings="TRUE")
# ncomps = len(summary['sdev'])
ncomps = len(summary.rx2("sdev"))

# if type(summary['sdev']) == type({}):
#    comps_unsorted = summary['sdev'].keys()
#    comps=[]
#    sd = summary['sdev'].values()
#    for i in range(ncomps):
#        sd[i] = summary['sdev'].values()[comps_unsorted.index('Comp.%s' %(i+1))]
#        comps.append('Comp.%s' %(i+1))
# elif type(summary['sdev']) == type([]):
#    comps=[]
#    for i in range(ncomps):
#        comps.append('Comp.%s' %(i+1))
#        sd = summary['sdev']

comps = []
for i in range(ncomps):
    comps.append("Comp.%s" % (i + 1))
sd = summary.rx2("sdev")

print("#Component\t%s" % (
    "\t".join(["%s" % el for el in range(1, ncomps + 1)])
), file=fout)
# print >>fout, "#Std. deviation\t%s" %("\t".join(["%.4g" % el for el in sd]))
print("#Std. deviation\t%s" % ("\t".join(["%.4g" % el for el in sd])), file=fout)
total_var = 0
vars = []
for s in sd:
    var = s * s
    total_var += var
    vars.append(var)
for i, var in enumerate(vars):
    vars[i] = vars[i] / total_var

print("#Proportion of variance explained\t%s" % (
    "\t".join(["%.4g" % el for el in vars])
), file=fout)

print("#Loadings\t%s" % ("\t".join(["%s" % el for el in range(1, ncomps + 1)])), file=fout)
xcolnames = ["c%d" % (el + 1) for el in x_cols]
# if 'loadings' in summary: #in case of princomp
if "loadings" in summary.names:  # in case of princomp
    loadings = "loadings"
# elif 'rotation' in summary: #in case of prcomp
elif "rotation" in summary.names:  # in case of prcomp
    loadings = "rotation"
# for i,val in enumerate(summary[loadings]):
#    print >>fout, "%s\t%s" %(xcolnames[i], "\t".join(["%.4g" % el for el in val]))
vm = summary.rx2(loadings)
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (xcolnames[i], "\t".join(vals)), file=fout)

print("#Scores\t%s" % ("\t".join(["%s" % el for el in range(1, ncomps + 1)])), file=fout)
# if 'scores' in summary: #in case of princomp
if "scores" in summary.names:  # in case of princomp
    scores = "scores"
# elif 'x' in summary: #in case of prcomp
elif "x" in summary.names:  # in case of prcomp
    scores = "x"
# for obs,sc in enumerate(summary[scores]):
#    print >>fout, "%s\t%s" %(obs+1, "\t".join(["%.4g" % el for el in sc]))
vm = summary.rx2(scores)
for i in range(vm.nrow):
    vals = []
    for j in range(vm.ncol):
        vals.append("%.4g" % vm.rx2(i + 1, j + 1)[0])
    print("%s\t%s" % (i + 1, "\t".join(vals)), file=fout)
r.pdf(outfile2, 8, 8)
r.biplot(pc)
# r.dev_off()
grdevices.dev_off()
