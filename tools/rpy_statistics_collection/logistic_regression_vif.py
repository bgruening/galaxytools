#!/usr/bin/env python

# from galaxy import eggs
import sys

import rpy2.rinterface as ri
import rpy2.rlike.container as rlc
# from rpy import *
import rpy2.robjects as robjects

r = robjects.r


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


infile = sys.argv[1]
y_col = int(sys.argv[2]) - 1
x_cols = sys.argv[3].split(",")
outfile = sys.argv[4]


print("Predictor columns: %s; Response column: %d" % (x_cols, y_col + 1))
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

y_vals = []
x_vals = []
x_vector = []
for k, col in enumerate(x_cols):
    x_cols[k] = int(col) - 1
    x_vals.append([])

NA = "NA"
for ind, line in enumerate(file(infile)):  # noqa F821
    if line and not line.startswith("#"):
        try:
            fields = line.split("\t")
            try:
                yval = float(fields[y_col])
            except Exception:
                yval = r("NA")
            y_vals.append(yval)
            for k, col in enumerate(x_cols):
                try:
                    xval = float(fields[col])
                except Exception:
                    xval = r("NA")
                x_vals[k].append(xval)
                x_vector.append(xval)
        except Exception as e:
            print(e)

# x_vals1 = numpy.asarray(x_vals).transpose()

check1 = 0
check0 = 0
for i in y_vals:
    if i == 1:
        check1 = 1
    if i == 0:
        check0 = 1
if check1 == 0 or check0 == 0:
    sys.exit("Warning: logistic regression must have at least two classes")

for i in y_vals:
    if i not in [1, 0, r("NA")]:
        print(str(i), file=fout)
        sys.exit(
            "Warning: the current version of this tool can run only with two classes and need to be labeled as 0 and 1."
        )


# dat= r.list(x=array(x_vals1), y=y_vals)
novif = 0
# set_default_mode(NO_CONVERSION)
# try:
#    linear_model = r.glm(r("y ~ x"), data = r.na_exclude(dat),family="binomial")
#    #r('library(car)')
#    #r.assign('dat',dat)
#    #r.assign('ncols',len(x_cols))
#    #r.vif(r('glm(dat$y ~ ., data = na.exclude(data.frame(as.matrix(dat$x,ncol=ncols))->datx),family="binomial")')).as_py()
#
# except Exception as rex:
#    stop_err("Error performing logistic regression on the input data.\nEither the response column or one of the predictor columns contain only non-numeric or invalid values.")

fv = robjects.FloatVector(x_vector)
m = r["matrix"](fv, ncol=len(x_cols), byrow=True)
# ensure order for generating formula
od = rlc.OrdDict([("y", robjects.FloatVector(y_vals)), ("x", m)])
dat = robjects.DataFrame(od)
# convert dat.names: ["y","x.1","x.2"] to formula string: 'y ~ x.1 + x.2'
formula = " + ".join(dat.names).replace("+", "~", 1)
print(formula)
try:
    linear_model = r.glm(formula, data=r["na.exclude"](dat), family="binomial")
except Exception:
    stop_err(
        "Error performing linear regression on the input data.\nEither the response column or one of the predictor columns contain only non-numeric or invalid values."
    )

if len(x_cols) > 1:
    try:
        r("library(car)")
        r.assign("dat", dat)
        r.assign("ncols", len(x_cols))
        # vif=r.vif(r('glm(dat$y ~ ., data = na.exclude(data.frame(as.matrix(dat$x,ncol=ncols))->datx),family="binomial")'))
        od2 = rlc.OrdDict([("datx", m)])
        glm_data_frame = robjects.DataFrame(od2)
        glm_result = r.glm(
            "dat$y ~ .", data=r["na.exclude"](glm_data_frame), family="binomial"
        )
        print("Have glm")
        vif = r.vif(glm_result)
    except Exception as rex:
        print(rex)
else:
    novif = 1

# set_default_mode(BASIC_CONVERSION)

# coeffs=linear_model.as_py()['coefficients']
coeffs = linear_model.rx2("coefficients")
# null_deviance=linear_model.as_py()['null.deviance']
null_deviance = linear_model.rx2("null.deviance")[0]
# residual_deviance=linear_model.as_py()['deviance']
residual_deviance = linear_model.rx2("deviance")[0]
# yintercept= coeffs['(Intercept)']
yintercept = coeffs.rx2("(Intercept)")[0]

summary = r.summary(linear_model)
# co = summary.get('coefficients', 'NA')
co = summary.rx2("coefficients")
print(co)
"""
if len(co) != len(x_vals)+1:
    stop_err("Stopped performing logistic regression on the input data, since one of the predictor columns contains only non-numeric or invalid values.")
"""

try:
    yintercept = r.round(float(yintercept), digits=10)[0]
    # pvaly = r.round(float(co[0][3]), digits=10)
    pvaly = r.round(float(co.rx(1, 4)[0]), digits=10)[0]
except Exception as e:
    print(str(e))
print("response column\tc%d" % (y_col + 1), file=fout)
tempP = []
for i in x_cols:
    tempP.append("c" + str(i + 1))
tempP = ",".join(tempP)
print("predictor column(s)\t%s" % (tempP), file=fout)
print("Y-intercept\t%s" % (yintercept), file=fout)
print("p-value (Y-intercept)\t%s" % (pvaly), file=fout)

print(coeffs)
if len(x_vals) == 1:  # Simple linear  regression case with 1 predictor variable
    try:
        # slope = r.round(float(coeffs['x']), digits=10)
        raw_slope = coeffs.rx2("x")[0]
        slope = r.round(float(raw_slope), digits=10)[0]
    except Exception:
        slope = "NA"
    try:
        # pval = r.round(float(co[1][3]), digits=10)
        pval = r.round(float(co.rx2(2, 4)[0]), digits=10)[0]
    except Exception:
        pval = "NA"
    print("Slope (c%d)\t%s" % (x_cols[0] + 1, slope), file=fout)
    print("p-value (c%d)\t%s" % (x_cols[0] + 1, pval), file=fout)
else:  # Multiple regression case with >1 predictors
    ind = 1
    # while ind < len(coeffs.keys()):
    print(len(coeffs.names))
    while ind < len(coeffs.names):
        try:
            # slope = r.round(float(coeffs['x'+str(ind)]), digits=10)
            raw_slope = coeffs.rx2("x." + str(ind))[0]
            slope = r.round(float(raw_slope), digits=10)[0]
        except Exception:
            slope = "NA"
        print("Slope (c%d)\t%s" % (x_cols[ind - 1] + 1, slope), file=fout)

        try:
            # pval = r.round(float(co[ind][3]), digits=10)
            pval = r.round(float(co.rx2(ind + 1, 4)[0]), digits=10)[0]
        except Exception:
            pval = "NA"
        print("p-value (c%d)\t%s" % (x_cols[ind - 1] + 1, pval), file=fout)
        ind += 1

# rsq = summary.get('r.squared','NA')
rsq = summary.rx2("r.squared")
if rsq == ri.RNULLType():
    rsq = "NA"
else:
    rsq = rsq[0]


try:
    # rsq= r.round(float((null_deviance-residual_deviance)/null_deviance), digits=5)
    rsq = r.round(float((null_deviance - residual_deviance) / null_deviance), digits=5)[
        0
    ]
    # null_deviance= r.round(float(null_deviance), digits=5)
    null_deviance = r.round(float(null_deviance), digits=5)[0]
    # residual_deviance= r.round(float(residual_deviance), digits=5)
    residual_deviance = r.round(float(residual_deviance), digits=5)[0]

except Exception:
    pass

print("Null deviance\t%s" % (null_deviance), file=fout)

print("Residual deviance\t%s" % (residual_deviance), file=fout)
print("pseudo R-squared\t%s" % (rsq), file=fout)
print("\n", file=fout)
print("vif", file=fout)

if novif == 0:
    # py_vif=vif.as_py()
    count = 0
    for i in sorted(vif.names):
        print("c" + str(x_cols[count] + 1), str(vif.rx2(i)[0]), file=fout)
        count += 1
elif novif == 1:
    print("vif can calculate only when model have more than 1 predictor", file=fout)
