#!/usr/bin/env python

import sys, string
import rpy2.robjects as robjects
import rpy2.rlike.container as rlc
from rpy2.robjects.packages import importr
r = robjects.r
grdevices = importr('grDevices')
#  from rpy import *
import numpy


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()

infile = sys.argv[1]
y_col = int(sys.argv[2])-1
x_cols = sys.argv[3].split(',')
outfile = sys.argv[4]
outfile2 = sys.argv[5]

print "Predictor columns: %s; Response column: %d" %(x_cols,y_col+1)
fout = open(outfile,'w')
elems = []
for i, line in enumerate( file ( infile )):
    line = line.rstrip('\r\n')
    if len( line )>0 and not line.startswith( '#' ):
        elems = line.split( '\t' )
        break 
    if i == 30:
        break # Hopefully we'll never get here...

if len( elems )<1:
    stop_err( "The data in your input dataset is either missing or not formatted properly." )

y_vals = []
x_vals = []

for k,col in enumerate(x_cols):
    x_cols[k] = int(col)-1
    # x_vals.append([])

NA = 'NA'
for ind,line in enumerate( file( infile )):
    if line and not line.startswith( '#' ):
        try:
            fields = line.split("\t")
            try:
                yval = float(fields[y_col])
            except:
                yval = r('NA')
            y_vals.append(yval)
            for k,col in enumerate(x_cols):
                try:
                    xval = float(fields[col])
                except:
                    xval = r('NA')
                # x_vals[k].append(xval)
                x_vals.append(xval)
        except:
            pass
# x_vals1 = numpy.asarray(x_vals).transpose()
# dat= r.list(x=array(x_vals1), y=y_vals)
fv = robjects.FloatVector(x_vals)
m = r['matrix'](fv, ncol=len(x_cols),byrow=True)
# ensure order for generating formula
od = rlc.OrdDict([('y',robjects.FloatVector(y_vals)),('x',m)])
dat = robjects.DataFrame(od)
# convert dat.names: ["y","x.1","x.2"] to formula string: 'y ~ x.1 + x.2'
formula = ' + '.join(dat.names).replace('+','~',1)

#set_default_mode(NO_CONVERSION)
try:
    #linear_model = r.lm(r("y ~ x"), data = r.na_exclude(dat))
    linear_model = r.lm(formula,  data =  r['na.exclude'](dat))
except RException, rex:
    stop_err("Error performing linear regression on the input data.\nEither the response column or one of the predictor columns contain only non-numeric or invalid values.")
#set_default_mode(BASIC_CONVERSION)

#coeffs=linear_model.as_py()['coefficients']
#yintercept= coeffs['(Intercept)']
coeffs=linear_model.rx2('coefficients')
yintercept= coeffs.rx2('(Intercept)')[0]
summary = r.summary(linear_model)

#co = summary.get('coefficients', 'NA')
co = summary.rx2("coefficients")

"""
if len(co) != len(x_vals)+1:
    stop_err("Stopped performing linear regression on the input data, since one of the predictor columns contains only non-numeric or invalid values.")
"""
#print >>fout, "p-value (Y-intercept)\t%s" %(co[0][3])
print >>fout, "p-value (Y-intercept)\t%s" %(co.rx(1,4)[0])

if len(x_vals) == 1:    #Simple linear  regression case with 1 predictor variable
    try:
        #slope = coeffs['x']
        slope = r.round(float(coeffs.rx2('x')[0]), digits=10)
    except:
        slope = 'NA'
    try:
        #pval = co[1][3]
        pval = r.round(float(co.rx(2,4)[0]), digits=10)
    except:
        pval = 'NA'
    print >>fout, "Slope (c%d)\t%s" %(x_cols[0]+1,slope)
    print >>fout, "p-value (c%d)\t%s" %(x_cols[0]+1,pval)
else:    #Multiple regression case with >1 predictors
    ind=1
    #while ind < len(coeffs.keys()):
    while ind < len(coeffs.names):
        # print >>fout, "Slope (c%d)\t%s" %(x_cols[ind-1]+1,coeffs['x'+str(ind)])
        print >>fout, "Slope (c%d)\t%s" %(x_cols[ind-1]+1,coeffs.rx2(coeffs.names[ind])[0])
        try:
            #pval = co[ind][3]
            pval = r.round(float(co.rx(ind+1,4)[0]), digits=10)
        except:
            pval = 'NA'
        print >>fout, "p-value (c%d)\t%s" %(x_cols[ind-1]+1,pval)
        ind+=1

rsq = summary.rx2('r.squared')[0]
adjrsq = summary.rx2('adj.r.squared')[0]
fstat = summary.rx2('fstatistic').rx2('value')[0]
sigma = summary.rx2('sigma')[0]

try:
    rsq = r.round(float(rsq), digits=5)
    adjrsq = r.round(float(adjrsq), digits=5)
    fval = r.round(fstat['value'], digits=5)
    fstat['value'] = str(fval)
    sigma = r.round(float(sigma), digits=10)
except:
    pass

print >>fout, "R-squared\t%s" %(rsq)
print >>fout, "Adjusted R-squared\t%s" %(adjrsq)
print >>fout, "F-statistic\t%s" %(fstat)
print >>fout, "Sigma\t%s" %(sigma)

r.pdf( outfile2, 8, 8 )
if len(x_vals) == 1:    #Simple linear  regression case with 1 predictor variable
    sub_title =  "Slope = %s; Y-int = %s" %(slope,yintercept)
    try:
        r.plot(x=x_vals[0], y=y_vals, xlab="X", ylab="Y", sub=sub_title, main="Scatterplot with regression")
        r.abline(a=yintercept, b=slope, col="red")
    except:
        pass
else:
    r.pairs(dat, main="Scatterplot Matrix", col="blue")
try:
    r.plot(linear_model)
except:
    pass
#r.dev_off()
grdevices.dev_off()
