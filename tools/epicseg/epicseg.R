#!/usr/lib/R/bin/Rscript
library(epicseg)
epicseg:::CLI(args=commandArgs(trailingOnly=TRUE), epicseg:::getProg())
