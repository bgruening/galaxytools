#!/usr/bin/env Rscript

library(epicseg)
epicseg:::CLI(args=commandArgs(trailingOnly=TRUE), epicseg:::getProg())
