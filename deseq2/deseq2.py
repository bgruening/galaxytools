#!/usr/bin/env python
import os, sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
deseq2 = importr("DESeq2")

sys.stdout.write('Works!')

