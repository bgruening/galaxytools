#!/usr/bin/env python
""""""

__author__ = "Tim te Beek"
__copyright__ = "Copyright 2013, Netherlands Bioinformatics Centre"
__license__ = "MIT"

#First argument contains fully qualified name of module to be imported
import sys
assert len(sys.argv) == 4, "Usage: ./gbk2ffn.py genbank.gbk label fasta.ffn"
genbank = sys.argv[1]
label = sys.argv[2]
fasta = sys.argv[3]

import tempfile
tmpfolder = tempfile.mkdtemp()

from divergence.translate import _extract_gene_and_protein
dnafile, aafile = _extract_gene_and_protein(tmpfolder, label, genbank, filetype='genbank')

import shutil
shutil.move(dnafile, fasta)
