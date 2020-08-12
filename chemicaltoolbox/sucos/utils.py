#!/usr/bin/env python
"""
Utility functions for SuCOS and other RDKit modules
"""

from __future__ import print_function
import sys, gzip
from rdkit import Chem

def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)

def open_file_for_reading(filename):
    """Open the file gunzipping it if it ends with .gz."""
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'rb')
    else:
        return open(filename, 'rb')

def open_file_for_writing(filename):
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'at')
    else:
        return open(filename, 'w+')

def read_single_molecule(filename, index=1, format=None):
    """Read a single molecule as a RDKit Mol object. This can come from a file in molfile or SDF format.
    If SDF then you can also specify an index of the molecule that is read (default is the first)
    """
    mol = None
    if format == 'mol' or filename.lower().endswith('.mol') or filename.lower().endswith('.mol.gz'):
        file = open_file_for_reading(filename)
        mol = Chem.MolFromMolBlock(file.read())
        file.close()
    elif format == 'sdf' or filename.lower().endswith('.sdf') or filename.lower().endswith('.sdf.gz'):
        file = open_file_for_reading(filename)
        supplier = Chem.ForwardSDMolSupplier(file)
        for i in range(0,index):
            if supplier.atEnd():
                break
            mol = next(supplier)
        file.close()

    if not mol:
        raise ValueError("Unable to read molecule")

    return mol