#!/usr/bin/env python
"""
    Input: tabular format file with one column storing the unique id for the compounds and any other with the Spectrophores(TM) descriptors.
    Output: parse the target file using the same protocol used to generate the databases in our servers. Physico-chemical properties are computed and stored as metadata in the sdf output file.
    Copyright 2012, Bjoern Gruening and Xavier Lucas
"""
import sys, os
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel
import math
import numpy as np

#TODO get rid of eval()

global spectrophore
spectrophore = pybel.ob.OBSpectrophore()

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', required=True, help='target file name in sdf format with Spectrophores(TM) descriptors stored as meta-data')
    parser.add_argument('--library', required=True, help='library of compounds with pre-computed physico-chemical properties, including Spectrophores(TM) in tabular format')
    parser.add_argument('-c', '--column', required=True, type=int, help='#column containing the Spectrophores(TM) descriptors in the library file')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    parser.add_argument('-n', '--normalization', default="ZeroMeanAndUnitStd", choices=['No', 'ZeroMean', 'UnitStd', 'ZeroMeanAndUnitStd'], help='Normalization method')
    parser.add_argument('-a', '--accuracy', default="20", choices=['1', '2', '5', '10', '15', '20', '30', '36', '45', '60'], help='Accuracy expressed as angular stepsize')
    parser.add_argument('-s', '--stereo', default="No", choices=['No', 'Unique', 'Mirror', 'All'], help='Stereospecificity of the cage')
    parser.add_argument('-r', '--resolution', type=float, default="3.0", help='Resolution')
    return parser.parse_args()

def set_parameters(args):
    if args.normalization == 'No':
        spectrophore.SetNormalization( spectrophore.NoNormalization )
    else:
        spectrophore.SetNormalization( eval('spectrophore.NormalizationTowards' + args.normalization) )
    spectrophore.SetAccuracy( eval('spectrophore.AngStepSize' + args.accuracy) )
    spectrophore.SetStereo( eval('spectrophore.' + args.stereo + 'StereoSpecificProbes') )
    spectrophore.SetResolution( args.resolution )
    return True

def Compute_Spectrophores_distance(target_spectrophore, args):
    outfile = open(args.output, 'w')
    for mol in open(args.library, 'r'):
        try:
            distance = ( ( np.asarray( target_spectrophore, dtype=float ) - np.asarray( mol.split('\t')[ args.column - 1 ].strip().split(', '), dtype=float) )**2).sum()
        except ValueError:
            distance = 0
        outfile.write( '%s\t%f\n' % (mol.strip(), distance ) )
    outfile.close()

def __main__():
    """
        Computation of Spectrophores(TM) distances to a target molecule.
    """
    args = parse_command_line()
    # This sets up the parameters for the Spectrophore generation. Parameters are set to fit those of our standard parsing tool
    set_parameters(args)

    mol = pybel.readfile('sdf', args.target).next()
    target_spectrophore = mol.data["Spectrophores(TM)"].strip().split(', ')
    # Compute the paired-distance between every molecule in the library and the target
    distances = Compute_Spectrophores_distance(target_spectrophore, args)

if __name__ == "__main__" :
    __main__()
