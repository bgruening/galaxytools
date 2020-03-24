#
# Copyright (C) 2015 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import argparse, sys
from builtins import range

from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from numpy import linalg

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


def write_out(mols,count,writer,file_format):
    for mol in mols:
        count += 1
        if mol is None: continue
        if file_format == 'sdf':
            writer.write(mol)
        elif file_format == 'json':
            writer.write(mol, format='mol')
    return count

def GetBestFitPlane(pts, weights=None):
  if weights is None:
    wSum = len(pts)
    origin = np.sum(pts, 0)
  origin /= wSum
  sums = np.zeros((3, 3), np.double)
  for pt in pts:
    dp = pt - origin
    for i in range(3):
      sums[i, i] += dp[i] * dp[i]
      for j in range(i + 1, 3):
        sums[i, j] += dp[i] * dp[j]
        sums[j, i] += dp[i] * dp[j]
  sums /= wSum
  vals, vects = linalg.eigh(sums)
  order = np.argsort(vals)
  normal = vects[:, order[0]]
  plane = np.zeros((4, ), np.double)
  plane[:3] = normal
  plane[3] = -1 * normal.dot(origin)
  return plane

def PBFRD(mol, confId=-1):
  conf = mol.GetConformer(confId)
  if not conf.Is3D():
    return 0

  pts = np.array([list(conf.GetAtomPosition(x)) for x in range(mol.GetNumAtoms())])
  plane = GetBestFitPlane(pts)
  denom = np.dot(plane[:3], plane[:3])
  denom = denom**0.5
  # add up the distance from the plane for each point:
  res = 0.0
  for pt in pts:
    res += np.abs(pt.dot(plane[:3]) + plane[3])
  res /= denom
  res /= len(pts)
  return res

def PBFev(mol):
    '''returns an array of exit vectors for this mol'''
    # Get murcko SMILES
    murcko = MurckoScaffold.GetScaffoldForMol(mol)

    # Get PBF plane for murcko scaffold only
    confId = -1
    conf = murcko.GetConformer(confId)
    if not conf.Is3D():
        print('This mol is not 3D - all PBFev angles will be 0 degrees')
        return [0]
    pts = np.array([list(conf.GetAtomPosition(i))  # Get atom coordinates
                    for i in range(murcko.GetNumAtoms())])
    # GetBestFitPlane is in the RDKit Contrib directory as part of PBF
    # Plane is xyz vector with a c intercept adjustment
    plane = GetBestFitPlane(pts)

    # Map onto parent structure coords (this func adds exit vectors [*])
    murckoEv = Chem.ReplaceSidechains(mol, murcko)

    confId = -1  # embed 3D conf object with EVs (atom indices do not change)
    conf = murckoEv.GetConformer(confId)

    # Where [#0] matches exit vector SMILES [*]
    patt = Chem.MolFromSmarts('[#0]-[*]')
    matches = murckoEv.GetSubstructMatches(patt)
    if len(matches) == 0:
        return None

    # Calculate angles between exit vectors and the murcko plane of best fit
    exitVectors = np.zeros(len(matches))
    denom = np.dot(plane[:3], plane[:3])
    denom = denom**0.5
    for n, match in enumerate(matches):
        evCoords = conf.GetAtomPosition(match[0])
        anchorCoords = conf.GetAtomPosition(match[1])
        v = np.array(((evCoords[0]-anchorCoords[0]),
                      (evCoords[1]-anchorCoords[1]),
                      (evCoords[2]-anchorCoords[2])))
        angle = np.arcsin((np.dot(v, plane[:3])) /
                          ((denom)*((np.dot(v, v))**0.5)))
        angle = np.abs(np.degrees(angle))
        exitVectors[n] = angle
    return exitVectors

def main():

    ### command line args defintions #########################################
    parser = argparse.ArgumentParser(description='Calculate plane of best fit for molecules')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    parameter_utils.add_default_io_args(parser)
    args = parser.parse_args()
    args.logfile.write("PBFEV args: %s\n\n" % str(args))
    input ,output ,suppl ,writer ,output_base = rdkit_utils.default_open_input_output(args.input, args.informat, args.output, 'PBFEV', args.outformat)

    i=0
    count=0
    errors=0
    writer = rdkit_utils.ThickSDWriter(args.output)

    out_results = []
    for mol in suppl:
        i +=1
        AllChem.EmbedMolecule(mol)
        if mol is None: continue
        out_vector = PBFev(mol)
        if out_vector is None: continue
        rd = PBFRD(mol)
        mol.SetDoubleProp("distance",rd)
        for j,angle in enumerate(out_vector):
            mol.SetDoubleProp("angle"+"_"+str(j), angle)
        out_results.append(mol)
    count = write_out(out_results,count,writer,args.outformat)
    args.logfile.write("Handled %d molecules, resulting in %d outputs" % (i, count))
    writer.flush()
    writer.close()
    input.close()
    output.close()

if __name__ == "__main__":
    main()