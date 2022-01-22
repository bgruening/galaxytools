#!/usr/bin/python3

import gzip
import os
import sys
from optparse import OptionParser

from rdkit.Chem import AllChem as Chem

"""
This script was originally written by David Koes, University of Pittsburgh:
https://github.com/dkoes/rdkit-scripts/blob/master/rdconf.py
It is licensed under the MIT licence.

Given a smiles file, generate 3D conformers in output sdf.
Energy minimizes and filters conformers to meet energy window and rms constraints.

Some time ago I compared this to alternative conformer generators and
it was quite competitive (especially after RDKit's UFF implementation
added OOP terms).
"""


# convert smiles to sdf
def getRMS(mol, c1, c2):
    rms = Chem.GetBestRMS(mol, mol, c1, c2)
    return rms


parser = OptionParser(usage="Usage: %prog [options] <input>.smi <output>.sdf")
parser.add_option(
    "--maxconfs",
    dest="maxconfs",
    action="store",
    help="maximum number of conformers to generate per a molecule (default 20)",
    default="20",
    type="int",
    metavar="CNT",
)
parser.add_option(
    "--sample_multiplier",
    dest="sample",
    action="store",
    help="sample N*maxconfs conformers and choose the maxconformers with lowest energy (default 1)",
    default="1",
    type="float",
    metavar="N",
)
parser.add_option(
    "--seed",
    dest="seed",
    action="store",
    help="random seed (default 9162006)",
    default="9162006",
    type="int",
    metavar="s",
)
parser.add_option(
    "--rms_threshold",
    dest="rms",
    action="store",
    help="filter based on rms (default 0.7)",
    default="0.7",
    type="float",
    metavar="R",
)
parser.add_option(
    "--energy_window",
    dest="energy",
    action="store",
    help="filter based on energy difference with lowest energy conformer",
    default="10",
    type="float",
    metavar="E",
)
parser.add_option(
    "-v",
    "--verbose",
    dest="verbose",
    action="store_true",
    default=False,
    help="verbose output",
)
parser.add_option(
    "--mmff",
    dest="mmff",
    action="store_true",
    default=False,
    help="use MMFF forcefield instead of UFF",
)
parser.add_option(
    "--nomin",
    dest="nomin",
    action="store_true",
    default=False,
    help="don't perform energy minimization (bad idea)",
)
parser.add_option(
    "--etkdg",
    dest="etkdg",
    action="store_true",
    default=False,
    help="use new ETKDG knowledge-based method instead of distance geometry",
)


(options, args) = parser.parse_args()

if len(args) < 2:
    parser.error("Need input and output")
    sys.exit(-1)

input = args[0]
output = args[1]
smifile = open(input)
if options.verbose:
    print("Generating a maximum of", options.maxconfs, "per a mol")

if options.etkdg and not Chem.ETKDG:
    print("ETKDB does not appear to be implemented.  Please upgrade RDKit.")
    sys.exit(1)

split = os.path.splitext(output)
if split[1] == ".gz":
    outf = gzip.open(output, "wt+")
    output = split[0]  # strip .gz
else:
    outf = open(output, "w+")


if os.path.splitext(output)[1] == ".pdb":
    sdwriter = Chem.PDBWriter(outf)
else:
    sdwriter = Chem.SDWriter(outf)

if sdwriter is None:
    print("Could not open ".output)
    sys.exit(-1)

for line in smifile:
    toks = line.split()
    smi = toks[0]
    name = " ".join(toks[1:])

    pieces = smi.split(".")
    if len(pieces) > 1:
        smi = max(pieces, key=len)  # take largest component by length
        print("Taking largest component: %s\t%s" % (smi, name))

    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        if options.verbose:
            print(smi)
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            mol.SetProp("_Name", name)

            if options.etkdg:
                cids = Chem.EmbedMultipleConfs(
                    mol, int(options.sample * options.maxconfs), Chem.ETKDG()
                )
            else:
                cids = Chem.EmbedMultipleConfs(
                    mol, int(options.sample * options.maxconfs), randomSeed=options.seed
                )
            if options.verbose:
                print(len(cids), "conformers found")
            cenergy = []
            for conf in cids:
                # not passing confID only minimizes the first conformer
                if options.nomin:
                    cenergy.append(conf)
                elif options.mmff:
                    converged = Chem.MMFFOptimizeMolecule(mol, confId=conf)
                    mp = Chem.MMFFGetMoleculeProperties(mol)
                    cenergy.append(
                        Chem.MMFFGetMoleculeForceField(
                            mol, mp, confId=conf
                        ).CalcEnergy()
                    )
                else:
                    converged = not Chem.UFFOptimizeMolecule(mol, confId=conf)
                    cenergy.append(
                        Chem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy()
                    )
                if options.verbose:
                    print("Convergence of conformer", conf, converged)

            mol = Chem.RemoveHs(mol)
            sortedcids = sorted(cids, key=lambda cid: cenergy[cid])
            if len(sortedcids) > 0:
                mine = cenergy[sortedcids[0]]
            else:
                mine = 0
            if options.rms == 0:
                cnt = 0
                for conf in sortedcids:
                    if cnt >= options.maxconfs:
                        break
                    if (options.energy < 0) or cenergy[conf] - mine <= options.energy:
                        sdwriter.write(mol, conf)
                        cnt += 1
            else:
                written = {}
                for conf in sortedcids:
                    if len(written) >= options.maxconfs:
                        break
                    # check rmsd
                    passed = True
                    for seenconf in written.keys():
                        rms = getRMS(mol, seenconf, conf)
                        if (rms < options.rms) or (
                            options.energy > 0 and cenergy[conf] - mine > options.energy
                        ):
                            passed = False
                            break
                    if passed:
                        written[conf] = True
                        sdwriter.write(mol, conf)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            print("Exception", e)
    else:
        print("ERROR:", smi)

sdwriter.close()
outf.close()
