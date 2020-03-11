#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse, logging, time, sys

from rdkit import DataStructs, rdBase, SimDivFilters
from rdkit.Chem import AllChem, MACCSkeys

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import mol_utils, rdkit_utils

descriptors = {
    #'atompairs':   lambda m: Pairs.GetAtomPairFingerprint(m),
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,2,1024),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3,1024)
}


### functions #########################################


def performPick(fpBitVector, howManyToPick, similarityThreshold, firstPicks):
    picker = SimDivFilters.MaxMinPicker()
    # LazyPickWithThreshold( (MaxMinPicker)self, (AtomPairsParameters)distFunc, (int)poolSize, (int)pickSize, (float)threshold [, (AtomPairsParameters)firstPicks=() [, (int)seed=-1]]) -> tuple :
    picks, thresh = picker.LazyBitVectorPickWithThreshold(fpBitVector, len(fpBitVector), howManyToPick, similarityThreshold, firstPicks=firstPicks)
    return picks, thresh


### start main execution #########################################

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit Butina Cluster')
    parser.add_argument('-t', '--threshold', type=float, default=0.0, help='similarity threshold (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='morgan2', help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-n', '--num', type=int, help='maximum number to pick for diverse subset selection')
    parser.add_argument('-s', '--seed-molecules', help='optional file containing any seed molecules that have already been picked')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'], default='hac', help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('--output-fragment', action='store_true', help='Output the biggest fragment rather than the original molecule')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'), 
        default=sys.stdout, help="Path to the log file, default it sdtout")
    parameter_utils.add_default_io_args(parser)

    args = parser.parse_args()

    #logfile = open(args.output+'.txt', 'w+')

    args.logfile.write("MaxMinPicker Args: %s\n" % str(args))

    descriptor = descriptors[args.descriptor]
    if descriptor is None:
        raise ValueError('No descriptor specified')

    if not args.num and not args.threshold:
        raise ValueError('--num or --threshold arguments must be specified, or both')

    # handle metadata
    source = "max_min_picker.py"
    datasetMetaProps = {"source":source, "description": "MaxMinPicker using RDKit " + rdBase.rdkitVersion}

    ### generate fingerprints
    fps = []
    mols = []
    errors = 0

    # first the initial seeds, if specified
    firstPicks = []
    num_seeds = 0
    if args.seed_molecules:
        seedsInput,seedsSuppl = rdkit_utils.default_open_input(args.seed_molecules, None)
        start = time.time()
        errors += mol_utils.fragmentAndFingerprint(seedsSuppl, mols, fps, descriptor, fragmentMethod=args.fragment_method, outputFragment=args.output_fragment, quiet=args.quiet)
        end = time.time()
        seedsInput.close()
        num_seeds = len(fps)
        args.logfile.write("Read %d fingerprints for seeds in %f secs, %d errors\n" % (len(fps), end-start, errors))
        firstPicks = list(range(num_seeds))

    # now the molecules to pick from
    input,output,suppl,writer,output_base = rdkit_utils.default_open_input_output(args.input, args.informat, args.output, 'max_min_picker',
                                                                            args.outformat, datasetMetaProps=datasetMetaProps)
    # reset the mols list as we don't need the seeds, only the candidates
    mols = []
    start = time.time()
    errs = mol_utils.fragmentAndFingerprint(suppl, mols, fps, descriptor, fragmentMethod=args.fragment_method, outputFragment=args.output_fragment, quiet=args.quiet)
    end = time.time()
    errors += errs

    input.close()
    num_fps = len(fps)
    num_candidates = num_fps - num_seeds
    args.logfile.write("Read %d fingerprints for candidates in %f secs, %d errors\n" % (num_candidates, end-start, errs))

    if not args.num:
        num_to_pick = num_candidates
    elif args.num > num_candidates:
        num_to_pick = num_candidates
        args.logfile.write("WARNING: --num argument (, %d, ) is larger than the total number of candidates (, %d, )\
         - resetting to %d\n" %(args.num, num_candidates, num_candidates))
    else:
        num_to_pick = args.num

    ### do picking
    args.logfile.write("MaxMinPicking with descriptor %s, and threshold %.1f, %d seeds, %d candidates, %d total\n" % (args.descriptor, args.threshold, num_seeds, num_candidates, num_fps))
    start = time.time()
    picks, thresh = performPick(fps, num_to_pick + num_seeds, args.threshold, firstPicks)
    end = time.time()
    num_picks = len(picks)

    args.logfile.write("Found %d molecules in %f secs, final threshold %.1f\n" % (num_picks, end-start, thresh ))
    args.logfile.write("Picks: %s\n" % str(list(picks[num_seeds:])))
    del fps

    # we want to return the results in the order they were in the input so first we record the order in the pick list
    indices = {}
    i = 0
    for idx in picks[num_seeds:]:
        indices[idx] = i
        i += 1
    # now do the sort
    sorted_picks = sorted(picks[num_seeds:])
    # now write out the mols in the correct order recording the value in the pick list as the PickIndex property
    i = 0
    writer = rdkit_utils.ThickSDWriter(args.output)
    for idx in sorted_picks:
        mol = mols[idx - num_seeds] # mols array only contains the candidates
        mol.SetIntProp("PickIndex", indices[idx] + 1)
        writer.write(mol)
        i += 1
    args.logfile.write("Output %d molecules\n" % i)

    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        metrics = {}
        status_str = "{} compounds picked. Final threshold was {}.".format(i, thresh)
        if errors > 0:
            metrics['__ErrorCount__'] = errors
            status_str = status_str + " {} errors.".format(errors)

        metrics['__StatusMessage__'] = status_str
        metrics['__InputCount__'] = num_fps
        metrics['__OutputCount__'] = i
        metrics['RDKitMaxMinPicker'] = num_picks

        utils.write_metrics(output_base, metrics)

if __name__ == "__main__":
    main()
