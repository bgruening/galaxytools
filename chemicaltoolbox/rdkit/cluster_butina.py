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

import argparse, sys

from rdkit import DataStructs, rdBase
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Cluster import Butina

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils, mol_utils

descriptors = {
    #'atompairs':   lambda m: Pairs.GetAtomPairFingerprint(m),
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,2,1024),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3,1024),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m),
    #'topo':        lambda m: Torsions.GetTopologicalTorsionFingerprint(m)
}

metrics = {
    'asymmetric':DataStructs.AsymmetricSimilarity,
    'braunblanquet':DataStructs.BulkBraunBlanquetSimilarity,
    'cosine':DataStructs.BulkCosineSimilarity,
    'dice': DataStructs.BulkDiceSimilarity,
    'kulczynski':DataStructs.BulkKulczynskiSimilarity,
    'mcconnaughey':DataStructs.BulkMcConnaugheySimilarity,
    #'onbit':DataStructs.OnBitSimilarity,
    'rogotgoldberg':DataStructs.BulkRogotGoldbergSimilarity,
    'russel':DataStructs.BulkRusselSimilarity,
    'sokal':DataStructs.BulkSokalSimilarity,
    'tanimoto': DataStructs.BulkTanimotoSimilarity
    #'tversky': DataStructs.TverskySimilarity
}

### start field name defintions #########################################

field_Cluster = "Cluster"

### functions #########################################

def ClusterFps(fps, metric, cutoff):

    # first generate the distance matrix:
    dists = []
    # dist is the part of the distance matrix below the diagonal as an array:
    # 1.0, 2.0, 2.1, 3.0, 3.1, 3.2 ...
    nfps = len(fps)
    matrix = []
    for i in range(1,nfps):

        func = metrics[metric]
        sims = func(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
        matrix.append(sims)

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs,dists,matrix

def ClustersToMap(clusters):
    d = {}
    i = 0
    for c in clusters:
        for id in c:
            d[id] = i
        i += 1
    return d

def FetchScore(idx, mols, field, maximise):
    if maximise:
        return 0 - mols[idx].GetDoubleProp(field)
    else:
        return mols[idx].GetDoubleProp(field)


def SelectDiverseSubset(mols, clusters, distances, count, field, maximise, score, quiet, logfile):
    total = len(mols)
    num_clusters = len(clusters)
    pickedList = []
    clustersList = []
    for i in range(0, num_clusters):
        pickedList.append([])
        if field:
            filteredByValue = [x for x in clusters[i] if mols[x].HasProp(field)]
            sortedByValue = sorted(filteredByValue, key=lambda idx: FetchScore(idx, mols, field, maximise))
            clustersList.append(sortedByValue)
        else:
            allRecords = [x for x in clusters[i]]
            clustersList.append(allRecords)

    totalIter = 0
    clusterIter = 0
    pickedCount = 0

    while totalIter < total and pickedCount < count:
        clusterNum = totalIter % num_clusters
        clus = clustersList[clusterNum]
        pick = pickedList[clusterNum]
        if len(clus) > 0:
            # remove that item from the cluster so that it's not tried again
            molIndex = clus.pop(0)
            if len(pick) == 0: # first time for this cluster
                pick.append(molIndex)
                pickedCount += 1
                clusterIter += 1
                if not quiet:
                    logfile.write("Cluster %d initialised with %d\n\n" % (clusterNum, molIndex))
            else:
                closestDist = GetClosestDistance(distances, molIndex, pick)
                if closestDist < score:
                    pick.append(molIndex)
                    pickedCount += 1
                    clusterIter += 1
                    if not quiet:
                        logfile.write("Cluster %d added %d with score %f\n\n" % (clusterNum, molIndex, closestDist))
                elif not quiet:
                    logfile.write("Cluster %d discarded %d with score %f\n\n" % (clusterNum, molIndex, closestDist))
        else: # cluster has been exhausted
            clusterIter += 1

        totalIter += 1

        logfile.write("Picked %d using %d iterations\n\n" % (pickedCount, totalIter))
    return pickedList

def GetDistance(idx1, idx2, distances):
    idx = 0
    for i in range(1, idx1):
        idx += i
    idx += idx2
    d = distances[idx]
    return d

def GetClosestDistance(distances, molIndex, compareTo):
    best = 0
    for i in compareTo:
        d = GetDistance(molIndex, i, distances)
        if best < d:
            best = d
    return best

### start main execution #########################################

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit Butina Cluster')
    parser.add_argument('-t', '--threshold', type=float, default=0.7, help='similarity clustering threshold (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit', help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto', help='similarity metric (default tanimoto)')
    parser.add_argument('-n', '--num', type=int, help='maximum number to pick for diverse subset selection')
    parser.add_argument('-e', '--exclude', type=float, default=0.9, help='threshold for excluding structures in diverse subset selection (1.0 means identical)')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'], default='hac', help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('--output-fragment', action='store_true', help='Output the biggest fragment rather than the original molecule')
    parser.add_argument('-f', '--field', help='field to use to optimise diverse subset selection')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--min', action='store_true', help='pick lowest value specified by the --field option')
    group.add_argument('--max', action='store_true', help='pick highest value specified by the --field option')

    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")

    args = parser.parse_args()
    args.logfile.write("Cluster Args: %s\n\n" % (str(args)))

    descriptor = descriptors[args.descriptor]
    if descriptor is None:
        raise ValueError('Invalid descriptor name ' + args.descriptor)

    if args.field and not args.num:
        raise ValueError('--num argument must be specified for diverse subset selection')
    if args.field and not (args.min or args.max):
        raise ValueError('--min or --max argument must be specified for diverse subset selection')

    # handle metadata
    source = "cluster_butina.py"
    datasetMetaProps = {"source":source, "description": "Butina clustering using RDKit " + rdBase.rdkitVersion}
    clsMappings = {"Cluster": "java.lang.Integer"}
    fieldMetaProps = [{"fieldName":"Cluster", "values": {"source":source, "description":"Cluster number"}}]

    input,output,suppl,writer,output_base = rdkit_utils.\
        default_open_input_output(args.input, args.informat, args.output,
                                  'cluster_butina', args.outformat,
                                  thinOutput=False, valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)

    ### fragment and generate fingerprints
    mols = []
    fps = []
    errs = mol_utils.fragmentAndFingerprint(suppl, mols, fps, descriptor, fragmentMethod=args.fragment_method, outputFragment=args.output_fragment, quiet=False)


    input.close()

    ### do clustering
    args.logfile.write("Clustering with descriptor %s metric %s and threshold %d\n\n" % (args.descriptor, 
                                                                                         args.metric, args.threshold))
    clusters, dists, matrix = ClusterFps(fps, args.metric, 1.0 - args.threshold)

    args.logfile.write("Found %d clusters\n\n" % (len(clusters)))

    ### generate diverse subset if specified
    if args.num:
        args.logfile.write("Generating diverse subset\n\n")
        # diverse subset selection is specified
        finalClusters = SelectDiverseSubset(mols, clusters, dists, args.num, args.field, args.max, args.exclude, 
                                            False, args.logfile)
    else:
        finalClusters = clusters

        args.logfile.write("Found %d clusters\n\n" % (len(finalClusters)))
    lookup = ClustersToMap(finalClusters)

    ### write the results
    i = 0
    result_count = 0
    writer = rdkit_utils.ThickSDWriter(args.output)
    
    for mol in mols:
        if i in lookup:
            if args.thin:
                rdkit_utils.clear_mol_props(mol, ["uuid"])
            cluster = lookup[i]
            mol.SetIntProp(field_Cluster, cluster)
            writer.write(mol)
            result_count += 1
        i += 1

    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        status_str = str(result_count) + ' results from ' + str(len(finalClusters)) + ' clusters'
        utils.write_metrics(output_base, {'__StatusMessage__':status_str, '__InputCount__':i, '__OutputCount__':result_count, 'RDKitCluster':i})

if __name__ == "__main__":
    main()