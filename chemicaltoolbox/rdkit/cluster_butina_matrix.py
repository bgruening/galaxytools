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
import collections
import uuid

from rdkit import rdBase

import cluster_butina
from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils
from pipelines_utils.TsvWriter import TsvWriter

### start field name defintions #########################################

field_Similarity = "Similarity"

### functions #########################################


def MapClusterToMols(clusters, mols):
    i = 0
    for cluster in clusters:
        print("Cluster:", i, cluster)
        for c in cluster:
            # print("Assigning mol",c,"to cluster",i)
            mols[c].SetIntProp("Cluster", i)
        i += 1
    j = 0
    for mol in mols:
        mol.SetIntProp("MolNum", j)
        j += 1
        # print(mol.GetPropsAsDict())

def MapClusterToMols(clusters, mols):
    i = 0
    for cluster in clusters:
        print("Cluster:", i, cluster)
        for c in cluster:
            # print("Assigning mol",c,"to cluster",i)
            mols[c].SetIntProp("Cluster", i)
        i += 1
    j = 0
    for mol in mols:
        mol.SetIntProp("MolNum", j)
        j += 1
        # print(mol.GetPropsAsDict())


def GetDistance(x, y, matrix):
    if x == y:
        return 1.0
    if x > y:
        x2 = y
        y2 = x
    else:
        x2 = x
        y2 = y
        # print("row",",".join(["%.2f" % x for x in matrix[y2-1]]))
    return matrix[y2 - 1][x2]


def GenerateId(cluster, structure):
    row = "%03i" % cluster
    row += "."
    row += "%04i" % structure
    return row


### start main execution #########################################

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit Butina Cluster Matrix')
    parameter_utils.add_default_input_args(parser)
    parser.add_argument('-o', '--output', help="Base name for output file (no extension). If not defined then SDTOUT is used for the structures and output is used as base name of the other files.")
    parser.add_argument('-of', '--outformat', choices=['tsv', 'json'], default='tsv', help="Output format. Defaults to 'tsv'.")
    parser.add_argument('--meta', action='store_true', help='Write metadata and metrics files')
    parser.add_argument('-t', '--threshold', type=float, default=0.7, help='Similarity clustering threshold (1.0 means identical)')
    parser.add_argument('-mt', '--matrixThreshold', type=float, default=0.5, help='Threshold for outputting values (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(cluster_butina.descriptors.keys()), default='rdkit', help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(cluster_butina.metrics.keys()), default='tanimoto', help='similarity metric (default tanimoto)')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    parser.add_argument('-mf', '--metafile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")

    args = parser.parse_args()
    args.logfile.write("Cluster Matrix Args: %s\n\n" % (str(args)))

    descriptor = cluster_butina.descriptors[args.descriptor]
    if descriptor is None:
        raise ValueError('Invalid descriptor name ' + args.descriptor)

    input,suppl = rdkit_utils.default_open_input(args.input, args.informat)

    # handle metadata
    source = "cluster_butina_matrix.py"
    datasetMetaProps = {"source":source, "description": "Butina clustering using RDKit " + rdBase.rdkitVersion}
    clsMappings = {
        "Cluster1": "java.lang.Integer",
        "Cluster2": "java.lang.Integer",
        "ID1": "java.lang.String",
        "ID2": "java.lang.String",
        "M1": "java.lang.String",
        "M2": "java.lang.String",
        "Similarity": "java.lang.Float"
    }
    fieldMetaProps = [{"fieldName":"Cluster", "values": {"source":source, "description":"Cluster number"}}]

    fieldNames = collections.OrderedDict()
    fieldNames['ID1'] = 'ID1'
    fieldNames['ID2'] ='ID2'
    fieldNames['Cluster1'] = 'Cluster1'
    fieldNames['Cluster2'] = 'Cluster2'
    fieldNames['Similarity'] = 'Similarity'
    fieldNames['M1'] = 'M1'
    fieldNames['M2'] = 'M2'


    ### generate fingerprints
    mols = [x for x in suppl if x is not None]
    fps = [descriptor(x) for x in mols]
    input.close()


    ### do clustering
    args.logfile.write("Clustering with descriptor %s metric %s and threshold %d\n\n" % (args.descriptor, args.metric, args.threshold))
    clusters, dists, matrix, = cluster_butina.ClusterFps(fps, args.metric, 1.0 - args.threshold)
    args.logfile.write("Found %d clusters\n\n" % (len(clusters)))

    MapClusterToMols(clusters, mols)

    size = len(matrix)
    args.logfile.write("len(matrix): %d" % (size))
    count = 0
    writer = TsvWriter(args.output, fieldNames)
    writer.writeHeader()

    for i in range(size ):
        args.logfile.write(" element %d has length %d" % (i, len(matrix[i])))
        writer.write(create_values(mols, i, i, 1.0))
        count += 1
        for j in range(len(matrix[i])):
            args.logfile.write(" writing %d %d\n\n" % (i, j))
            dist = matrix[i][j]
            if dist > args.matrixThreshold:
                # the matrix is the lower left segment without the diagonal
                x = j
                y = i + 1
                writer.write(create_values(mols, x, y, dist))
                writer.write(create_values(mols, y, x, dist))
                count += 2
    writer.write(create_values(mols, size, size, 1.0))

    writer.writeFooter()
    writer.close()

    if args.meta:
        values = {'__InputCount__':i, '__OutputCount__':count, 'RDKitCluster':i}
        for key in values:
            args.metafile.write("%s = %s\n\n" % (key, str(values[key])))


def create_values(mols, x, y, dist):
    c1 = mols[x].GetIntProp("Cluster")
    c2 = mols[y].GetIntProp("Cluster")
    bo = collections.OrderedDict()
    bo["uuid"] = str(uuid.uuid4())
    props = {}
    props["Cluster1"] = c1 + 1
    props["Cluster2"] = c2 + 1
    props["ID1"] = GenerateId(c1 + 1, x + 1)
    props["ID2"] = GenerateId(c2 + 1, y + 1)
    props[field_Similarity] = dist
    if mols[x].HasProp("uuid"):
        props["M1"] = mols[x].GetProp("uuid")
    if mols[y].HasProp("uuid"):
        props["M2"] = mols[y].GetProp("uuid")

    return props

if __name__ == "__main__":
    main()