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

import argparse, sys, collections

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

### start field name defintions #########################################

field_StructureNum = 'StructureNum'
field_ConformerNum = 'ConformerNum'
field_RMSToCentroid = 'RMSToCentroid'
field_ClusterNum = 'ClusterNum'
field_ClusterCentroid = 'ClusterCentroid'
field_EnergyAbs = 'EnergyAbs'
field_EnergyDelta = 'EnergyDelta'
field_MinimizationConverged = 'MinimizationConverged'


### start function defintions #########################################

def process_mol_conformers(args, mol, i, numConfs, maxAttempts, pruneRmsThresh, clusterMethod, clusterThreshold,
                           minimizeIterations):
    args.logfile.write("Generating conformers for molecule %d\n" % (i+1))
    # generate the conformers
    conformerIds = gen_conformers(mol, numConfs, maxAttempts, pruneRmsThresh, True, True, True)
    conformerPropsDict = {}
    minEnergy = 9999999999999
    for conformerId in conformerIds:
        # energy minimise (optional) and energy calculation
        props = collections.OrderedDict()
        energy = calc_energy(mol, conformerId, minimizeIterations, props)
        if energy and energy < minEnergy:
            minEnergy = energy
        conformerPropsDict[conformerId] = props
    # cluster the conformers
    if clusterMethod:
        rmsClusters = cluster_conformers(mol, clusterMethod, clusterThreshold)
        args.logfile.write("Molecule %d generated %d conformers and %d clusters\n\n" % ((i+1), len(conformerIds), len(rmsClusters)))
        rmsClustersPerCluster = []
        clusterNumber = 0

        for cluster in rmsClusters:
            clusterNumber = clusterNumber + 1
            rmsWithinCluster = align_conformers(mol, cluster)
            for conformerId in cluster:
                props = conformerPropsDict[conformerId]
                props[field_ClusterNum] = clusterNumber
                props[field_ClusterCentroid] = cluster[0] + 1
                idx = cluster.index(conformerId)
                if idx > 0:
                    props[field_RMSToCentroid] = rmsWithinCluster[idx - 1]
                else:
                    props[field_RMSToCentroid] = 0.0
    else:
        args.logfile.write("Molecule %d generated %d conformers\n\n" % ((i+1), len(conformerIds)))

    return conformerPropsDict, minEnergy


def write_conformers(mol, i, conformerPropsDict, minEnergy, writer):
    if mol.HasProp("uuid"):
        parentUuid = mol.GetProp("uuid")
    else:
        parentUuid = None

    for id in range(mol.GetNumConformers()):
        # utils.log("Writing",i,id)

        if mol.HasProp("uuid"):
            parentUuid = mol.GetProp("uuid")
        else:
            parentUuid = None

        for name in mol.GetPropNames():
            mol.ClearProp(name)
        if parentUuid:
            mol.SetProp("SourceMolUUID", parentUuid)
        mol.SetIntProp(field_StructureNum, i + 1)
        mol.SetIntProp(field_ConformerNum, id + 1)
        props = conformerPropsDict[id]
        for key in props:
            mol.SetProp(key, str(props[key]))
        if field_EnergyAbs in props:
            energy = props[field_EnergyAbs]
            if energy:
                mol.SetDoubleProp(field_EnergyAbs, energy)
                mol.SetDoubleProp(field_EnergyDelta, energy - minEnergy)
        writer.write(mol, confId=id)


def gen_conformers(mol, numConfs=1, maxAttempts=1, pruneRmsThresh=0.1, useExpTorsionAnglePrefs=True,
                   useBasicKnowledge=True, enforceChirality=True):
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh,
                                     useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
                                     useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality,
                                     numThreads=0)
    # utils.log("generated",len(ids),"conformers")
    return list(ids)


def calc_energy(mol, conformerId, minimizeIts, props):
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
    if ff:
        # creating ff sometimes fails dues to bad atom types
        ff.Initialize()
        if minimizeIts > 0:
            props[field_MinimizationConverged] = ff.Minimize(maxIts=minimizeIts)
        e = ff.CalcEnergy()
        props[field_EnergyAbs] = e
        return e
    else:
        return None


def cluster_conformers(mol, mode="RMSD", threshold=2.0):
    if mode == "TFD":
        dmat = TorsionFingerprints.GetTFDMatrix(mol)
    else:
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
    # utils.log("generated",len(rms_clusters),"clusters")
    return rms_clusters


def align_conformers(mol, clust_ids):
    rmslist = []
    AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
    return rmslist


### start main execution #########################################

def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit conformers')
    parser.add_argument('-n', '--num', type=int, default=1, help='number of conformers to generate')
    parser.add_argument('-a', '--attempts', type=int, default=0, help='number of attempts')
    parser.add_argument('-r', '--rmsd', type=float, default=1.0, help='prune RMSD threshold')
    parser.add_argument('-c', '--cluster', type=str.lower, choices=['rmsd', 'tfd'],
                        help='Cluster method (RMSD or TFD). If None then no clustering')
    parser.add_argument('-t', '--threshold', type=float,
                        help='cluster threshold (default of 2.0 for RMSD and 0.3 for TFD)')
    parser.add_argument('-e', '--emin', type=int, default=0,
                        help='energy minimisation iterations (default of 0 means none)')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('--smiles', help='input structure as smiles (incompatible with using files or stdin for input)')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")

    args = parser.parse_args()

    if not args.threshold:
        if args.cluster == 'tfd':
            args.threshold = 0.3
        elif args.cluster == 'rmsd':
            args.threshold = 2.0


    args.logfile.write("Conformers Args: %s\n\n" % str(args))

    source = "conformers.py"
    datasetMetaProps = {"source": source, "description": "Conformer generation using RDKit " + rdBase.rdkitVersion}
    clsMappings = {
        "RMSToCentroid": "java.lang.Float",
        "EnergyDelta": "java.lang.Float",
        "EnergyAbs": "java.lang.Float",
        "ConformerNum": "java.lang.Integer",
        "ClusterCentroid": "java.lang.Integer",
        "ClusterNum": "java.lang.Integer",
        "StructureNum": "java.lang.Integer"}
    fieldMetaProps = [
        {"fieldName": "RMSToCentroid",
         "values": {"source": source, "description": "RMS distance to the cluster centroid"}},
        {"fieldName": "EnergyDelta",
         "values": {"source": source, "description": "Energy difference to lowest energy structure"}},
        {"fieldName": "EnergyAbs", "values": {"source": source, "description": "Absolute energy"}},
        {"fieldName": "ConformerNum", "values": {"source": source, "description": "Conformer number"}},
        {"fieldName": "ClusterCentroid",
         "values": {"source": source, "description": "Conformer number of the cluster centroid"}},
        {"fieldName": "ClusterNum", "values": {"source": source, "description": "Cluster number"}},
        {"fieldName": "StructureNum",
         "values": {"source": source, "description": "Structure number this conformer was generated from"}}
    ]

    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        suppl = [mol]
        input = None
        output, writer, output_base = rdkit_utils. \
            default_open_output(args.output, 'conformers', args.outformat,
                                valueClassMappings=clsMappings,
                                datasetMetaProps=datasetMetaProps,
                                fieldMetaProps=fieldMetaProps)
    else:
        input, output, suppl, writer, output_base = rdkit_utils. \
            default_open_input_output(args.input, args.informat, args.output,
                                      'conformers', args.outformat,
                                      valueClassMappings=clsMappings,
                                      datasetMetaProps=datasetMetaProps,
                                      fieldMetaProps=fieldMetaProps)

    # OK, all looks good so we can hope that things will run OK.
    # But before we start lets write the metadata so that the results can be handled.
    # if args.meta:
    #    t = open(output_base + '_types.txt', 'w')
    #    t.write(field_StructureNum + '=integer\n')
    #    t.write(field_StructureNum + '=integer\n')
    #    t.write(field_ConformerNum + '=integer\n')
    #    t.write(field_EnergyAbs + '=double\n')
    #    t.write(field_EnergyDelta + '=double\n')
    #    if args.emin > 0:
    #        t.write(field_MinimizationConverged + '=boolean\n')
    #    if args.cluster:
    #        t.write(field_RMSToCentroid + '=double\n')
    #        t.write(field_ClusterNum + '=integer\n')
    #        t.write(field_ClusterCentroid + '=integer\n')
    #    t.flush()
    #    t.close()

    i = 0
    count = 0
    writer = rdkit_utils.ThickSDWriter(args.output)
    
    for mol in suppl:
        if mol is None: continue
        m = Chem.AddHs(mol)
        conformerPropsDict, minEnergy = process_mol_conformers(args, m, i, args.num, args.attempts, args.rmsd, args.cluster,
                                                               args.threshold, args.emin)
        m = Chem.RemoveHs(m)
        write_conformers(m, i, conformerPropsDict, minEnergy, writer)
        count = count + m.GetNumConformers()
        i += 1

    args.logfile.write("Generated conformers for %d molecules\n" % i)

    if input:
        input.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RDKitConformer': count})


if __name__ == "__main__":
    main()