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

import argparse, sys, json

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign

import conformers
from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

### start field name defintions #########################################

field_O3DAScore = "O3DAScore"


### start function definitions #########################################

def doO3Dalign(args, i, mol, qmol, use_crippen, threshold, perfect_score, writer, conformerProps=None, minEnergy=None):
    if use_crippen:
        pyO3As = rdMolAlign.GetCrippenO3AForProbeConfs(mol, qmol)
    else:
        pyO3As = rdMolAlign.GetO3AForProbeConfs(mol, qmol)

    if len(pyO3As) == 0:
        return 0

    best_score = 0
    j = 0
    conf_id = -1

    for pyO3A in pyO3As:
        align = pyO3A.Align()
        score = pyO3A.Score()
        if score > best_score:
            best_score = score
            conf_id = j
        j += 1

    if not threshold or perfect_score - best_score < threshold:
        args.logfile.write("Molecule %d %s %.4f" % (i, align, score))
        mol.SetDoubleProp(field_O3DAScore, score)
        if conformerProps and minEnergy:
            eAbs = conformerProps[conf_id][(conformers.field_EnergyAbs)]
            eDelta = eAbs - minEnergy
            if eAbs:
                mol.SetDoubleProp(conformers.field_EnergyAbs, eAbs)
            if eDelta:
                mol.SetDoubleProp(conformers.field_EnergyDelta, eDelta)
        writer.write(mol, confId=conf_id)
        return 1
    return 0

def write_squonk_datasetmetadata(thinOutput, valueClassMappings, datasetMetaProps, fieldMetaProps, size=None):

    meta = {}
    props = {}
    # TODO add created property - how to handle date formats?
    if datasetMetaProps:
        props.update(datasetMetaProps)

    if fieldMetaProps:
        meta["fieldMetaProps"] = fieldMetaProps

    if len(props) > 0:
        meta["properties"] = props

    if valueClassMappings:
        meta["valueClassMappings"] = valueClassMappings
    if thinOutput:
        meta['type'] = 'org.squonk.types.BasicObject'
    else:
        meta['type'] = 'org.squonk.types.MoleculeObject'
    if size != None:
        meta['size'] = size
    s = json.dumps(meta)
    return s
    #metadataFile = open(outputBase + 'check' + '.metadata', 'w')
    #metadataFile.write(s)
    #metadataFile.close()

### start main execution #########################################

def main():
    parser = argparse.ArgumentParser(description='Open3DAlign with RDKit')
    parser.add_argument('--query', help='query molfile')
    parser.add_argument('--qmolidx', help="Query molecule index in SD file if not the first", type=int, default=1)
    parser.add_argument('--crippen', action='store_true', help='Use Crippen (logP) contributions')
    parser.add_argument('-t', '--threshold', type=float, help='score cuttoff relative to alignment of query to itself')
    parser.add_argument('-n', '--num', default=0, type=int,
                        help='number of conformers to generate, if None then input structures are assumed to already be 3D')
    parser.add_argument('-a', '--attempts', default=0, type=int, help='number of attempts to generate conformers')
    parser.add_argument('-r', '--rmsd', type=float, default=1.0, help='prune RMSD threshold for excluding conformers')
    parser.add_argument('-e', '--emin', type=int, default=0,
                        help='energy minimisation iterations for generated conformers (default of 0 means none)')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    parser.add_argument('-m', '--metadata', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the metadata file, when output type is json")
    parameter_utils.add_default_io_args(parser)

    args = parser.parse_args()
    args.logfile.write("o3dAlign Args: %s\n\n" % str(args))

    # TODO - handle molecules with multiple fragments
    # TODO - allow to specify threshold as fraction of perfect score?

    qmol = rdkit_utils.read_single_molecule(args.query, index=args.qmolidx)
    qmol = Chem.RemoveHs(qmol)
    qmol2 = Chem.Mol(qmol)

    source = "conformers.py"
    datasetMetaProps = {"source": source, "description": "Open3DAlign using RDKit " + rdBase.rdkitVersion}
    clsMappings = {"O3DAScore": "java.lang.Float"}
    fieldMetaProps = [
        {"fieldName": "O3DAScore", "values": {"source": source, "description": "Open3DAlign alignment score"}}
    ]
    if args.num > 0:
        # we generate the conformers so will add energy info
        clsMappings["EnergyDelta"] = "java.lang.Float"
        clsMappings["EnergyAbs"] = "java.lang.Float"
        fieldMetaProps.append({"fieldName": "EnergyDelta", "values": {"source": source,
                                                                      "description": "Energy difference to lowest energy conformer"}})
        fieldMetaProps.append(
            {"fieldName": "EnergyAbs", "values": {"source": source, "description": "Absolute energy"}})

    input, output, suppl, writer, output_base = rdkit_utils. \
        default_open_input_output(args.input, args.informat, args.output,
                                  'o3dAlign', args.outformat,
                                  valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)
    if args.outformat == 'json':
        metadata = write_squonk_datasetmetadata(thinOutput=False, valueClassMappings=clsMappings,
                                     datasetMetaProps=datasetMetaProps,
                                     fieldMetaProps=fieldMetaProps, size=None)
        args.metadata.write(metadata)

    if args.crippen:
        pyO3A = rdMolAlign.GetCrippenO3A(qmol2, qmol)
    else:
        pyO3A = rdMolAlign.GetO3A(qmol2, qmol)

    perfect_align = pyO3A.Align()
    perfect_score = pyO3A.Score()
    args.logfile.write('Perfect score: %.3f %.3f %s %d\n\n' % (perfect_align, perfect_score,
                                                               str(Chem.MolToSmiles(qmol, isomericSmiles=True)),
                                                               qmol.GetNumAtoms()))

    i = 0
    count = 0
    total = 0
    errors = 0
    writer = rdkit_utils.ThickSDWriter(args.output)

    for mol in suppl:
        if mol is None:
            i += 1
            continue
        try:
            if args.num > 0:
                mol.RemoveAllConformers()
                conformerProps, minEnergy = conformers.process_mol_conformers(args, mol, i, args.num, args.attempts,
                                                                              args.rmsd, None, None, 0)
                mol = Chem.RemoveHs(mol)
                count += doO3Dalign(args, i, mol, qmol, args.crippen, args.threshold, perfect_score, writer,
                                    conformerProps=conformerProps, minEnergy=minEnergy)
            else:
                mol = Chem.RemoveHs(mol)
                count += doO3Dalign(i, mol, qmol, args.crippen, args.threshold, perfect_score, writer)
            total += mol.GetNumConformers()
        except ValueError as e:
            errors += 1
            args.logfile.write("Molecule %d failed to align: %s\n\n" % (i, e.message))
        i += 1

    input.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, '__ErrorCount__': errors,
                                          'RDKitO3DAlign': total})


if __name__ == "__main__":
    main()