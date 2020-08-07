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
from copy import copy

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils



def combine_conformers(mols, args):
    """
    Combine the molecules in the input as confomers of a single molecule
    :param mols: A supplier of molecules
    :return: The single molecule with all the conformers
    """
    i = 0
    basemol = None
    writer = rdkit_utils.ThickSDWriter(args.output)

    for mol in mols:
        i += 1
        if mol is None:
            args.logfile.write("WARNING: Molecule %d appears to be invalid\n\n" % i)
            continue
        m = Chem.RemoveHs(mol)
        if not basemol:
            basemol = copy(m)
        else:
            basemol.AddConformer(m.GetConformer())
        writer.write(m)

    args.logfile.write("Number processed: %d New conformer count: %d\n\n" % (i, basemol.GetNumConformers()))


### start main execution #########################################

def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit cluster 3D')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    parameter_utils.add_default_io_args(parser)

    args = parser.parse_args()

    args.logfile.write("Cluster_3d Args: %s\n\n" % str(args))

    source = "cluster_3d.py"
    datasetMetaProps = {"source": source, "description": "Cluster 3D using RDKit " + rdBase.rdkitVersion}
    clsMappings = {
        # "RMSToCentroid": "java.lang.Float",
        # "EnergyDelta": "java.lang.Float",
        # "EnergyAbs": "java.lang.Float",
        # "ConformerNum": "java.lang.Integer",
        # "ClusterCentroid": "java.lang.Integer",
        # "ClusterNum": "java.lang.Integer",
        # "StructureNum": "java.lang.Integer"
    }
    fieldMetaProps = [
        # {"fieldName":"RMSToCentroid",   "values": {"source":source, "description":"RMS distance to the cluster centroid"}},
        # {"fieldName":"EnergyDelta",     "values": {"source":source, "description":"Energy difference to lowest energy structure"}},
        # {"fieldName":"EnergyAbs",       "values": {"source":source, "description":"Absolute energy"}},
        # {"fieldName":"ConformerNum",    "values": {"source":source, "description":"Conformer number"}},
        # {"fieldName":"ClusterCentroid", "values": {"source":source, "description":"Conformer number of the cluster centroid"}},
        # {"fieldName":"ClusterNum",      "values": {"source":source, "description":"Cluster number"}},
        # {"fieldName":"StructureNum",    "values": {"source":source, "description":"Structure number this conformer was generated from"}}
    ]

    input, output, suppl, writer, output_base = rdkit_utils. \
        default_open_input_output(args.input, args.informat, args.output,
                                  'conformers', args.outformat,
                                  valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)


    basemol = combine_conformers(suppl, args)


    if input:
        input.close()
    writer.flush()
    writer.close()
    output.close()

    # if args.meta:
    #     utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RDKitConformer': count})


if __name__ == "__main__":
    main()