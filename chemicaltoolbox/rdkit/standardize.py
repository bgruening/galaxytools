#!/usr/bin/env python

# Copyright 2018 Informatics Matters Ltd.
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

from rdkit import rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils, mol_utils


### functions #########################################


uncharger = rdMolStandardize.Uncharger()


def standardize(mol, neutralize, fragment):
    """
    :param mol: The molecule to standardize
    :param neutralize: Boolean for whether to neutralize the molecule
    :param fragment: The approach for choosing the largest fragment. Either 'hac' or 'mw'. If not specified the whole
    molecule is used.
    :return: The standardized molecule
    """
    mol = rdMolStandardize.Cleanup(mol)

    # We use our own largest fragment picker as the RDKit one behaves slightly differently
    if fragment:
        mol = mol_utils.fragment(mol, fragment)
    if neutralize:
        mol = uncharger.uncharge(mol)
    return mol


### start main execution #########################################

def main():

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='RDKit Standardize')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'], help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('-n', '--neutralize', action='store_true', help='Neutralize the molecule')

    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-t', '--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    
    args = parser.parse_args()
    args.logfile.write("Standardize Args: %s\n" % str(args))

    # handle metadata
    source = "standardize.py"
    datasetMetaProps = {"source":source, "description": "Standardize using RDKit " + rdBase.rdkitVersion}
    clsMappings = {}
    fieldMetaProps = []


    input,output,suppl,writer,output_base = rdkit_utils.\
        default_open_input_output(args.input, args.informat, args.output,
                                  'standardize', args.outformat,
                                  thinOutput=False, valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)

    count = 0
    total = 0
    errors = 0
    writer = rdkit_utils.ThickSDWriter(args.output)
    
    for mol in suppl:
        count += 1
        if mol is None:
            errors += 1
            continue
        m = standardize(mol, args.neutralize, args.fragment_method)
        writer.write(m)
        total += 1

    input.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':count, '__OutputCount__':total, '__ErrorCount__':errors, 'RDKitStandardize':total})

if __name__ == "__main__":
    main()