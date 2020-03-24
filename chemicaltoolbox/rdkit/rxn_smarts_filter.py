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
import os

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

### start main execution #########################################

def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit rxn smarts filter')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-m', '--multi', action='store_true', help='Output one file for each reaction')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('--fer', help='molecule to compare against')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")

    args = parser.parse_args()
    args.logfile.write("Screen Args: %s\n\n" % str(args))

    args.fer = os.path.dirname(args.output)
    if not args.output and args.multi:
        raise ValueError("Must specify output location when writing individual result files")

    ### Define the filter chooser - lots of logic possible
    # SMARTS patterns are defined in poised_filter.py. Currently this is hardcoded.
    # Should make this configurable so that this can be specified by the user at some stage.
    poised_filter = True
    if poised_filter == True:
        from poised_filter import Filter
        filter_to_use = Filter()
    rxn_names = filter_to_use.get_rxn_names()
    args.logfile.write("Using %d reaction filters\n\n" % (len(rxn_names)))


    # handle metadata
    source = "rxn_smarts_filter.py"
    datasetMetaProps = {"source":source, "description": "Reaction SMARTS filter"}
    clsMappings = {}
    fieldMetaProps = []

    for name in rxn_names:
        # this is the Java class type for an array of MoleculeObjects
        clsMappings[name] = "[Lorg.squonk.types.MoleculeObject;"
        fieldMetaProps.append({"fieldName":name, "values": {"source":source, "description":"Sythons from " + name + " reaction"}})

    input, output, suppl, writer, output_base = rdkit_utils.default_open_input_output(
        args.input, args.informat, args.output,
        'rxn_smarts_filter', args.outformat, thinOutput=args.thin,
        valueClassMappings=clsMappings, datasetMetaProps=datasetMetaProps,
        fieldMetaProps=fieldMetaProps)

    i = 0
    count = 0
    writer = rdkit_utils.ThickSDWriter(args.output)

    if args.multi:
        dir_base = os.path.dirname(args.output)
        writer_dict = filter_to_use.get_writers(dir_base)
    else:
        writer_dict = None
        dir_base = None

    for mol in suppl:
        i += 1
        if mol is None: continue
        # Return a dict/class here - indicating which filters passed
        filter_pass = filter_to_use.pass_filter(mol)
        args.logfile.write("Found %d matches\n\n" % (len(filter_pass)))

        if filter_pass:
            props = {}
            count += 1
            for reaction in filter_pass:
                molObjList = []
                # Write the reaction name as a newline separated list of the synthons to the mol object
                # this is used in SDF output
                mol.SetProp(reaction, "\n".join(filter_pass[reaction]))
                # now write to the props is a way that can be used for the JSON output
                for smiles in filter_pass[reaction]:
                    # generate a dict that generates MoleculeObject JSON
                    mo = utils.generate_molecule_object_dict(smiles, "smiles", None)
                    molObjList.append(mo)
                props[reaction] = molObjList

                if args.multi:
                    writer_dict[reaction].write(mol)
                    writer_dict[reaction].flush()
            # write the output.
            # In JSON format the props will override values set on the mol
            # In SDF format the props are ignored so the values in the mol are used
            writer.write(mol, props)
            writer.flush()
    args.logfile.write("Matched %d molecules from a total of %d\n\n" % (count, i))
    if dir_base:
        args.logfile.write("Individual SD files found in: %s\n\n" % (dir_base))

    writer.flush()
    writer.close()
    if input:
        input.close()
    if output:
        output.close()
    # close the individual writers
    if writer_dict:
        for key in writer_dict:
            writer_dict[key].close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RxnSmartsFilter': count})


if __name__ == "__main__":
    main()