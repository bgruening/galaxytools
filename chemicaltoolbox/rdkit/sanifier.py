#!/usr/bin/env python

# Copyright 2019 Informatics Matters Ltd.
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

### Use MolVS to do tautomer enumeration, stereochemistry enumeration, charge neutralisation.

import sys, argparse

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils
from rdkit import Chem
from rdkit.Chem import AllChem
from sanify_utils import enumerateStereoIsomers, enumerateTautomers, STANDARD_MOL_METHODS


def write_out(mols,count,writer,mol_format,file_format):
    for mol in mols:
        count += 1
        if mol is None: continue
        if mol_format == 'mol_3d':
            AllChem.EmbedMolecule(mol,AllChem.ETKDG())
            fmt = 'mol'
        elif mol_format == 'mol_2d':
            AllChem.Compute2DCoords(mol)
            fmt = 'mol'
        else:
            fmt = 'smiles'

        if file_format == 'sdf':
            writer.write(mol)
        elif file_format == 'json':
            writer.write(mol, format=fmt)

    return count

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit molecule standardizer / enumerator')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-et', '--enumerate_tauts', action='store_true', help='Enumerate all tautomers')
    parser.add_argument('-es', '--enumerate_stereo', action='store_true', help='Enumerate all stereoisomers')
    parser.add_argument('-st', '--standardize', action='store_true', help='Standardize molecules. Cannot  be true if enumerate is on.')
    parser.add_argument('-stm','--standardize_method', default="molvs", choices=STANDARD_MOL_METHODS.keys(), help="Choose the method to standardize.")
    parser.add_argument('-mf', '--mol_format', choices=['smiles', 'mol_2d', 'mol_3d'], help="Format for molecules.")
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    
    args = parser.parse_args()

    args.logfile.write("Sanifier Args: %s\n\n" % str(args))

    if args.standardize and args.enumerate_tauts:
        raise ValueError("Cannot Enumerate Tautomers and Standardize")

    if args.standardize and args.enumerate_stereo:
        raise ValueError("Cannot Enumerate Stereo and Standardize")

    if args.outformat == 'sdf' and args.mol_format == 'smiles':
        raise ValueError("Smiles cannot be used when outputting as SDF")

    if args.standardize:
        getStandardMolecule = STANDARD_MOL_METHODS[args.standardize_method]

    # handle metadata
    source = "sanifier.py"
    datasetMetaProps = {"source":source, "description": "Enumerate tautomers and stereoisomers"}
    clsMappings = {
        "EnumTautIsoSourceMolUUID": "java.lang.String",
        "EnumTautIsoSourceMolIdx": "java.lang.Integer"
    }
    fieldMetaProps = [
        {"fieldName":"EnumTautIsoSourceMolUUID", "values": {"source":source, "description":"UUID of source molecule"}},
        {"fieldName":"EnumTautIsoSourceMolIdx", "values": {"source":source, "description":"Index of source molecule"}}
    ]

    oformat = utils.determine_output_format(args.outformat)

    input,output,suppl,writer,output_base = rdkit_utils. \
        default_open_input_output(args.input, args.informat, args.output,
                                  'sanifier', args.outformat,
                                  thinOutput=False, valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)


    i=0
    count=0
    errors=0
    writer = rdkit_utils.ThickSDWriter(args.output)
    
    for mol in suppl:
        i +=1
        if mol is None: continue

        if args.standardize:
            # we keep the original UUID as there is still a 1-to-1 relationship between the input and outputs
            oldUUID = mol.GetProp("uuid")
            inputCanSmiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            try:
                std = getStandardMolecule(mol)
                outputCanSmiles = Chem.MolToSmiles(std, isomericSmiles=True, canonical=True)
                if oldUUID:
                    std.SetProp("uuid", oldUUID)
                if inputCanSmiles == outputCanSmiles:
                    std.SetProp("Standardized", "False")
                else:
                    std.SetProp("Standardized", "True")
            except:
                errors += 1
                args.logfile.write("Error standardizing %s\n\n" % str(sys.exc_info()[0]))
                std = mol
                std.SetProp("Standardized", "Error")

            count = write_out([std],count,writer,args.mol_format,args.outformat)
        else:
            # we want a new UUID generating as we are generating new molecules
            if mol.HasProp('uuid'):
                parentUuid = mol.GetProp("uuid")
            else:
                parentUuid = None

            results = []


            if args.enumerate_tauts:
                args.logfile.write("Enumerating tautomers\n\n")
                results = enumerateTautomers(mol)
            else:
                results.append(mol)

            if args.enumerate_stereo:
                args.logfile.write("Enumerating steroisomers\n\n")
                mols = results
                results = []
                for m in mols:
                    enumerated = enumerateStereoIsomers(m)
                    results.extend(enumerated)

            for m in results:
                # copy the src mol props
                for name in mol.GetPropNames():
                    m.SetProp(name, mol.GetProp(name))
                # add our new props
                m.ClearProp("uuid")
                m.SetIntProp("EnumTautIsoSourceMolIdx", i)
                if parentUuid:
                    m.SetProp("EnumTautIsoSourceMolUUID", parentUuid)

            count = write_out(results,count,writer,args.mol_format,args.outformat)

    args.logfile.write("Handled %d molecules, resulting in %d outputs" % (i, count))

    writer.flush()
    writer.close()
    input.close()
    output.close()


    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count, '__ErrorCount__':errors , 'RDKitSanify':count })

    return count

if __name__ == "__main__":
    main()