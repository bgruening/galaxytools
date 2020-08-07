#!/usr/bin/env python

import argparse
import os

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


### start main execution #########################################

def main():
    ### command line args defintions #########################################

    ### Define the reactions available
    poised_filter = True
    if poised_filter == True:
        from .poised_filter import Filter
        filter_to_use = Filter()


    parser = argparse.ArgumentParser(description='RDKit rxn process')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-m', '--multi', action='store_true', help='Output one file for each reaction')
    parser.add_argument('-r', '--reaction', choices=filter_to_use.poised_reactions.keys(), help='Name of reaction to be run')
    parser.add_argument('-rl', '--reagent_lib', help="Input SD file, if not defined the STDIN is used")
    parser.add_argument('-rlf', '--reagent_lib_format', choices=['sdf', 'json'], help="Input format. When using STDIN this must be specified.")
    
    args = parser.parse_args()
    utils.log("Screen Args: ", args)

    if not args.output and args.multi:
        raise ValueError("Must specify output location when writing individual result files")

    input, suppl = rdkit_utils.default_open_input(args.input, args.informat)
    reagent_input, reagent_suppl = rdkit_utils.default_open_input(args.reagent_lib, args.reagent_lib_format)
    output, writer, output_base = rdkit_utils.default_open_output(args.output, "rxn_maker", args.outformat)

    i = 0
    count = 0

    if args.multi:
        dir_base = os.path.dirname(args.output)
        writer_dict = filter_to_use.get_writers(dir_base)
    else:
        writer_dict = None
        dir_base = None

    for mol in suppl:
        i+=1
        if mol is None: continue
        # Return a dict/class here - indicating which filters passed
        count = filter_to_use.perform_reaction(mol,args.reaction,reagent_suppl,writer,count)

    utils.log("Created", count, "molecules from a total of ", i, "input molecules")

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