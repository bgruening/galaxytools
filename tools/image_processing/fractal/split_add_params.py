import argparse
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Split entries from 'parallelization_list' into individual JSON files"
    )
    parser.add_argument(
        '--init_json', '-i', 
        dest='init_json',
        required=True,
        help='Path to the input JSON file created by cellvoyager_to_ome_zarr_init'
    )
    parser.add_argument(
       '--extra_params','-e',
       dest='extra_params',
       default=None,
       required=False,
       help='Extra parameters as a JSON string'
    )
    parser.add_argument(
        '--outdir','-o',
        dest='outdir',
        required=True,
        help='Directory to write split JSON files'
    )
    args = parser.parse_args()

    with open(args.init_json, 'r') as f:
        data = json.load(f)

    if args.extra_params:
        with open(args.extra_params,'r') as f:
            extra_params=json.load(f)

    plist = data.get('parallelization_list')
    
    for i, entry in enumerate(plist):
        out_data = entry.copy()
        out_data['chunk_sizes']= extra_params
        outfile = os.path.join(args.outdir, f'entry_{i}.json')
        with open(outfile, 'w') as out_f:
            json.dump(out_data, out_f, indent=2)
        
if __name__ == '__main__':
    main()

