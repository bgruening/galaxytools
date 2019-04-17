import numpy as np
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
import argparse


def get_params(options):
    # make sure we have a mol file by initiating rdkit mol object from input
    mol = Chem.MolFromMolFile(options.ligand_path)
    if not mol:
        raise IOError

    # get rdkit conformer and compute x,y,z of top and bottom corner of confounding cuboid
    conf = mol.GetConformer()
    params = rdShapeHelpers.ComputeConfBox(conf)

    # change tuples to arrays
    coords1 = np.array(params[0])
    coords2 = np.array(params[1])

    # get the centre of the box
    center = np.mean((coords1, coords2), axis=0)

    # calculate box dimensions
    dims = np.abs(coords1 - coords2)

    # optionally add buffers in each direction - expansion
    box_dims = [dims[0] + options.bufx, dims[1] + options.bufy, dims[2] + options.bufz]

    with open(options.output, 'w') as f:
        f.write(
            """
size_x =  {}
size_y =  {}
size_z =  {}
center_x =  {}
center_y =  {}
center_z =  {}
num_modes = 9999
energy_range = 9999
exhaustiveness = {}
cpu = 4
seed = 1
            """.format(box_dims[0], box_dims[1], box_dims[2], center[0], center[1], center[2], options.exhaustiveness)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand', dest='ligand_path')
    parser.add_argument('--config', dest='output')
    parser.add_argument('--exh', dest='exhaustiveness', default=10, type=float)
    parser.add_argument('--bufx', dest='bufx', default=0, type=float)
    parser.add_argument('--bufy', dest='bufy', default=0, type=float)
    parser.add_argument('--bufz', dest='bufz', default=0, type=float)

    options = parser.parse_args()
    get_params(options)
