import numpy as np
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
import argparse
from random import randint

def get_mol_from_file(fname, ftype):
    if ftype in ['mol', 'sdf']:
        mol = Chem.MolFromMolFile(options.ligand_path)
    elif ftype == 'pdb':
        mol = Chem.MolFromPDBFile(options.ligand_path)
    elif ftype == 'mol2':
        mol = Chem.MolFromMol2File(options.ligand_path)
    else:
        raise IOError
    if not mol:
        raise IOError
    return mol

def get_params(options):
    mol = get_mol_from_file(options.ligand_path, options.ftype)

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

    optionalvals = ""


    if options.seed != None:
        optionalvals += "seed = " + str(options.seed) + "\n"
    if options.exhaustiveness != None:
        optionalvals += "exhaustiveness = " + str(options.exhaustiveness) + "\n"

    with open(options.output, 'w') as f:
        f.write(
            """
size_x =  {}
size_y =  {}
size_z =  {}
center_x =  {}
center_y =  {}
center_z =  {}
{}""".format(box_dims[0], box_dims[1], box_dims[2], center[0], center[1], center[2], optionalvals)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This tool calculates a confounding box around an input ligand, and uses it to
    generate the input parameters for an autodock vina job. The output file can be fed into
    the autodock vina tool as an alternative to creating the parameter file manually. 
    
    Optionally, you can include a 'buffer' in each of the x,y and z directions (in Å),
    which will be added to the confounding box in the appropriate direction.
    """)
    parser.add_argument('--ligand', dest='ligand_path', help='The input ligand filepath.')
    parser.add_argument('--ftype', dest='ftype', help='Filetype of the input ligand (mol, sdf, pdb, mol2)')
    parser.add_argument('--config', dest='output', help='The output file containing calculated params (txt)')
    parser.add_argument('--exh', dest='exhaustiveness', type=int, help='Exhaustiveness of global search')
    parser.add_argument('--bufx', dest='bufx', default=0, type=float, help='the buffer in the x direction '
                                                                           '(float - in angs.)')
    parser.add_argument('--bufy', dest='bufy', default=0, type=float, help='the buffer in the y direction '
                                                                           '(float - in angs.)')
    parser.add_argument('--bufz', dest='bufz', default=0, type=float, help='the buffer in the z direction '
                                                                           '(float - in angs.)')
    parser.add_argument('--seed', dest='seed', default=None, type=int, help='Random seed for reproducibility')

    options = parser.parse_args()
    get_params(options)
