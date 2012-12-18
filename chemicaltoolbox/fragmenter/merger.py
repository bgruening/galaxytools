#!/usr/bin/env python
"""
Description of the sticky end concept and what does the marking mean. what is the reaction matrix and so on ....

Copyright 2012 B. Gruening and Hitesh Patel
"""
import pybel
import argparse

original_atomic_num_mapping = {89:6,90:7,91:6,92:8,93:6,94:7,95:6,96:7,97:6,98:8,99:6,100:6,101:7,102:6,103:7,104:6,105:7,106:6,107:7,108:16}

def read_reaction_matrix():
    """
        Reads the reaction matrix from a external file.
        This matrix specifies which sticky ends are allowed to be merged together.
        We read a comma separated matrix and stores them into a dictionary of dictionaries to access every cell easyly.
    """

    reaction_matrix_raw = """
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0
    1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0
    1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0
    1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0
    """
    reaction_matrix = {}
    row_name = 89
    for line in reaction_matrix_raw.split('\n'):
        line = line.strip()
        if not line:
            continue
        b = dict(zip(range(89,109), map(int, line.split(',') ) ))
        reaction_matrix[row_name] = b
        row_name += 1
    return reaction_matrix

def get_replaced_atoms(molecule):
    """
        Extracts a tuple of (atomic_number, idx) from each atom and creates a dictionary out of it.
        Key = atomic_number and Value = idx of the molecule.
        Only atoms with a atomic_number between 89 and 109 are processed. Theses are the marked atoms, the sticky ends.
    """
    return dict( [(a.atomicnum,  a.idx) for a in molecule.atoms if a.atomicnum in range(89,109) ] )

def merge(mol_one, mol_two):
    replaced_atoms_1 = get_replaced_atoms(mol_one)
    replaced_atoms_2 = get_replaced_atoms(mol_two)
    reaction_matrix = read_reaction_matrix()

    possible_bonds = []
    for atomicnum_atom1, idx_atom1 in replaced_atoms_1.items():
        for atomicnum_atom2, idx_atom2 in replaced_atoms_2.items():
            if reaction_matrix[ atomicnum_atom1 ][ atomicnum_atom2 ] == 1:
                possible_bonds.append( [idx_atom1, idx_atom2] )

    merged_molecules = []
    # Create each possible bond and build a new molecule
    for bond_atom1, bond_atom2 in possible_bonds:
        # merge two molecules together with the concatination of SMILES
        concat_mol = pybel.readstring('smi', "%s.%s" % (mol_one.write('smi').split()[0], mol_two.write('smi').split()[0]) )
        mol_one_atom_count = len(mol_one.atoms)
        # build the bond between two atoms
        # the assumption is, when we cancatenate whith the SMILES trick, the atom-order remains the same.
        # For the second molecule all atomnums are increased by the number of atoms in molecule one.
        concat_mol.OBMol.AddBond( bond_atom1, bond_atom2 + mol_one_atom_count, 1)
        # replace all markers with the original atom type
        for atom in concat_mol.atoms:
            if atom.atomicnum in range(89,109):
                atom.OBAtom.SetAtomicNum( original_atomic_num_mapping[ atom.atomicnum ] )
        concat_mol.title = 'Fragment1: %s Fragment2: %s' % (mol_one.write('can').split()[0], mol_two.write('can').split()[0])
        merged_molecules.append( concat_mol.write('can') )
    return merged_molecules

def test():
    mol_one = pybel.readstring('can','[Th]c1ccc(cc1)[Ac]=O')
    mol_two = pybel.readstring('can','[Th]c1ccc(cc1)[Ac]=O')
    result = merge(mol_one, mol_two)
    assert result[0] == 'O=Cc1ccc(cc1)NC(=O)c1ccc(cc1)N\tFragment1: O=[Ac]c1ccc(cc1)[Th] Fragment2: O=[Ac]c1ccc(cc1)[Th]\n'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Splits a molecule to several fragments.')

    parser.add_argument("-i", "--input", dest="input_path",
                    required=True,
                    help="Path to the input file.")

    parser.add_argument("-o", "--output", dest="output_path",
                    required=True,
                    help="Path to the output file.")

    parser.add_argument("--input-format", dest="iformat",
                    help="Input format. It must be supported by openbabel.")

    parser.add_argument("--output-format", dest="oformat",
                    default="can",
                    help="Output format. It must be supported by openbabel.")

    options = parser.parse_args()
    #output = pybel.Outputfile(options.oformat, options.output_path, overwrite=True)
    output = open(options.output_path, 'w+')
    
    if not options.iformat:
        options.iformat = check_filetype(options.input_path)

    unique_compounds = set()
    no_merge = 0
    for i, mol_one in enumerate(pybel.readfile( options.iformat, options.input_path )):
        print '--',i
        for y, mol_two in enumerate( pybel.readfile( options.iformat, options.input_path ) ):
            if y <= i:
                continue
            #print y,
            merged_molecules = merge(mol_one, mol_two)
            for molecule_string in merged_molecules:
                if not molecule_string.split()[0] in unique_compounds:
                    unique_compounds.add( molecule_string.split()[0] ) #unique_compounds[molecule_string.split()[0]] = molecule_string
                    output.write( molecule_string )
            else:
                no_merge += 1

    print no_merge
    print len(unique_compounds)

    test()

