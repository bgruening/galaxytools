#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import shutil
import networkx as nx
import pybel
import openbabel as ob
import matplotlib.pyplot as plt
import argparse
import datetime

isoprene_smarts = pybel.Smarts('C~C(~C)~C~C')

def mol_to_networkxgraph(mol):
    edges = []
    bondorders = []
    for bond in ob.OBMolBondIter(mol.OBMol):
        #bondorders.append(bond.GetBO())
        if bond.GetBeginAtom().IsCarbon() and bond.GetEndAtom().IsCarbon():
            edges.append( (bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, {'BondOrder': bond.GetBO()}) )
    g = nx.Graph()
    g.add_edges_from(edges)
    return g


def get_largest_component_length(graph):
    # get all connected_components from our molecule graph, its ordered from largest to smallest
    #a = nx.connected_components(g)
    comp = nx.connected_component_subgraphs(graph)
    if comp:
        return len(comp[0])
    else:
        return 0

def consisting_of_isoprens(mol):
    """
        ok we have a real problem here, smarts matching all atoms that satisfy
        the pattern
        but we have overlapping sets of matches, it can be that overlapping
        matches covers the whole molecule, but overlapping isoprens are not
        really a chain of isoprens (the difinition of terpens)
        So: We need a set of __non__ overlapping sets of matches that covers 
        all carbons. It should/can be a terpen and that is what these function
        tries to calculate.
        Example: CC(C)(C)C1=CC(=C2C(=C3C=C(C(=O)C(=C3)C(C)(C)C)C(C)(C)C)C2=C4C=C(C(=O)C(=C4)C(C)(C)C)C(C)(C)C)C=C(C1=O)C(C)(C)C
        (Hitesh sad its not a terpen, yet it has a set of overlapping matches that covers all carbons)
    
    """

    ids = [x.idx for x in mol.atoms if x.atomicnum != 1]
    #match_id_set = set()
    # gehe alle isoprene matches durch und merke dir die Atom-IDs
    #[match_id_set.update(x) for x in isoprene_smarts.findall(mol)]
    g = Grouper()
    # isoprene_smarts.findall(mol) -> [ [2,3,4,5,6], [12,13,14,15,16], [6,7,8,9,19] ]
    [g.join(i) for i in isoprene_smarts.findall(mol)]
    if not g.sets:
        # Wenn nichts gefunden wurde ist es auch kein Terpen, da kein Isprenoid
        # vorhanden ist
        return False
    
    for non_overlapping_match_id_set in g.sets:
        remaining_atom_ids = set(ids).difference( non_overlapping_match_id_set )
        """
            assume you remaining_atom_ids are all non_carbons
            if your assumption is false, change the value to false but if our
            assumption hold for all atoms, than we can return True right now,
            because we only need one set with non-carbons (non-overlapping)
        """
        ret = True
        for id in remaining_atom_ids:
            # Wenn eines der übrigen atome Kohlenstoff ist, handelt es sich nicht
            # um ein Terpene

            # index -1, potential errorr, check later
            if mol.atoms[id-1].atomicnum == 6:
                ret = False
                #break
        if ret:
            return True

    return False


def run( (input_path, id_string, is_dir, inchi, images, write_mol, consisting_of_isoprens_required) ):
    print 'Processing:\t%s' % input_path, id_string
    file_handles = create_file_handles(id_string, is_dir)

    for mol in pybel.readfile('mol', input_path):
        g = mol_to_networkxgraph(mol)
        largest_component = get_largest_component_length(g)
        cid = mol.data[ id_string ]
        #print "%s\t%s" % (cid, largest_component)
        rings = len(mol.sssr)
        if inchi:
            out = "%s\t%s\t%s\n" % (cid, rings, mol.write('inchi').strip())
        else:
            out = "%s\t%s\n" % (cid, rings)

        if largest_component % 5 == 0:
            # a multiple of 5 carbon atoms
            if consisting_of_isoprens_required:
                if not consisting_of_isoprens(mol):
                    continue

            # all remaining atoms are none carbon atoms
            if largest_component > 40:
                file_handles['polyterpene'].write(out)
                try:
                    if images:
                        mol.draw(show=False, update=True, usecoords=True, filename = os.path.join( os.path.dirname( file_handles['polyterpene'].name ), cid + '.png') )
                    if write_mol:
                        mol.write('sdf', os.path.join( os.path.dirname( file_handles['polyterpene'].name ), cid + '.sdf') )
                except:
                    print 'no image for %s' % cid
            elif largest_component == 0:
                #print 'Largest Component is Zero:', cid
                pass
            else:
                file_handles[largest_component].write(out)
                try:
                    if images:
                        mol.draw(show=False, update=True, usecoords=True, filename = os.path.join( os.path.dirname( file_handles[largest_component].name ), cid + '.png') )
                    if write_mol:
                        mol.write('sdf', os.path.join( os.path.dirname( file_handles[largest_component].name ), cid + '.sdf') )
                except:
                    pass#print 'no image for %s' % cid
        else:
            # __not__ a multiple of 5 carbon atoms
            file_handles['potential_terpenoids'].write(out)
            """
            try:
                mol.draw(show=False, update=True, usecoords=True, filename = os.path.join( os.path.dirname( file_handles['potential_terpenoids'].name ), cid + '.png') )
            except:
                print 'no image for %s' % cid
            """

class Grouper():
    """
    
    """
    def __init__(self, init = set()):
        self.sets = [init]
    
    def join(self, *args):
        temp = []
        for arg in args:
            arg_inserted = False
            for pos, myset in enumerate(self.sets):
                #print 'dies ist ein set:', myset
                # if one set is has no intersection with the set arg-set, join the two sets
                if not myset.intersection(arg):
                    temp.append( myset.union(arg) )
                else:
                    if not arg_inserted:
                        temp.append(set(arg))
                        arg_inserted = True
                    temp.append(myset)
        self.sets = temp


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def create_file_handles(key_name, is_dir):
    """
        clean up existing directory structure

        if we should create images or SDF Files we need to create a nested
        directory structure
    """

    if is_dir:
        dir_prefix = "%s_%s" % (key_name, datetime.date.today())
        if not os.path.exists('./%s/%s' % (dir_prefix, 'hemiterpene')):
            for folder_name in ['hemiterpene', 'monoterpene', 'sesquinterpene', 'diterpene', 'sesterterpene', 'triterpene', 'tetraterpene', 'polyterpene', 'potential_terpenoids', 'noname']:
                folder_name = './%s/%s' % (dir_prefix, folder_name)
                #shutil.rmtree(folder_name, ignore_errors = True)
                os.makedirs(folder_name)

        file_handles = {
            5: open('./%s/hemiterpene/hemiterpene.txt' % (dir_prefix), 'a'),
            10: open('./%s/monoterpene/monoterpene.txt' % (dir_prefix), 'a'),
            15: open('./%s/sesquinterpene/sesquinterpene.txt' % (dir_prefix), 'a'),
            20: open('./%s/diterpene/diterpene.txt' % (dir_prefix), 'a'),
            25: open('./%s/sesterterpene/sesterterpene.txt' % (dir_prefix), 'a'),
            30: open('./%s/triterpene/triterpene.txt' % (dir_prefix), 'a'),
            35: open('./%s/noname/noname.txt' % (dir_prefix), 'a'),
            40: open('./%s/tetraterpene/tetraterpene.txt' % (dir_prefix), 'a'),
            'polyterpene': open('./%s/polyterpene/polyterpene.txt' % (dir_prefix), 'a'),
            'potential_terpenoids': open('./%s/potential_terpenoids/potential_terpenoids.txt' % (dir_prefix), 'a'),
        }
    else:
        file_handles = {
            5: open('hemiterpene.txt', 'a'),
            10: open('monoterpene.txt', 'a'),
            15: open('sesquinterpene.txt', 'a'),
            20: open('diterpene.txt', 'a'),
            25: open('sesterterpene.txt', 'a'),
            30: open('triterpene.txt', 'a'),
            35: open('noname.txt', 'a'),
            40: open('tetraterpene.txt', 'a'),
            'polyterpene': open('polyterpene.txt', 'a'),
            'potential_terpenoids': open('potential_terpenoids.txt', 'a'),
        }

    return file_handles

def args_wrapper():

    parser = argparse.ArgumentParser(
        description='Filters a compound library and returns terpens.',
        epilog="TODO!")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')

    parser.add_argument('--infile', nargs='*', type=str)
    """parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'),
        default=sys.stdout)
    """
    parser.add_argument('--processes', metavar='N', type=int,
        help='Number of processes to use. default: 4', default=4)

    parser.add_argument('--key_name', metavar='META_KEY', type=str,
        required=True, help='Specify the key name of the SDF metadata, \
        that contains the molecule identifier')

    parser.add_argument('--isoprens', action='store_true', 
        default=True, help='The structure must be consisting of non overlapping isopren units. default: True')

    parser.add_argument('--images', action='store_true', 
        default=False, help='Creates images for each found terpen. default: False')
    parser.add_argument('--inchi', action='store_true', 
        default=False, help="Creates InChI's for each found terpen. default: False")
    parser.add_argument('--mol', action='store_true', 
        default=False, help="Creates SDF-Files for each found terpen. default: False")

    return parser

def main():
    from multiprocessing import Pool
    args = args_wrapper().parse_args()

    pool = Pool(processes = args.processes)
    collection = [
        #('/home/hitesh/Documents/TerpenesDB/All7STLs.sdf', 'PUBCHEM_COMPOUND_CID', False),
        #('/media/data/pubchem/compound_current/', 'PUBCHEM_COMPOUND_CID', False),
        #('/media/data/pubchem/compound_current/Compound_001250001_001275000.sdf.gz', 'PUBCHEM_COMPOUND_CID', False)
    ]

    args.key_name = args.key_name.replace(' ', '')

    if args.images or args.mol:
        output_should_be_dir = True
    else:
        output_should_be_dir = False

    for inpath in args.infile:
        if os.path.isdir(inpath):
            for filename in os.listdir(inpath):
                #path = os.path.join( '/media/data/pubchem/compound_current/', filename )
                path = os.path.join( inpath, filename )
                if os.path.splitext(path)[-1].lower() in ['.sdf', '.mol', '.gz']:
                    #collection.append( (path, 'PUBCHEM_COMPOUND_CID', False, False, True) )
                    collection.append( (path, args.key_name, output_should_be_dir, args.inchi, args.images, args.mol, args.isoprens) )
        elif os.path.exists(inpath):
            collection.append( (inpath, args.key_name, output_should_be_dir, args.inchi, args.images, args.mol, args.isoprens) )


    #result = pool.apply_async(run, collection[1])
    #print result.get()

    pool.map(run, collection)
    #run('/media/data/chebi/ChEBI_complete.sdf', 'ChEBI ID')

    #run('/media/data/pubchem/compound_current/', 'PUBCHEM_COMPOUND_CID', terpens)


if __name__ == '__main__':
    main()

