#!/usr/bin/env python

"""
Description of the sticky end concept and what does the marking mean. what is the reaction matrix and so on ....

Copyright 2013 B. Gruening and Hitesh Patel

version: 0.3
"""

import openbabel
openbabel.obErrorLog.StopLogging()
import os
import re
import sys
import pybel
import shutil
import logging
import argparse
import tempfile
import cStringIO
import subprocess
import multiprocessing

from cheminfolib import split_smi_library, check_filetype, CountLines


original_atomic_num_mapping = {89:6,90:7,91:6,92:8,93:6,94:7,95:6,96:7,97:6,98:8,99:6,100:6,101:7,102:6,103:7,104:6,105:7,106:6,107:7,108:16}
atomic_num_regex = re.compile(r'(Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs)')


def unique_files( file_paths, unique_file, temp = False ):
    """
        Concatenate files and makes them unique afterwarts, utilising the GNU `sort` and `cat` tools.
        file_paths is a list of at least one file path.
        unique_file is a file-object, ready to write into it.
    """
    try:
        if not file_paths:
            sys.stdout.write('You need to specify at least one file to run that function.\n')
            sys.exit('FileNotFound %s' % temp)

        # concat the files with python, because cat on the commandline has a limited argument list
        concat = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
        for filename in file_paths:
            shutil.copyfileobj(open(filename, 'rb'), concat)
        concat.close()


        #cat_command = 'cat %s' % ' '.join(file_paths)
        #p1 = subprocess.Popen(cat_command.split(), stdout=subprocess.PIPE) #Set up the echo command and direct the output to a pipe
        p2 = subprocess.Popen(['sort', '-u', '-k', '1,1'], stdin=open(concat.name), stdout=unique_file) #send p1's output to p2
        #p1.stdout.close() #make sure we close the output so p2 doesn't hang waiting for more input
    except Exception, err:
        sys.stderr.write("Error invoking command:\n%s\n" % (err))
        raise
    stdout, stderr = p2.communicate()
    return_code = p2.returncode

    if return_code:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        sys.stderr.write("Return error code %i from command:\n" % return_code)
        sys.stderr.write("%s\n" % cmd)
    unique_file.close()
    os.remove(concat.name)
    return return_code


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
    return list( [(a.atomicnum,  a.idx) for a in molecule.atoms if a.atomicnum in range(89,109) ] )


def replace_markers( mol ):
    """
        Replace all markers with the corresponding original atomicnumbers
    """
    #[atom.OBAtom.SetAtomicNum( original_atomic_num_mapping[ atom.atomicnum ] ) for atom in mol.atoms if atom.atomicnum in range(89,109)]
    for atom in mol.atoms:
        if atom.atomicnum in range(89,109):
            atom.OBAtom.SetAtomicNum( original_atomic_num_mapping[ atom.atomicnum ] )
    return mol


def check_constraint( mol, options ):
    """
        Check if a molecule satisfy the given constraints.
    """
    if mol.molwt > options.molwt:
        return False
    return True


def is_fragment( smiles ):
    """
        Checks if the molecule has a marker defined. If so return True, otherweise False.
        If the mol is a string we assume a SMILES input and convert it to a real pybel molecule
    """
    if type(smiles) != type(''):
        smiles = smiles.write('smi').split()[0]
    hits = len(atomic_num_regex.findall( smiles ))
    return hits


reaction_matrix = read_reaction_matrix()
def merge(mol_one, mol_two, options, mark_fragments = False):
    """
        Main merge function.
        Two pybel molecules are passed to that function and the possible bonds that can be created are looked up.
        Merge mol_two inside mol_one.
    """
    temp_results = list()
    temp_fragments = list()

    replaced_atoms_1 = get_replaced_atoms(mol_one)
    replaced_atoms_2 = get_replaced_atoms(mol_two)

    mol_one_smiles = mol_one.write('smi').split()[0]
    mol_two_smiles = mol_two.write('smi').split()[0]
    possible_bonds = []
    for atomicnum_atom1, idx_atom1 in replaced_atoms_1:
        for atomicnum_atom2, idx_atom2 in replaced_atoms_2:
            if reaction_matrix[ atomicnum_atom1 ][ atomicnum_atom2 ] == 1:
                """ in the special case of alkene we need to take the bond order into account, to create a double bond """
                if atomicnum_atom1 == 99:
                    bond_order = 2
                else:
                    bond_order = 1
                possible_bonds.append( [idx_atom1, idx_atom2, bond_order] )

    no_fragments = True
    no_result = True
    # Create each possible bond and build a new molecule
    for bond_atom1, bond_atom2, bond_order in possible_bonds:
        # merge two molecules together with the concatination of SMILES
        # that step is necessary to reserve the atom-order
        concat_mol = pybel.readstring('smi', "%s.%s" % (mol_one_smiles, mol_two_smiles) )
        mol_one_atom_count = len(mol_one.atoms)

        # build the bond between two atoms
        # the assumption is, when we cancatenate whith the SMILES trick, the atom-order remains the same.
        # For the second molecule all atomnums are increased by the number of atoms in molecule one.
        try:
            concat_mol.OBMol.AddBond( bond_atom1, bond_atom2 + mol_one_atom_count, bond_order)
            # replace only bond making markers with the original atom type
            concat_mol.OBMol.GetAtom(bond_atom1).SetAtomicNum( original_atomic_num_mapping[ concat_mol.OBMol.GetAtom(bond_atom1).GetAtomicNum() ] )
            concat_mol.OBMol.GetAtom(bond_atom2 + mol_one_atom_count).SetAtomicNum( original_atomic_num_mapping[ concat_mol.OBMol.GetAtom(bond_atom2 + mol_one_atom_count).GetAtomicNum() ] )
        except:
            loggin.debug( "Bond could not be created for the following molecule.\nBond1: %s - Bond2: %s - SMILES: %s.%s\n" % (bond_atom1, bond_atom2, mol_one_smiles, mol_two_smiles) )
            return [], []

        """
            An initial fragment has no title or no "Fragment1:" included. At that point add a counter to the header, called 'sticky_ends:'
            We iterate over all fragments that still have that header set and decrease it everytime we see that fragment.
            If the counter is gone or null, we will not process these compound further.
        """
        sticky_ends = False

        if mark_fragments:
            if mol_one.title.find('sticky_ends:') != -1:
                tokens = mol_one.title.split(':')
                sticky_ends = int(tokens[1]) -1
                concat_mol.title = re.sub('sticky_ends:\d*:', 'sticky_ends:%s:' % ( sticky_ends ), mol_one.title)
                #logging.debug('decrease sticky ends to: %s' % sticky_ends)
            elif mol_one.title.strip() or mol_one.title.find('Fragment1:') == -1:
                """
                    set initial fragment counter to the count of all possible replacements -1, 
                    because we already did one replacement
                """
                sticky_ends = max(len(replaced_atoms_1), len(replaced_atoms_2))
                concat_mol.title = 'sticky_ends:%s:' % sticky_ends
                #logging.debug('set sticky end count: %s' % sticky_ends)

        # save the parent fragments to the SMILES header
        #concat_mol.title += 'Fragment1: %s Fragment2: %s' % (mol_one_smiles, mol_two_smiles)


        # concat_mol_smiles = concat_mol.write('can')
        concat_mol_smiles = smi2can( concat_mol )
        true_mol = replace_markers(concat_mol)
        true_mol.title = ''

        if check_constraint( true_mol, options ):
            # temp_results.append( true_mol.write('can') )
            temp_results.append( smi2can( true_mol ) )
            no_result = False
            if is_fragment( concat_mol_smiles.split()[0] ):
                if mark_fragments and sticky_ends == 0:
                    continue
                temp_fragments.append( concat_mol_smiles )
                no_fragments = False

    return temp_results, temp_fragments


def mp_helper(file_one, file_two, mark_fragments = False):
    """
        Helper function for the multiprocessing library.
        Two fragment files gets passed and we merge all against all molecules in that two files.
    """
    results = list()
    fragments = list()
    for mol_one in pybel.readfile( 'smi', file_one ):
        for i,mol_two in enumerate(pybel.readfile( 'smi', file_two )):
            result, fragment = merge(mol_two,mol_one, options, mark_fragments)
            if result:
                results.extend( result )
            if fragment:
                fragments.extend( fragment )

    fragment_return, molecule_return = None, None
    if fragments:
        fragment_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
        fragment_file.writelines( fragments )
        fragment_file.close()
        fragment_return = fragment_file.name
    if results:
        result_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
        result_file.writelines( results )
        result_file.close()
        molecule_return = result_file.name
    return molecule_return, fragment_return


def test(options):
    """
        Small test case.
    """
    
    # change the molwt temporary to a higher level, reset it afterwards
    options_temp = options.molwt
    options.molwt = 10000
    print options_temp
    mol_one = pybel.readstring('can','[Th]c1ccc(cc1)[Ac]=O')
    mol_two = pybel.readstring('can','[Th]c1ccc(cc1)[Ac]=O')

    result, temp = merge(mol_one, mol_two, options)
    result = temp.pop().split('\t')[0]
    assert result == 'O=[Ac]c1ccc(cc1)NC(=O)c1ccc(cc1)[Th]'

    options.molwt = options_temp


def fill_temp_file( file_path, pos ):
    """
        Create temp files for the NxN clustering, so that we can run in parallel.
        We will return two paths. The first with one molecule 
        (the pos-th molecule from the input file) and
        the second with the first pos-molecules.
    """
    out_multiple_molecule = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    out_single_molecule = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    with open( file_path ) as handle:
        for counter, line in enumerate( handle, start = 1):
            out_multiple_molecule.write(line)
            if counter == pos:
                out_single_molecule.write(line)
                break

    out_multiple_molecule.close()
    out_single_molecule.close()
    return out_single_molecule.name, out_multiple_molecule.name


def smi2can( mol ):

    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("can")
    conv.SetOptions("i", conv.OUTOPTIONS)
    conv.SetOptions("c", conv.OUTOPTIONS)
    conv.SetOptions("d", conv.OUTOPTIONS)
    conv.SetOptions("b", conv.OUTOPTIONS)
    can_smiles = conv.WriteString(mol.OBMol)
    return can_smiles


def clean( trash_list ):
    """
        remove temporary files from hard drive
        trash_list is a list with paths, either dirs or files
    """
    while trash_list:
        path = trash_list.pop()
        if os.path.exists( path ):
            if os.path.isdir( path ):
                shutil.rmtree( path )
            else:
                os.remove( path )


fragment_temp_files = list()
result_compounds = list()
def log_result( results ):
    """
        log all results generated by the multiprocessing pool workers
    """
    compounds, fragments = results
    if fragments:
        fragment_temp_files.append( fragments )
    if compounds:
        result_compounds.append( compounds )

result_files_sepcial_mode = list()
def log_result_special_mode( results ):
    if results:
        result_files_sepcial_mode.append(results)


def mp_helper_special_mode(one_fragment, fragment_file_chunks, molecule_dependent_iter_depth = False, chunk_size = 100):

    results = list()
    # one_fragment in a list of one fragment
    fragments = one_fragment

    def mp_callback(res):
        r, f = res
        if r:
            results.append(r)
        if f:
            fragments.append(f)


    for counter, fragment_file_two in enumerate( fragment_file_chunks ):
        combined_fragments = tempfile.NamedTemporaryFile(dir=temp_dir, prefix='tmp_fragments_', delete=False)
        unique_files( fragments, combined_fragments, '1' )

        #logging.debug('fragments (%s)' % ( open(combined_fragments.name).read().strip() ))
        if counter > 0:
            clean( fragments )
        fragments = [combined_fragments.name]
        #logging.debug('Linecount: %s (%s) - loop-counter (%s)' % ( CountLines( combined_fragments.name ), combined_fragments.name, counter))

        splitted_files = split_smi_library( combined_fragments.name, chunk_size )
        logging.info('Fragments to process: %s (%s); Files to process: %s (%s)' % ( CountLines(combined_fragments.name), combined_fragments.name, len(splitted_files), counter) )

        pool = multiprocessing.Pool( options.processors )
        for fragment_file in splitted_files:
            pool.apply_async(mp_helper, args=(fragment_file_two, fragment_file, options.molecule_dependent_iter_depth), callback=mp_callback)
            #log_result( mp_helper(unique_input.name, fragment_file, options.molecule_dependent_iter_depth) )
        pool.close()
        pool.join()
        clean( splitted_files )

    if results:
        combined_results = tempfile.NamedTemporaryFile(dir=temp_dir, prefix='tmp_results_', delete=False)
        unique_files( results, combined_results, '2' )
        clean( results )
        return combined_results.name
    else:
        return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Splits a molecule to several fragments.')

    parser.add_argument("-i", "--input", dest="input_path",
                    required=True,
                    help="Path to the input file.")

    parser.add_argument("-o", "--output", type=argparse.FileType('w'),
                    default=sys.stdout,
                    help="Path to the output file.")

    parser.add_argument("--tempdir", default=tempfile.gettempdir(),
                    help="Specify a temporary directory to use.")

    parser.add_argument("--molwt-cutoff", dest="molwt", type=int,
                    default=1000,
                    help="Limiting the molecular weight of produced molecules. [1000]")

    parser.add_argument("--max-iteration-depth", dest="max_depth", type=int,
                    default=2,
                    help="Max iteration depth. Only used when --mdid-off is specified. [2]")

    parser.add_argument("--mdid-off", dest="molecule_dependent_iter_depth", action="store_false",
                    default=True,
                    help="Disable molecule dependent iteration depth. If specified the options is deactivated and --max-iteration-depth is taken into account.")

    parser.add_argument('-v', action='append_const', const=1, help="Verbose output")

    parser.add_argument('-p', '--processors', type=int, 
        default=multiprocessing.cpu_count())

    options = parser.parse_args()

    log_level = logging.WARNING

    verb_level = 0
    if options.v:
        verb_level = sum(options.v)
    if verb_level >= 2:
        log_level = logging.DEBUG
    elif verb_level == 1:
        log_level = logging.INFO

    logging.basicConfig(level=log_level)
    test(options)

    temp_dir = tempfile.mkdtemp(dir=options.tempdir)
    unique_compounds = set()
    multiple_merge_compounds = set()
    trash_list = list()

    unique_input_raw = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    unique_files( [options.input_path], unique_input_raw )
    unique_input_raw.close()

    unique_input = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    unique_input_non_fragments = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    non_fragment_counter = 0
    for mol in open( unique_input_raw.name ):
        mol = mol.strip()
        if mol:
            if is_fragment(mol):
                unique_input.write('%s\n' % mol)
            else:
                non_fragment_counter += 1
                unique_input_non_fragments.write('%s\n' % mol)
    logging.info('Your input file contains %s non-fragment molecules. These molecules will be merged back in the end-results.' % non_fragment_counter)
    trash_list.append( unique_input_raw.name )
    unique_input.close()
    unique_input_non_fragments.close()

    ##### get a clue out of the input file
    #for line in open(unique_input.name):
    #    if is_fragment(line) >= 4:
    #        print line.strip()
    #        #print is_fragment(line)
    #sys.exit()

    result_compounds.append( unique_input_non_fragments.name )
    """
        special mode
    """
    splitted_files = split_smi_library( unique_input.name, 1 )
    trash_list.extend( splitted_files )

    for counter, fragment_file_one in enumerate( splitted_files ):
        logging.debug('Fragment-file content %s (%s)' % (open(fragment_file_one).read().strip(), counter))
        #pool.apply_async(mp_helper_special_mode, args=( [fragment_file_one], splitted_files, options.molecule_dependent_iter_depth ), callback=log_result_special_mode)
        log_result_special_mode(mp_helper_special_mode( [fragment_file_one], splitted_files, options.molecule_dependent_iter_depth ))


    temp_final_result = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    trash_list.append(temp_final_result.name)
    unique_files( result_files_sepcial_mode, temp_final_result )
    temp_final_result.close()
    for line in open(temp_final_result.name):
        options.output.write( '%s\t\n' % line.strip().split('\t')[0] )

    if not type(options.output.name) != str():
        logging.info('Results are written to %s (%s).' % (options.output.name, CountLines( options.output.name )))
    trash_list.extend( result_compounds )
    trash_list.append( temp_dir )
    logging.info('Cleaning temporary files.')
    clean( trash_list )
    options.output.close()





    '''
    sys.exit('exit normal, after special mode')
    """
        NxN Merge with the same input file
    """
    logging.info('Fragments to start the NxN merge: %s (%s)' % ( CountLines( unique_input.name ), unique_input.name) )
    pool = multiprocessing.Pool( options.processors )
    for counter in range( 1, CountLines( unique_input.name ) + 1):
        single_molecule, multiple_molecule = fill_temp_file( unique_input.name, counter )
        logging.debug('NxN-iteration: %s (%s molecule vs %s molecules)' % (counter, CountLines( single_molecule ), CountLines( multiple_molecule ) ) )
        #log_result( mp_helper(single_molecule, multiple_molecule, options.molecule_dependent_iter_depth) )
        pool.apply_async(mp_helper, args=(single_molecule, multiple_molecule, options.molecule_dependent_iter_depth), callback=log_result)
        trash_list.extend( [single_molecule, multiple_molecule] )

    pool.close()
    pool.join()
    if not fragment_temp_files:
        if not result_compounds:
            sys.exit('No molecules created. Program exits.')
        unique_files( result_compounds, options.output )
        trash_list.extend( result_compounds )
        clean( trash_list )
        options.output.close()
        sys.exit('No further fragments generated. Program exits.')

    combined_fragments = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    unique_files( fragment_temp_files, combined_fragments )
    logging.info('Unique fragments after NxN merge: %s (%s)' % ( CountLines( combined_fragments.name ), combined_fragments.name) )
    trash_list.extend(fragment_temp_files)
    fragment_temp_files = list()
    clean( trash_list )

    """
        #NxM with two different files
    """
    iteration_no = 0
    while True:
        iteration_no += 1
        if iteration_no > options.max_depth and not options.molecule_dependent_iter_depth:
            sys.stdout.write('Max iteration depth is reached. If you want more results please increase --max-depth \n')
            break
        splitted_files = split_smi_library( combined_fragments.name, 50 )
        logging.info('Fragments to process: %s (%s); Files to process: %s' % ( CountLines(combined_fragments.name), combined_fragments.name, len(splitted_files)) )
        trash_list.extend( splitted_files )
        pool = multiprocessing.Pool( options.processors )
        for fragment_file in splitted_files:
            pool.apply_async(mp_helper, args=(unique_input.name, fragment_file, options.molecule_dependent_iter_depth), callback=log_result)
            #log_result( mp_helper(unique_input.name, fragment_file, options.molecule_dependent_iter_depth) )
        pool.close()
        pool.join()

        if options.v > 0:
            # logging is enabled and we can process these intermediate results
            intermediate_results = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
            unique_files( result_compounds, intermediate_results )
            logging.info('Intermediate results: %s (%s)' % ( intermediate_results.name, CountLines(intermediate_results.name) ) )

        if not fragment_temp_files:
            # Leave the while loop, when no further fragments remains.
            break
        combined_fragments = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
        unique_files( fragment_temp_files, combined_fragments )
        trash_list.extend(fragment_temp_files)
        fragment_temp_files = list()
        clean( trash_list )

    unique_files( result_compounds, options.output )
    logging.info('Results are written to %s (%s).' % (options.output.name, CountLines( options.output.name )))
    trash_list.extend( result_compounds )
    trash_list.append( temp_dir )
    trash_list.append( combined_fragments.name )
    logging.info('Cleaning temporary files.')
    clean( trash_list )
    options.output.close()
    '''
