#!/usr/bin/env python
"""
    Small library with cheminformatic functions based on openbabel and pgchem.
    Copyright 2012, Bjoern Gruening and Xavier Lucas
"""

import os, sys

try:
    from galaxy import eggs
    eggs.require('psycopg2')
except:
    print 'psycopg2 is not available. It is currently used in the pgchem wrappers, that are not shipped with default CTB'

try:
    import pybel
    import openbabel
except:
    print 'OpenBabel could not be found. A few functions are not available without OpenBabel.'

from multiprocessing import Pool
import glob, tempfile, re
import subprocess

def CountLines( path ):
    out = subprocess.Popen(['wc', '-l', path],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

def grep(pattern, file_obj):
    grepper = re.compile(pattern)
    for line in file_obj:
        if grepper.search(line):
            return True
    return False

def check_filetype(filepath):
    mol = False
    possible_inchi = True
    for line_counter, line in enumerate(open(filepath)):
        if line_counter > 10000:
            break
        if line.find('$$$$') != -1:
            return 'sdf'
        elif line.find('@<TRIPOS>MOLECULE') != -1:
            return 'mol2'
        elif line.find('ligand id') != -1:
            return 'drf'
        elif possible_inchi and re.findall('^InChI=', line):
            return 'inchi'
        elif re.findall('^M\s+END', line):
            mol = True
        # first line is not an InChI, so it can't be an InChI file
        possible_inchi = False

    if mol:
        # END can occures before $$$$, so and SDF file will 
        # be recognised as mol, if you not using this hack'
        return 'mol'
    return 'smi'

def db_connect(args):
    try:
        db_conn = psycopg2.connect("dbname=%s user=%s host=%s password=%s" % (args.dbname, args.dbuser, args.dbhost, args.dbpasswd));
        return db_conn
    except:
        sys.exit('Unable to connect to the db')

ColumnNames = {
    'can_smiles' : 'Canonical SMILES',
    'can' : 'Canonical SMILES',
    'inchi' : 'InChI',
    'inchi_key' : 'InChI key',
    'inchi_key_first' : 'InChI key first',
    'inchi_key_last' : 'InChI key last',
    'molwt' : 'Molecular weight',
    'hbd' : 'Hydrogen-bond donors',
    'donors' : 'Hydrogen-bond donors',
    'hba' : 'Hydrogen-bond acceptors',
    'acceptors' : 'Hydrogen-bond acceptors',
    'rotbonds' : 'Rotatable bonds',
    'logp' : 'logP',
    'psa' : 'Polar surface area',
    'mr' : 'Molecular refractivity',
    'atoms' : 'Number of heavy atoms',
    'rings' : 'Number of rings',
    'set_bits' : 'FP2 bits',
    'id' : 'Internal identifier',
    'tani' : 'Tanimoto coefficient',
    'spectrophore' : 'Spectrophores(TM)',
    'dist_spectrophore' : 'Spectrophores(TM) distance to target',
    'synonym' : 'Entry id',
}

OBDescriptor = {
    'atoms': ["atoms","Number of atoms"],
    'hatoms': ["hatoms","Number of heavy atoms"], # self defined tag hatoms in plugindefines.txt
    'can_smiles' : ["cansmi","Canonical SMILES"],
    'can_smilesNS' : ["cansmiNS","Canonical SMILES without isotopes or stereo"],
    #["abonds","Number of aromatic bonds"],
    #["bonds","Number of bonds"],
    #["dbonds","Number of double bonds"],
    #["formula","Chemical formula"],
    'hba': ["HBA1","Number of Hydrogen Bond Acceptors 1 (JoelLib)"],
    'hba2': ["HBA2","Number of Hydrogen Bond Acceptors 2 (JoelLib)"],
    'hbd': ["HBD","Number of Hydrogen Bond Donors (JoelLib)"],
    'inchi': ["InChI","IUPAC InChI identifier"],
    'inchi_key': ["InChIKey","InChIKey"],
    #["L5","Lipinski Rule of Five"],
    'logp': ["logP","octanol/water partition coefficient"],
    'mr': ["MR","molar refractivity"],
    'molwt': ["MW","Molecular Weight filter"],
    #["nF","Number of Fluorine Atoms"],
    #["s","SMARTS filter"],
    #["sbonds","Number of single bonds"],
    #["smarts","SMARTS filter"],
    #["tbonds","Number of triple bonds"],
    #["title","For comparing a molecule's title"],
    'psa': ["TPSA","topological polar surface area"],
    'rotbonds' : ['ROTATABLE_BOND', 'rotatable bonds'],
}


def print_output(args, rows):
    if args.oformat == 'table':
        outfile = open(args.output, 'w')
        requested_fields = (filter(lambda x: x not in ["[", "]", "'"], args.fetch)).split(', ')
        if args.header:
            outfile.write( 'Identifier\t' + '\t'.join( [ColumnNames[key] for key in requested_fields] ) + '\n' )
        for row in rows:
            outfile.write( row['synonym'] + '\t' + '\t'.join( [str(row[key]) for key in requested_fields] ) + '\n' )

    elif args.oformat in ['sdf', 'mol2']:
        outfile = pybel.Outputfile(args.oformat, args.output, overwrite=True)
        for row in rows:
            try:
                mol = pybel.readstring('sdf', row['mol'])
                if args.oformat == 'sdf':
                    keys = filter(lambda x: x not in ["[", "]", "'"], args.fetch).split(', ')
                    mol.data.update( { ColumnNames['synonym'] : row['synonym'] } )
                    if 'inchi_key' in keys:
                        keys = (', '.join(keys).replace( "inchi_key", "inchi_key_first, inchi_key_last" )).split(', ')
                    [ mol.data.update( { ColumnNames[key] : row[key] } ) for key in keys if key]
                outfile.write(mol)
            except:
                pass
    else:
        outfile = open(args.output, 'w')
        outfile.write( '\n'.join( [ '%s\t%s' % (row[args.oformat], row['synonym'] ) for row in rows ] ) )
    outfile.close()

def pybel_stop_logging():
    openbabel.obErrorLog.StopLogging()

def get_properties_ext(mol):

    HBD = pybel.Smarts("[!#6;!H0]")
    HBA = pybel.Smarts("[$([$([#8,#16]);!$(*=N~O);" +
                       "!$(*~N=O);X1,X2]),$([#7;v3;" +
                       "!$([nH]);!$(*(-a)-a)])]"
                      )
    calc_desc_dict = mol.calcdesc()

    try:
        logp = calc_desc_dict['logP']
    except:
        logp = calc_desc_dict['LogP']

    return {"molwt": mol.molwt,
            "logp": logp,
            "donors": len(HBD.findall(mol)),
            "acceptors": len(HBA.findall(mol)), 
            "psa": calc_desc_dict['TPSA'],
            "mr": calc_desc_dict['MR'],
            "rotbonds": mol.OBMol.NumRotors(),
            "can": mol.write("can").split()[0].strip(), ### tthis one works fine for both zinc and chembl (no ZINC code added after can descriptor string)
            "inchi": mol.write("inchi").strip(),
            "inchi_key": get_inchikey(mol).strip(),
            "rings": len(mol.sssr),
            "atoms": mol.OBMol.NumHvyAtoms(),
            "spectrophore" : OBspectrophore(mol),
           }

def get_inchikey(mol):
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("mol", "inchi")
    conv.SetOptions("K", conv.OUTOPTIONS)
    inchikey = conv.WriteString( mol.OBMol )
    return inchikey

def OBspectrophore(mol):
    spectrophore = pybel.ob.OBSpectrophore()
    # Parameters: rotation angle = 20, normalization for mean and sd, accuracy = 3.0 A and non-stereospecific cages.
    spectrophore.SetNormalization( spectrophore.NormalizationTowardsZeroMeanAndUnitStd )
    return ', '.join( [ "%.3f" % value for value in spectrophore.GetSpectrophore( mol.OBMol ) ] )

def squared_euclidean_distance(a, b):
    try:
        return ((np.asarray( a ) - np.asarray( b ))**2).sum()
    except ValueError:
        return 0

def split_library( lib_path, lib_format = 'sdf', package_size = None ):
    """
        Split a library of compounds. Usage: split_library( lib_path, lib_format, package_size )
        IT currently ONLY WORKS FOR SD-Files
    """
    pack = 1
    mol_counter = 0

    outfile = open('/%s/%s_pack_%i.%s' % ( '/'.join(lib_path.split('/')[:-1]), lib_path.split('/')[-1].split('.')[0], pack, 'sdf'), 'w' )

    for line in open(lib_path, 'r'):
        outfile.write( line )
        if line.strip() == '$$$$':
            mol_counter += 1
            if mol_counter % package_size == 0:
                outfile.close()
                pack += 1
                outfile = open('/%s/%s_pack_%i.%s' % ( '/'.join(lib_path.split('/')[:-1]), lib_path.split('/')[-1].split('.')[0], pack, 'sdf'), 'w' )
                if mol_counter*10 % package_size == 0:
                    print '%i molecules parsed, starting pack nr. %i' % ( mol_counter, pack - 1 )
    outfile.close()

    return True

def split_smi_library( smiles_file, structures_in_one_file ):
    """
        Split a file with SMILES to several files for multiprocessing usage. 
        Usage: split_smi_library( smiles_file, 10 )
    """
    output_files = []
    tfile = tempfile.NamedTemporaryFile(delete=False)

    smiles_handle = open(smiles_file, 'r')
    for count, line in enumerate( smiles_handle ):
        if count % structures_in_one_file == 0 and count != 0:
            tfile.close()
            output_files.append(tfile.name)
            tfile = tempfile.NamedTemporaryFile(delete=False)
        tfile.write(line)
    tfile.close()
    output_files.append(tfile.name)
    smiles_handle.close()
    return output_files


def mp_run(input_path, regex, PROCESSES, function_to_call ):
    paths = []
    [ paths.append(compound_file) for compound_file in glob.glob(str(input_path) + str(regex)) ]
    paths.sort()

    pool = Pool(processes=PROCESSES)
    print 'Process initialized with', PROCESSES, 'processors'
    result = pool.map_async(function_to_call, paths)
    result.get()

    return paths

if __name__ == '__main__':
    print check_filetype(sys.argv[1])

