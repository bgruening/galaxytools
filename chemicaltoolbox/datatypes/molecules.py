# -*- coding: utf-8 -*-

from galaxy.datatypes import data
import logging
from galaxy.datatypes.sniff import get_headers, get_test_fname
from galaxy.datatypes.data import get_file_peek
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.binary import Binary
import subprocess
import os
#import pybel
#import openbabel
#openbabel.obErrorLog.StopLogging()

from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes import metadata

log = logging.getLogger(__name__)

def count_special_lines( word, filename, invert = False ):
    """
        searching for special 'words' using the grep tool
        grep is used to speed up the searching and counting
        The number of hits is returned.
    """
    try:
        cmd = ["grep", "-c"]
        if invert:
            cmd.append('-v')
        cmd.extend([word, filename])
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        return int(out.communicate()[0].split()[0])
    except:
        pass
    return 0

def count_lines( filename, non_empty = False):
    """
        counting the number of lines from the 'filename' file
    """
    try:
        if non_empty:
            out = subprocess.Popen(['grep', '-cve', '^\s*$', filename], stdout=subprocess.PIPE)
        else:
            out = subprocess.Popen(['wc', '-l', filename], stdout=subprocess.PIPE)
        return int(out.communicate()[0].split()[0])
    except:
        pass
    return 0


class GenericMolFile( data.Text ):
    """
        abstract class for most of the molecule files
    """
    MetadataElement( name="number_of_molecules", default=0, desc="Number of molecules", readonly=True, visible=True, optional=True, no_value=0 )

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            if (dataset.metadata.number_of_molecules == 1):
                dataset.blurb = "1 molecule"
            else:
                dataset.blurb = "%s molecules" % dataset.metadata.number_of_molecules
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def get_mime(self):
        return 'text/plain'



class SDF( GenericMolFile ):
    file_ext = "sdf"
    def sniff( self, filename ):
        if count_special_lines("^\$\$\$\$", filename) > 0:
            return True
        else:
            return False

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines("^\$\$\$\$", dataset.file_name)#self.count_data_lines(dataset.file_name)

    def split( cls, input_datasets, subdir_generator_function, split_params):
        """
        Split the input files by molecule records.
        """
        if split_params is None:
            return None

        if len(input_datasets) > 1:
            raise Exception("SD-file splitting does not support multiple files")
        input_files = [ds.file_name for ds in input_datasets]

        chunk_size = None
        if split_params['split_mode'] == 'number_of_parts':
            raise Exception('Split mode "%s" is currently not implemented for SD-files.' % split_params['split_mode'])
        elif split_params['split_mode'] == 'to_size':
            chunk_size = int(split_params['split_size'])
        else:
            raise Exception('Unsupported split mode %s' % split_params['split_mode'])

        def _read_sdf_records( filename ):
            lines = []
            with open(filename) as handle:
                for line in handle:
                    lines.append( line )
                    if line.startswith("$$$$"):
                        yield lines
                        lines = []

        def _write_part_sdf_file( accumulated_lines ):
            part_dir = subdir_generator_function()
            part_path = os.path.join(part_dir, os.path.basename(input_files[0]))
            part_file = open(part_path, 'w')
            part_file.writelines( sdf_lines_accumulated )
            part_file.close()

        try:
            sdf_records = _read_sdf_records( input_files[0] )
            sdf_lines_accumulated = []
            for counter, sdf_record in enumerate( sdf_records, start = 1):
                sdf_lines_accumulated.extend( sdf_record )
                if counter % chunk_size == 0:
                    _write_part_sdf_file( sdf_lines_accumulated )
                    sdf_lines_accumulated = []
            if sdf_lines_accumulated:
                _write_part_sdf_file( sdf_lines_accumulated )
        except Exception,  e:
            log.error('Unable to split files: %s' % str(e))
            raise
    split = classmethod(split)


class MOL2( GenericMolFile ):
    file_ext = "mol2"
    def sniff( self, filename ):
        if count_special_lines("@\<TRIPOS\>MOLECULE", filename) > 0:
            return True
        else:
            return False

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines("@<TRIPOS>MOLECULE", dataset.file_name)#self.count_data_lines(dataset)

    def split( cls, input_datasets, subdir_generator_function, split_params):
        """
        Split the input files by molecule records.
        """
        if split_params is None:
            return None

        if len(input_datasets) > 1:
            raise Exception("MOL2-file splitting does not support multiple files")
        input_files = [ds.file_name for ds in input_datasets]

        chunk_size = None
        if split_params['split_mode'] == 'number_of_parts':
            raise Exception('Split mode "%s" is currently not implemented for MOL2-files.' % split_params['split_mode'])
        elif split_params['split_mode'] == 'to_size':
            chunk_size = int(split_params['split_size'])
        else:
            raise Exception('Unsupported split mode %s' % split_params['split_mode'])

        def _read_sdf_records( filename ):
            lines = []
            start = True
            with open(filename) as handle:
                for line in handle:
                    if line.startswith("@<TRIPOS>MOLECULE"):
                        if start:
                            start = False
                        else:
                            yield lines
                            lines = []
                    lines.append( line )

        def _write_part_sdf_file( accumulated_lines ):
            part_dir = subdir_generator_function()
            part_path = os.path.join(part_dir, os.path.basename(input_files[0]))
            part_file = open(part_path, 'w')
            part_file.writelines( sdf_lines_accumulated )
            part_file.close()

        try:
            sdf_records = _read_sdf_records( input_files[0] )
            sdf_lines_accumulated = []
            for counter, sdf_record in enumerate( sdf_records, start = 1):
                sdf_lines_accumulated.extend( sdf_record )
                if counter % chunk_size == 0:
                    _write_part_sdf_file( sdf_lines_accumulated )
                    sdf_lines_accumulated = []
            if sdf_lines_accumulated:
                _write_part_sdf_file( sdf_lines_accumulated )
        except Exception,  e:
            log.error('Unable to split files: %s' % str(e))
            raise
    split = classmethod(split)



class FPS( GenericMolFile ):
    """
        chemfp fingerprint file: http://code.google.com/p/chem-fingerprints/wiki/FPS
    """
    file_ext = "fps"
    def sniff( self, filename ):
        header = get_headers( filename, sep='\t', count=1 )
        if header[0][0].strip() == '#FPS1':
            return True
        else:
            return False

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines('^#', dataset.file_name, invert = True)#self.count_data_lines(dataset)


class DRF( GenericMolFile ):
    file_ext = "drf"

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines('\"ligand id\"', dataset.file_name, invert = True)#self.count_data_lines(dataset)


class PHAR( GenericMolFile ):
    file_ext = "phar"
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            dataset.blurb = "pharmacophore"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


class PDB( GenericMolFile ):
    file_ext = "pdb"
    def sniff( self, filename ):
        headers = get_headers( filename, sep=' ', count=300 )
        h = t = c = s = k = e = False
        for line in headers:
            section_name = line[0].strip()
            if section_name == 'HEADER':
                h = True
            elif section_name == 'TITLE':
                t = True
            elif section_name == 'COMPND':
                c = True
            elif section_name == 'SOURCE':
                s = True
            elif section_name == 'KEYWDS':
                k = True
            elif section_name == 'EXPDTA':
                e = True

        if h*t*c*s*k*e == True:
            return True
        else:
            return False

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            atom_numbers = count_special_lines("^ATOM", dataset.file_name)
            hetatm_numbers = count_special_lines("^HETATM", dataset.file_name)
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            dataset.blurb = "%s atoms and %s HET-atoms" % (atom_numbers, hetatm_numbers)
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


class grd( data.Text ):
    file_ext = "grd"
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            dataset.blurb = "grids for docking"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


class grdtgz( Binary ):
    file_ext = "grd.tgz"
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = 'binary data'
            dataset.blurb = "compressed grids for docking"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


class InChI( Tabular ):
    file_ext = "inchi"
    column_names = [ 'InChI' ]
    MetadataElement( name="columns", default=2, desc="Number of columns", readonly=True, visible=False )
    MetadataElement( name="column_types", default=['str'], param=metadata.ColumnTypesParameter, desc="Column types", readonly=True, visible=False )
    MetadataElement( name="number_of_molecules", default=0, desc="Number of molecules", readonly=True, visible=True, optional=True, no_value=0 )

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = self.count_data_lines(dataset)

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            if (dataset.metadata.number_of_molecules == 1):
                dataset.blurb = "1 molecule"
            else:
                dataset.blurb = "%s molecules" % dataset.metadata.number_of_molecules
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def sniff( self, filename ):
        """
            InChI files starts with 'InChI='
        """
        inchi_lines = get_headers( filename, sep=' ', count=10 )
        for inchi in inchi_lines:
            if not inchi[0].startswith('InChI='):
                return False
        return True


class SMILES( Tabular ):
    file_ext = "smi"
    column_names = [ 'SMILES', 'TITLE' ]
    MetadataElement( name="columns", default=2, desc="Number of columns", readonly=True, visible=False )
    MetadataElement( name="column_types", default=['str','str'], param=metadata.ColumnTypesParameter, desc="Column types", readonly=True, visible=False )
    MetadataElement( name="number_of_molecules", default=0, desc="Number of molecules", readonly=True, visible=True, optional=True, no_value=0 )

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = self.count_data_lines(dataset)

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            if (dataset.metadata.number_of_molecules == 1):
                dataset.blurb = "1 molecule"
            else:
                dataset.blurb = "%s molecules" % dataset.metadata.number_of_molecules
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


    '''
    def sniff( self, filename ):
        """
        Its hard or impossible to sniff a SMILES File. We can
        try to import the first SMILES and check if it is a molecule, but 
        currently its not possible to use external libraries from the toolshed
        in datatype definition files. TODO
        """
        self.molecule_number = count_lines( filename, non_empty = True )
        word_count = count_lines( filename )

        if self.molecule_number != word_count:
            return False

        if self.molecule_number > 0:
            # test first 3 SMILES
            smiles_lines = get_headers( filename, sep='\t', count=3 )
            for smiles_line in smiles_lines:
                if len(smiles_line) > 2:
                    return False
                smiles = smiles_line[0]
                try:
                    # if we have atoms, we have a molecule
                    if not len( pybel.readstring('smi', smiles).atoms ) > 0:
                        return False
                except:
                    # if convert fails its not a smiles string
                    return False
            return True
        else:
            return False
    '''


if __name__ == '__main__':
    """
        TODO: We need to figure out, how to put example files under /lib/galaxy/datatypes/test/ from a toolshed, so that doctest can work properly.
    """
    inchi = get_test_fname('drugbank_drugs.inchi')
    smiles = get_test_fname('drugbank_drugs.smi')
    sdf = get_test_fname('drugbank_drugs.sdf')
    fps = get_test_fname('50_chemfp_fingerprints_FPS1.fps')
    pdb = get_test_fname('2zbz.pdb')

    print 'SMILES test'
    print SMILES().sniff(smiles), 'smi'
    print SMILES().sniff(inchi)
    print SMILES().sniff(pdb)

    print 'InChI test'
    print InChI().sniff(smiles)
    print InChI().sniff(sdf)
    print InChI().sniff(inchi), 'inchi'

    print 'FPS test'
    print FPS().sniff(smiles)
    print FPS().sniff(sdf)
    f = FPS()
    print f.sniff(fps)

    print 'SDF test'
    print SDF().sniff(smiles)
    print SDF().sniff(sdf), 'sdf'
    print SDF().sniff(fps)

    print 'PDB test'
    print PDB().sniff(smiles)
    print PDB().sniff(sdf)
    print PDB().sniff(fps)
    print PDB().sniff(pdb), 'pdb'
