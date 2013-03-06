# -*- coding: utf-8 -*-

from galaxy.datatypes import data
import logging
from galaxy.datatypes.sniff import *
from galaxy.datatypes.data import get_file_peek
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.binary import Binary
import subprocess
import pybel
import openbabel
openbabel.obErrorLog.StopLogging()

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
        return int(out[0].strip())
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
        return int(out[0].strip())
    except:
        pass
    return 0


class GenericMolFile( data.Text ):
    """
        abstract class for most of the molecule files
    """
    MetadataElement( name="molecules", default=0, desc="Number of molecules", readonly=True, visible=False, optional=True, no_value=0 )

    self.molecule_number = 0

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            if(self.check_filetype(dataset.file_name)):
                dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
                if (self.molecule_number == 1):
                    dataset.blurb = "1 molecule"
                else:
                    dataset.blurb = "%s molecules" % self.molecule_number
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def get_mime(self):
        return 'text/plain'


def read_sdf_records( filename ):
    lines = []
    with open(filename) as handle:
        for line in handle:
            lines.append( line )
            if line.startswith("$$$$"):
                yield lines
                lines = []

class SDF( GenericMolFile ):
    file_ext = "sdf"
    def sniff( self, filename ):
        self.molecule_number = count_special_lines("^\$\$\$\$", filename)
        if self.molecule_number > 0:
            return True
        else:
            return False

    def merge(split_files, output_file):
        """Merge SD files (not yet implemented)."""
        raise NotImplementedError("Merging SD files is not supported at the moment.")

    def split( cls, input_datasets, subdir_generator_function, split_params):
        """Split a SD file (not implemented for now). See read_sdf_records() functions as a start"""
        if split_params is None:
            return None
        raise NotImplementedError("Splitting SD files is not supported at the moment.")


class MOL2( GenericMolFile ):
    file_ext = "mol2"
    def sniff( self, filename ):
        self.molecule_number = count_special_lines("@\<TRIPOS\>MOLECULE", filename)
        if self.molecule_number > 0
            return True
        else:
            return False


class FPS( GenericMolFile ):
    """
        chemfp fingerprint file: http://code.google.com/p/chem-fingerprints/wiki/FPS
    """
    file_ext = "fps"
    def sniff( self, filename ):
        self.molecule_number = count_special_lines('^#', filename, invert = True)
        header = get_headers( filename, sep='\t', count=1 )
        if header[0].strip() == '#FPS1':
            return True
        else:
            return False


class DRF( GenericMolFile ):
    file_ext = "drf"
    def sniff( self, filename ):
        self.molecule_number = count_special_lines('\"ligand id\"', filename)
        if self.molecule_number > 0:
            return True
        else:
            return False


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
        headers = get_headers( filename, sep='\t', count=300 )
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
            residue_numbers = count_special_lines("^HETATM", filename) + count_special_lines("^ATOM", filename)
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            if (residue_numbers == 1):
                dataset.blurb = "1 residue"
            else:
                dataset.blurb = "%s residues" % self.residue_numbers
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

    def sniff( self, filename ):
        self.molecule_number = count_special_lines("^InChI=", filename)
        word_count = count_lines( filename )

        if self.molecule_number != word_count:
            return False

        if self.molecule_number > 0:
            return True
        else:
            return False


class SMILES( Tabular ):
    file_ext = "smi"
    column_names = [ 'SMILES', 'HEADER' ]
    MetadataElement( name="columns", default=2, desc="Number of columns", readonly=True, visible=False )
    MetadataElement( name="column_types", default=['str','str'], param=metadata.ColumnTypesParameter, desc="Column types", readonly=True, visible=False )


    def sniff( self, filename ):
        """
        Its hard or impossible to sniff a SMILES File. All what i know is the
        word_count must be the same as the non-empty line count. And that i can
        try to import the first SMILES and check if it is a molecule.
        """

        self.molecule_number = count_lines( filename, non_empty = True )
        word_count = count_lines( filename )

        if int(self.molecule_number) != word_count:
            return False

        if self.molecule_number > 0:
            for line in open(filename):
                line = line.strip()
                if line:
                    # if we have atoms, we have a molecule
                    try:
                        if len(pybel.readstring('smi', line).atoms) > 0:
                            return True
                        else:
                            return False
                    except:
                        # if convert fails its not a smiles string
                        return False
            return True
        else:
            return False
