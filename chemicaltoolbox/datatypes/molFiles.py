# -*- coding: utf-8 -*-

from galaxy.datatypes import data
import logging
from galaxy.datatypes.sniff import *
import commands
import pybel
import openbabel
openbabel.obErrorLog.StopLogging()

from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes import metadata

log = logging.getLogger(__name__)

class GenericMolFile( data.Text ):

    MetadataElement( name="molecules", default=0, desc="Number of molecules", readonly=True, visible=False, optional=True, no_value=0 )

    file_ext = "mol2/sdf/drf"
    def check_filetype( self,filename ):
        self.no_mols = commands.getstatusoutput("grep -c \\$\\$\\$\\$ "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            self.file_ext="sdf"
            return True
        self.no_mols = commands.getstatusoutput("grep -c @\<TRIPOS\>MOLECULE "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            self.file_ext="mol2"
            return True
        self.no_mols = commands.getstatusoutput("grep -c \"ligand id\" "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            self.file_ext="drf"
            return True
        self.no_mols = commands.getstatusoutput("grep -c HEADER "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            self.file_ext="pdb"
            return True
        return False

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            if(self.check_filetype(dataset.file_name)) :
                if (self.no_mols[1] == '1'):
                    dataset.blurb = "1 molecule"
                else:
                    dataset.blurb = "%s molecules" % self.no_mols[1]
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def get_mime(self):
        return 'text/plain'


class GenericMultiMolFile( GenericMolFile ):
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            self.sniff(dataset.file_name)
            if (self.no_mols[1] == '1'):
                dataset.blurb = "1 molecule"
            else:
                dataset.blurb = "%s molecules" % self.no_mols[1]
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

class SDF( GenericMultiMolFile ):
    file_ext = "sdf"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c \\$\\$\\$\\$ "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            return True
        else:
            return False

class MOL2( GenericMultiMolFile ):
    file_ext = "mol2"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c @\<TRIPOS\>MOLECULE "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            return True
        else:
            return False

class FPS( GenericMultiMolFile ):
    file_ext = "fps"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c -v '^#' "+filename)
        with open(filename) as in_handle:
            for line_counter, line in enumerate(in_handle):
                line = line.strip()
                if line.startswith('#FPS1'):
                    return True
                if line_counter > 10:
                    return False

class DRF( GenericMultiMolFile ):
    file_ext = "drf"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c \"ligand id\" "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            return True
        else:
            return False


class PHAR( GenericMultiMolFile ):
    MetadataElement( name="base_name", desc="base name", default='Phar',
    readonly=True, set_in_upload=True)
    file_ext = "phar"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c -v '^#' "+filename)
        return False

class PDB( GenericMolFile ):
    file_ext = "pdb"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c HEADER "+filename)
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            return True
        else:
            return False
    def set_peek( self, dataset, is_multi_byte=False ):
    #def set_peek( self, dataset, line_count=None, is_multi_byte=False ):
        if not dataset.dataset.purged:
            res = commands.getstatusoutput("lib/galaxy/datatypes/countResidues.sh "+dataset.file_name)
            dataset.peek = res[1]
            self.sniff(dataset.file_name)
            if (self.no_mols[1] == '1'):
                dataset.blurb = "1 protein structure"
            else:
                dataset.blurb = "%s protein structures"%self.no_mols[1]
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

class grd ( data.Text ) :
    file_ext = "grd"
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            #dataset.peek = ""
            dataset.blurb = "score-grids for docking"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

class grdtgz ( data.Text ) :
    file_ext = "grd.tgz"
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            #dataset.peek = ""
            dataset.blurb = "compressed score-grids for docking"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


class InChI( GenericMultiMolFile ):
    file_ext = "inchi"
    def sniff( self, filename ):
        self.no_mols = commands.getstatusoutput("grep -c '^InChI=' "+filename)
        word_count = commands.getoutput("wc -w "+filename).split()[0]
        
        if self.no_mols[1] != word_count:
            return False
        
        if (self.no_mols[0] == 0) & (self.no_mols[1] > 0):
            return True
        else:
            return False

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of sequences and the number of data lines
        in dataset.
        """
        if self.max_optional_metadata_filesize >= 0 and dataset.get_size() > self.max_optional_metadata_filesize:
            dataset.metadata.data_lines = None
            dataset.metadata.sequences = None
            return
        #word_count = commands.getoutput("wc -w "+filename).split()[0]
        # word_count are the lines of the file, if word_count and molecule count
        # are the same, that must hold to be an InChI File, then that should be
        # the same number as all non-empty lines
        #dataset.metadata.data_lines = word_count
        #int(commands.getoutput("grep -cve '^\s*$' "+filename))
        #dataset.metadata.molecules = word_count


class SMILES( GenericMultiMolFile ):
    file_ext = "smi"
    def sniff( self, filename ):
        """
        Its hard or impossible to sniff a SMILES File. All what i know is the
        word_count must be the same as the non-empty line count. And that i can
        try to import the first SMILES and check if it is a molecule.
        """
        
        # that corresponds to non-empty line count
        self.no_mols = commands.getstatusoutput("grep -cve '^\s*$' "+filename)
        word_count = int(commands.getoutput("wc -w "+filename).split()[0])
        
        if int(self.no_mols[1]) != word_count:
            return False
        
        if (self.no_mols[0] == 0) & (int(self.no_mols[1]) > 0):
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

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of sequences and the number of data lines
        in dataset.
        """
        if self.max_optional_metadata_filesize >= 0 and dataset.get_size() > self.max_optional_metadata_filesize:
            dataset.metadata.data_lines = None
            dataset.metadata.sequences = None
            return

        #word_count = int(commands.getoutput("wc -w "+filename).split()[0])
        # word_count are the lines of the file, if word_count and molecule count
        # are the same, that must hold to be an InChI File, then that should be
        # the same number as all non-empty lines
        #dataset.metadata.data_lines = word_count
        #dataset.metadata.molecules = word_count


