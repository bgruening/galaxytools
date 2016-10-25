# -*- coding: utf-8 -*-

from galaxy.datatypes import data
import logging
from galaxy.datatypes.sniff import get_headers, get_test_fname
from galaxy.datatypes.data import get_file_peek
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.binary import Binary
from galaxy.datatypes.xml import GenericXml
import subprocess
import os

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


class MOL( GenericMolFile ):
    file_ext = "mol"
    def sniff( self, filename ):
        if count_special_lines("^M\s*END", filename) == 1:
            return True
        else:
            return False

    def set_meta( self, dataset, **kwd ):
        """
        Set the number molecules, in the case of MOL its always one.
        """
        dataset.metadata.number_of_molecules = 1


class SDF( GenericMolFile ):
    file_ext = "sdf"
    def sniff( self, filename ):
        if count_special_lines("^\$\$\$\$", filename) > 0:
            return True
        else:
            return False

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of molecules in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines("^\$\$\$\$", dataset.file_name)

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
            part_file.writelines( accumulated_lines )
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

        def _read_mol2_records( filename ):
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

        def _write_part_mol2_file( accumulated_lines ):
            part_dir = subdir_generator_function()
            part_path = os.path.join(part_dir, os.path.basename(input_files[0]))
            part_file = open(part_path, 'w')
            part_file.writelines( accumulated_lines )
            part_file.close()

        try:
            mol2_records = _read_mol2_records( input_files[0] )
            mol2_lines_accumulated = []
            for counter, mol2_record in enumerate( mol2_records, start = 1):
                mol2_lines_accumulated.extend( mol2_record )
                if counter % chunk_size == 0:
                    _write_part_mol2_file( mol2_lines_accumulated )
                    mol2_lines_accumulated = []
            if mol2_lines_accumulated:
                _write_part_mol2_file( mol2_lines_accumulated )
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
        dataset.metadata.number_of_molecules = count_special_lines('^#', dataset.file_name, invert = True)


    def split( cls, input_datasets, subdir_generator_function, split_params):
        """
        Split the input files by fingerprint records.
        """
        if split_params is None:
            return None

        if len(input_datasets) > 1:
            raise Exception("FPS-file splitting does not support multiple files")
        input_files = [ds.file_name for ds in input_datasets]

        chunk_size = None
        if split_params['split_mode'] == 'number_of_parts':
            raise Exception('Split mode "%s" is currently not implemented for MOL2-files.' % split_params['split_mode'])
        elif split_params['split_mode'] == 'to_size':
            chunk_size = int(split_params['split_size'])
        else:
            raise Exception('Unsupported split mode %s' % split_params['split_mode'])


        def _write_part_fingerprint_file( accumulated_lines ):
            part_dir = subdir_generator_function()
            part_path = os.path.join(part_dir, os.path.basename(input_files[0]))
            part_file = open(part_path, 'w')
            part_file.writelines( accumulated_lines )
            part_file.close()

        try:
            header_lines = []
            lines_accumulated = []
            fingerprint_counter = 0
            for line in open( input_files[0] ):
                if not line.strip():
                    continue
                if line.startswith('#'):
                    header_lines.append( line )
                else:
                    fingerprint_counter += 1
                    lines_accumulated.append( line )
                if fingerprint_counter != 0 and fingerprint_counter % chunk_size == 0:
                    _write_part_fingerprint_file( header_lines + lines_accumulated )
                    lines_accumulated = []
            if lines_accumulated:
                _write_part_fingerprint_file( header_lines + lines_accumulated )
        except Exception,  e:
            log.error('Unable to split files: %s' % str(e))
            raise
    split = classmethod(split)


    def merge(split_files, output_file):
        """
        Merging fps files requires merging the header manually.
        We take the header from the first file.
        """
        if len(split_files) == 1:
            #For one file only, use base class method (move/copy)
            return data.Text.merge(split_files, output_file)
        if not split_files:
            raise ValueError("No fps files given, %r, to merge into %s" \
                             % (split_files, output_file))
        out = open(output_file, "w")
        first = True
        for filename in split_files:
            with open(filename) as handle:
                for line in handle:
                    if line.startswith('#'):
                        if first:
                            out.write(line)
                    else:
                        # line is no header and not a comment, we assume the first header is written to out and we set 'first' to False
                        first = False
                        out.write(line)
        out.close()
    merge = staticmethod(merge)



class OBFS( Binary ):
    """OpenBabel Fastsearch format (fs)."""
    file_ext = 'fs'
    composite_type ='basic'
    allow_datatype_change = False

    MetadataElement( name="base_name", default='OpenBabel Fastsearch Index',
        readonly=True, visible=True, optional=True,)

    def __init__(self,**kwd):
        """
            A Fastsearch Index consists of a binary file with the fingerprints
            and a pointer the actual molecule file.
        """
        Binary.__init__(self, **kwd)
        self.add_composite_file('molecule.fs', is_binary = True,
            description = 'OpenBabel Fastsearch Index' )
        self.add_composite_file('molecule.sdf', optional=True,
            is_binary = False, description = 'Molecule File' )
        self.add_composite_file('molecule.smi', optional=True,
            is_binary = False, description = 'Molecule File' )
        self.add_composite_file('molecule.inchi', optional=True,
            is_binary = False, description = 'Molecule File' )
        self.add_composite_file('molecule.mol2', optional=True,
            is_binary = False, description = 'Molecule File' )
        self.add_composite_file('molecule.cml', optional=True,
            is_binary = False, description = 'Molecule File' )

    def set_peek( self, dataset, is_multi_byte=False ):
        """Set the peek and blurb text."""
        if not dataset.dataset.purged:
            dataset.peek  = "OpenBabel Fastsearch Index"
            dataset.blurb = "OpenBabel Fastsearch Index"
        else:
            dataset.peek = "file does not exist"
            dataset.blurb = "file purged from disk"

    def display_peek( self, dataset ):
        """Create HTML content, used for displaying peek."""
        try:
            return dataset.peek
        except:
            return "OpenBabel Fastsearch Index"

    def display_data(self, trans, data, preview=False, filename=None,
                     to_ext=None, size=None, offset=None, **kwd):
        """Apparently an old display method, but still gets called.

        This allows us to format the data shown in the central pane via the "eye" icon.
        """
        return "This is a OpenBabel Fastsearch format. You can speed up your similarity and substructure search with it."

    def get_mime(self):
        """Returns the mime type of the datatype (pretend it is text for peek)"""
        return 'text/plain'

    def merge(split_files, output_file, extra_merge_args):
        """Merging Fastsearch indices is not supported."""
        raise NotImplementedError("Merging Fastsearch indices is not supported.")

    def split( cls, input_datasets, subdir_generator_function, split_params):
        """Splitting Fastsearch indices is not supported."""
        if split_params is None:
            return None
        raise NotImplementedError("Splitting Fastsearch indices is not possible.")



class DRF( GenericMolFile ):
    file_ext = "drf"

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines('\"ligand id\"', dataset.file_name, invert = True)


class PHAR( GenericMolFile ):
    """
    Pharmacophore database format from silicos-it.
    """
    file_ext = "phar"
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            dataset.blurb = "pharmacophore"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'


class PDB( GenericMolFile ):
    """
    Protein Databank format.
    http://www.wwpdb.org/documentation/format33/v3.3.html
    """
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

class PDBQT(GenericMolFile):
    """
    Autodock Format
    The PDBQT format stores the atomic coordinates, partial charges and AutoDock atom types,
    for both the receptor and the ligand.
    http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
    """
    file_ext = "pdbqt"

    def sniff(self, filename):
        if count_special_lines('COMPND', filename) > 0:
            return True
        elif count_special_lines('REMARK', filename) > 0:
            return True
        elif count_special_lines('ATOM') > 0:
            return True
        else:
            return False

    def set_peek(self, dataset, is_multi_byte=False):
        if not dataset.dataset.purged:
            self.sniff(dataset.filename)
            dataset.blurb = "protein structure file used by autodock vina"
        else:
            dataset.peek = "file does not exist"
            dataset.blurb = "file purged from disk"


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
        currently its not possible to use external libraries in datatype definition files.
        Moreover it seems mpossible to inlcude OpenBabel as python library because OpenBabel
        is GPL licensed.
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



class CML( GenericXml ):
    """
    Chemical Markup Language
    http://cml.sourceforge.net/
    """
    file_ext = "cml"
    MetadataElement( name="number_of_molecules", default=0, desc="Number of molecules", readonly=True, visible=True, optional=True, no_value=0 )


    def set_meta( self, dataset, **kwd ):
        """
        Set the number of lines of data in dataset.
        """
        dataset.metadata.number_of_molecules = count_special_lines( '^\s*<molecule', dataset.file_name )

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
        Try to guess if the file is a CML file.
        TODO: add true positive test, need to submit a CML example

        >>> fname = get_test_fname( 'interval.interval' )
        >>> CML().sniff( fname )
        False
        """
        handle = open(filename)
        line = handle.readline()
        if line.strip() != '<?xml version="1.0"?>':
            handle.close()
            return False
        line = handle.readline()
        if line.strip().find('http://www.xml-cml.org/schema') == -1:
            handle.close()
            return False
        handle.close()
        return True


    def split( cls, input_datasets, subdir_generator_function, split_params):
        """
        Split the input files by molecule records.
        """
        if split_params is None:
            return None

        if len(input_datasets) > 1:
            raise Exception("CML-file splitting does not support multiple files")
        input_files = [ds.file_name for ds in input_datasets]

        chunk_size = None
        if split_params['split_mode'] == 'number_of_parts':
            raise Exception('Split mode "%s" is currently not implemented for CML-files.' % split_params['split_mode'])
        elif split_params['split_mode'] == 'to_size':
            chunk_size = int(split_params['split_size'])
        else:
            raise Exception('Unsupported split mode %s' % split_params['split_mode'])

        def _read_cml_records( filename ):
            lines = []
            with open(filename) as handle:
                for line in handle:
                    if line.lstrip().startswith('<?xml version="1.0"?>') or \
                        line.lstrip().startswith('<cml xmlns="http://www.xml-cml.org/schema') or \
                        line.lstrip().startswith('</cml>'):
                        continue
                    lines.append( line )
                    if line.lstrip().startswith('</molecule>'):
                        yield lines
                        lines = []

        header_lines = ['<?xml version="1.0"?>\n', '<cml xmlns="http://www.xml-cml.org/schema">\n']
        footer_line = ['</cml>\n']
        def _write_part_cml_file( accumulated_lines ):
            part_dir = subdir_generator_function()
            part_path = os.path.join(part_dir, os.path.basename(input_files[0]))
            part_file = open(part_path, 'w')
            part_file.writelines( header_lines )
            part_file.writelines( accumulated_lines )
            part_file.writelines( footer_line )
            part_file.close()

        try:
            cml_records = _read_cml_records( input_files[0] )
            cml_lines_accumulated = []
            for counter, cml_record in enumerate( cml_records, start = 1):
                cml_lines_accumulated.extend( cml_record )
                if counter % chunk_size == 0:
                    _write_part_cml_file( cml_lines_accumulated )
                    cml_lines_accumulated = []
            if cml_lines_accumulated:
                _write_part_cml_file( cml_lines_accumulated )
        except Exception,  e:
            log.error('Unable to split files: %s' % str(e))
            raise
    split = classmethod(split)


    def merge(split_files, output_file):
        """
        Merging CML files.
        """
        if len(split_files) == 1:
            #For one file only, use base class method (move/copy)
            return Text.merge(split_files, output_file)
        if not split_files:
            raise ValueError("Given no CML files, %r, to merge into %s" \
                             % (split_files, output_file))
        with open(output_file, "w") as out:
            for filename in split_files:
                with open( filename ) as handle:
                    header = handle.readline()
                    if not header:
                        raise ValueError("CML file %s was empty" % f)
                    if not header.lstrip().startswith('<?xml version="1.0"?>'):
                        out.write(header)
                        raise ValueError("%s is not a valid XML file!" % f)
                    line = handle.readline()
                    header += line
                    if not line.lstrip().startswith('<cml xmlns="http://www.xml-cml.org/schema'):
                        out.write(header)
                        raise ValueError("%s is not a CML file!" % f)
                    molecule_found = False
                    for line in handle.readlines():
                        # we found two required header lines, the next line should start with <molecule >
                        if line.lstrip().startswith('</cml>'):
                            continue
                        if line.lstrip().startswith('<molecule'):
                            molecule_found = True
                        if molecule_found:
                            out.write( line )
            out.write("</cml>\n")
    merge = staticmethod(merge)


if __name__ == '__main__':
    """
    TODO: We need to figure out, how to put example files under /lib/galaxy/datatypes/test/ from a toolshed, so that doctest can work properly.
    """
    inchi = get_test_fname('drugbank_drugs.inchi')
    smiles = get_test_fname('drugbank_drugs.smi')
    sdf = get_test_fname('drugbank_drugs.sdf')
    fps = get_test_fname('50_chemfp_fingerprints_FPS1.fps')
    pdb = get_test_fname('2zbz.pdb')
    cml = get_test_fname('/home/bag/Downloads/approved.cml')

    print 'CML test'
    print CML().sniff(cml), 'cml'
    print CML().sniff(inchi)
    print CML().sniff(pdb)
    CML().split()
    """
    print 'SMILES test'
    print SMILES().sniff(smiles), 'smi'
    print SMILES().sniff(inchi)
    print SMILES().sniff(pdb)
    """
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
