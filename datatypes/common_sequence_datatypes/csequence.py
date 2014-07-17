"""
Classes for all common sequence formats
"""

from galaxy.datatypes.data import get_file_peek
from galaxy.datatypes import data
from galaxy.datatypes.metadata import MetadataElement
import logging
log = logging.getLogger(__name__)


class GenBank( data.Text ):
    """
        abstract class for most of the sequence files
    """
    file_ext = "genbank"

    MetadataElement( name="number_of_sequences", default=0, desc="Number of sequences", readonly=True, visible=True, optional=True, no_value=0 )
    MetadataElement( name="name_of_sequences", default=[], desc="Name of sequences", readonly=True, visible=True, optional=True, no_value=0 )

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            if (dataset.metadata.number_of_sequences == 1):
                dataset.blurb = "1 sequence [%s]" % dataset.metadata.name_of_sequences[0]
            else:
                seq_list = ', '.join(dataset.metadata.name_of_sequences[0:10])
                if len(dataset.metadata.name_of_sequences) > 10:
                    seq_list += ', ...'
                dataset.blurb = "%s sequences [%s]" % ( dataset.metadata.number_of_sequences,
                                                       seq_list )
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def get_mime(self):
        return 'text/plain'

    def sniff( self, filename ):
        header = open(filename).read(5)
        return header == 'LOCUS'

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of sequences in dataset.
        """
        (count, seq_names) = self._count_genbank_sequences( dataset.file_name )
        dataset.metadata.number_of_sequences = count
        dataset.metadata.name_of_sequences = seq_names
    def split( cls, input_datasets, subdir_generator_function, split_params):
        """
        Split the input files by sequence records.
        """
        pass
    split = classmethod(split)

    def _count_genbank_sequences( self, filename ):
        """
        This is not a perfect definition, but should suffice for general usage. It fails to detect any
        errors that would result in parsing errors like incomplete files.
        """
        # Specification for the genbank file format can be found in
        # ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
        # in section 3.4.4 LOCUS Format
        count = 0
        seq_names = []
        with open( filename ) as gbk:
            for line in gbk:
                if line[0:5] == 'LOCUS':
                    seq_names.append(line[12:28].strip())
                    count += 1
        return (count, seq_names)
