"""
Classes for all common sequence formats
"""

from galaxy.datatypes.data import get_file_peek
from galaxy.datatypes import data
from galaxy.datatypes.xml import GenericXml
from galaxy.datatypes.metadata import MetadataElement

import os
import logging

log = logging.getLogger(__name__)

def count_genbank_sequences( filename ):
    """
    """
    return 0

class GenBank( data.Text ):
    """
        abstract class for most of the molecule files
    """
    file_ext = "genbank"

    MetadataElement( name="number_of_sequences", default=0, desc="Number of sequences", readonly=True, visible=True, optional=True, no_value=0 )

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            if (dataset.metadata.number_of_molecules == 1):
                dataset.blurb = "1 sequence"
            else:
                dataset.blurb = "%s sequences" % dataset.metadata.number_of_molecules
            dataset.peek = data.get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def get_mime(self):
        return 'text/plain'

    def sniff( self, filename ):
        pass

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of molecules in dataset.
        """
        dataset.metadata.number_of_sequences = count_genbank_sequences( dataset.file_name )

    def split( cls, input_datasets, subdir_generator_function, split_params):
        """
        Split the input files by sequence records.
        """
        pass
    split = classmethod(split)
