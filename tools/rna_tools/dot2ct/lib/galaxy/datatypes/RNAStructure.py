import logging
log = logging.getLogger(__name__)

from galaxy import util
import galaxy
import galaxy.model
import galaxy.datatypes
import galaxy.datatypes.data

from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.sequence import Sequence

import re

class DotBracket ( Sequence ):
    edam_format = "format_1457"
    file_ext = "dot-bracket"
    
    sequence_regexp = re.compile( "^[ACGTURYKMSWBDHVN]*" )
    structure_regexp = re.compile( "^[\(\)\.]*" )
    
    def set_meta( self, dataset, **kwd ):
        """
        Set the number of sequences and the number of data lines
        in dataset.
        """
        if self.max_optional_metadata_filesize >= 0 and dataset.get_size() > self.max_optional_metadata_filesize:
            dataset.metadata.data_lines = None
            dataset.metadata.sequences = None
            dataset.metadata.seconday_structures = None
            return
        
        data_lines = 0
        sequences = 0
        
        for line in file( dataset.file_name ):
            line = line.strip()
            data_lines += 1
            
            if line and line.startswith( '>' ):
                sequences += 1
        
        dataset.metadata.data_lines = data_lines
        dataset.metadata.sequences = sequences
    
    def sniff(self, filename):
        """
        The format is as follows, although it remains unclear whether
        the Dot-Bracket format may contain multiple sequences per file:
        
        >sequenceName1
        CCCaaaGGG
        (((...)))
        >sequenceName2
        GGGuuuCCC
        (((...)))
        """
        
        i = 0
        pairs = False
        
        with open( filename ) as handle:
            for line in handle:
                line = line.strip()
                
                state = i % 3
                
                if state == 0:#header line
                    if(line[0] != '>'):
                        return False
                elif state == 1:#sequence line
                    if not sequence_regexp.match(line.upper()):
                        return False
                    else:
                        sequence_size = len(line)
                elif state == 2:#dot-bracket structure line
                    if (sequence_size != len(line)) or (not structure_regexp.match(line)):
                        return False
                
                i += 1
