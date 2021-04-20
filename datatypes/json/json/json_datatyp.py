# -*- coding: utf-8 -*-

from galaxy.datatypes.data import Text
from galaxy.datatypes.data import get_file_peek
import json
import os

class Json( Text ):
    file_ext = "json"

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            dataset.blurb = "JavaScript Object Notation (JSON)"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disc'

    def sniff( self, filename ):
        """
            Try to load the string with the json module. If successful it's a json file.
        """
        try:
            json.load( open(filename) )
            return True
        except Exception:
            return False

    def display_peek( self, dataset ):
        try:
            return dataset.peek
        except Exception:
            return "JSON file (%s)" % ( data.nice_size( dataset.get_size() ) )



