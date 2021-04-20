# -*- coding: utf-8 -*-

#from galaxy.datatypes.json_datatyp import Json as JsonClass
from galaxy.datatypes.data import Text
from galaxy.datatypes.data import get_file_peek
from galaxy import util
import subprocess
import tempfile
import logging
import json
import os

log = logging.getLogger(__name__)

#class Ipynb( JsonClass ):
class Ipynb( Text ):
    file_ext = "ipynb"

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek( dataset.file_name, is_multi_byte=is_multi_byte )
            dataset.blurb = "IPython Notebook"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disc'

    def sniff( self, filename ):
        """
            Try to load the string with the json module. If successful it's a json file.
        """
        try:
            ipynb = json.load( open(filename) )
            if ipynb.get('nbformat', False) != False and ipynb.get('metadata', False):
                return True
            else:
                return False
        except Exception:
            return False

    def display_data(self, trans, dataset, preview=False, filename=None, to_ext=None, chunk=None, **kwd):
        preview = util.string_as_bool( preview )
        if chunk:
            return self.get_chunk(trans, dataset, chunk)
        elif to_ext or not preview:
            return self._serve_raw(trans, dataset, to_ext)
        else:
            ofile_handle = tempfile.NamedTemporaryFile(delete=False)
            ofilename = ofile_handle.name
            ofile_handle.close()
            try:
                cmd = 'ipython nbconvert --to html --template basic %s --output %s' % (dataset.file_name, ofilename)
                subprocess.call(cmd, shell=True)
                ofilename = '%s.html' % ofilename
            except Exception:
                ofilename = dataset.file_name
                log.exception( 'Command "%s" failed. Could not convert the IPython Notebook to HTML, defaulting to plain text.' % cmd )
            return open( ofilename )

    def set_meta( self, dataset, **kwd ):
        """
        Set the number of models in dataset.
        """
        pass


