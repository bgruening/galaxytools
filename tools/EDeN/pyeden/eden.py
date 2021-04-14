"""
    EDeN filetypes
"""

from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes import data


class Gspan(Tabular):
    """Class describing an gSpan file"""

    file_ext = "gspan"

    def set_peek(self, dataset, is_multi_byte=False):
        if not dataset.dataset.purged:
            dataset.peek = "gSpan"
            dataset.blurb = data.nice_size(dataset.get_size())
        else:
            dataset.peek = "file does not exist"
            dataset.blurb = "file purged from disk"

    def display_peek(self, dataset):
        try:
            return dataset.peek
        except Exception:
            return "Tabular gSpan file (%s)" % (data.nice_size(dataset.get_size()))


class SparseVector(Tabular):
    """Class describing an SparseVector file"""

    file_ext = "sparsevector"

    def set_peek(self, dataset, is_multi_byte=False):
        if not dataset.dataset.purged:
            dataset.peek = "SparseVector"
            dataset.blurb = data.nice_size(dataset.get_size())
        else:
            dataset.peek = "file does not exist"
            dataset.blurb = "file purged from disk"

    def display_peek(self, dataset):
        try:
            return dataset.peek
        except Exception:
            return "Tabular SparseVector file (%s)" % (
                data.nice_size(dataset.get_size())
            )
