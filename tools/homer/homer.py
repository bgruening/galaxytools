"""
HOMER special datatypes
"""
import os

from galaxy.datatypes.data import Data, Text, get_file_peek
from galaxy.datatypes.images import Html
from galaxy.datatypes.metadata import MetadataElement


class TagDirectory(Html):
    """Base class for HOMER's Tag Directory datatype."""

    file_ext = "homer_tagdir"
    composite_type = "auto_primary_file"
    allow_datatype_change = False

    def __init__(self, **kwd):
        Html.__init__(self, **kwd)
        # self.add_composite_file('tagInfo.txt', description = 'basic configuration information', mimetype = 'text/html') # Contains basic configuration information
        self.add_composite_file(
            "tagLengthDistribution.txt",
            description="histogram of read lengths used for alignment",
            mimetype="text/html",
        )  # File contains a histogram of read lengths used for alignment.
        self.add_composite_file(
            "tagCountDistribution.txt",
            description="histogram of clonal read depth, showing the number of reads per unique position",
            mimetype="text/html",
        )  # File contains a histogram of clonal read depth, showing the number of reads per unique position.
        self.add_composite_file(
            "tagAutocorrelation.txt",
            description="distribution of distances between adjacent reads in the genome",
            mimetype="text/html",
        )  # The autocorrelation routine creates a distribution of distances between adjacent reads in the genome.
        self.add_composite_file(
            "tagFreq.txt",
            description="nucleotide and dinucleotide frequencies as a function of distance from the 5' end of all reads",
            mimetype="text/html",
            optional=True,
        )  # Calculates the nucleotide and dinucleotide frequencies as a function of distance from the 5' end of all reads.
        self.add_composite_file(
            "tagFreqUniq.txt",
            description="nucleotide and dinucleotide frequencies as a function of distance from the 5' end of all reads (counted only once)",
            mimetype="text/html",
            optional=True,
        )  # Same as tagFreq.txt, however individual genomic positions are only counted once.
        self.add_composite_file(
            "tagGCcontent.txt",
            description="Distribution of fragment GC%-content",
            mimetype="text/html",
            optional=True,
        )  # Distribution of fragment GC%-content.
        self.add_composite_file(
            "genomeGCcontent.txt",
            description="Distribution of fragment GC%-content at each location in the genome",
            mimetype="text/html",
            optional=True,
        )  # Distribution of fragment GC%-content at each location in the genome.

    def regenerate_primary_file(self, dataset):
        """
        regenerate the index file after metadata generation
        """
        rval = ["<html><head><title>HOMER database files</title></head>"]
        rval.append("<body>")
        rval.append("<p/>CuffDiff Outputs:<p/><ul>")
        for fname in os.listdir(dataset.files_path):
            sfname = os.path.split(fname)[-1]
            rval.append('<li><a href="%s" type="text/html">%s</a>' % (sfname, sfname))
        rval.append("</ul></body></html>")
        f = file(dataset.file_name, "w")
        f.write("%s\n" % "\n".join(rval))
        f.close()
        if not dataset.info:
            dataset.info = "HOMER datatype object"
        if not dataset.blurb:
            dataset.blurb = "Composite file - HOMER"
        return True

    def generate_primary_file(self, dataset=None):
        rval = ["<html><head><title>HOMER database files</title></head><ul>"]
        for composite_name, composite_file in self.get_composite_files(
            dataset=dataset
        ).iteritems():
            opt_text = ""
            if composite_file.optional:
                opt_text = " (optional)"
            rval.append(
                '<li><a href="%s">%s</a>%s' % (composite_name, composite_name, opt_text)
            )
        rval.append("</ul></html>")
        return "\n".join(rval)

    def set_meta(self, dataset, **kwd):
        Html.set_meta(self, dataset, **kwd)
        self.regenerate_primary_file(dataset)

    def display_data(
        self,
        trans,
        data,
        preview=False,
        filename=None,
        to_ext=None,
        size=None,
        offset=None,
        **kwd
    ):
        """Apparently an old display method, but still gets called.

        This allows us to format the data shown in the central pane via the "eye" icon.
        """
        return "This is a HOMER database."

    def set_peek(self, dataset, is_multi_byte=False):
        """Set the peek and blurb text."""
        if not dataset.dataset.purged:
            dataset.peek = "HOMER database (multiple files)"
            dataset.blurb = "HOMER database (multiple files)"
        else:
            dataset.peek = "file does not exist"
            dataset.blurb = "file purged from disk"

    def display_peek(self, dataset):
        """Create HTML content, used for displaying peek."""
        try:
            return dataset.peek
        except Exception:
            return "HOMER database (multiple files)"

    def get_mime(self):
        """Returns the mime type of the datatype (pretend it is text for peek)"""
        return "text/plain"

    def merge(split_files, output_file):
        """Merge HOMER databases (not implemented)."""
        raise NotImplementedError("Merging HOMER databases is not supported")

    def split(cls, input_datasets, subdir_generator_function, split_params):
        """Split a HOMER database (not implemented)."""
        if split_params is None:
            return None
        raise NotImplementedError("Can't split HOMER databases")
