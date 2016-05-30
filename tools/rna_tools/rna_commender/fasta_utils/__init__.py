"""Util functions for FASTA format."""

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


def import_fasta(fasta_file):
    """Import a fasta file as a dictionary."""
    dic = {}
    f = open(fasta_file)
    fasta = f.read().strip()
    f.close()
    for a in fasta.split('>'):
        k = a.split('\n')[0]
        v = ''.join(a.split('\n')[1:])
        if k != '':
            dic[k] = v
    return dic


def export_fasta(dic):
    """Export a dictionary."""
    fasta = ""
    for (k, v) in dic.iteritems():
        fasta += ">%s\n%s\n" % (k, v)
    return fasta


def seq_names(fasta_file):
    """Get sequence names from fasta file."""
    names = []
    f = open(fasta_file)
    fasta = f.read()
    f.close()
    for a in fasta.split('>'):
        names.append(a.split('\n')[0])
    return [a for a in names if a != '']


def stockholm2fasta(stockholm):
    """Convert alignment in stockholm format to fasta format."""
    fasta = ""
    for line in stockholm.split("\n"):
        # comment line
        if line[0] == "#":
            continue
        # termination line
        elif line == "//":
            return fasta
        # alignment line
        else:
            name, align = line.split()
            seq = align.replace(".", "")
            fasta += ">%s\n%s\n" % (name, seq)
