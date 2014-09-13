#!/usr/bin/env python

"""
    antiSMASH 2.0 output plugin to write all cluster proteins to a file (*_genecluster_proteins.fa)
"""
import logging
import textwrap
from os import path
from antismash import utils

name = "genecluster_proteins"
short_description = "Ouptut gene clusters as FASTA sequences"
# Output plugins are sorted by priority, lower numbers get run first
priority = 9

def write(seq_records, options):
    """Write all cluster proteins to a file

    Args:
        seq_records (iterable): An iterable containing Bio.SeqRecords
        options (argparse.Namespace): The options passed to the program
    """
    basename = seq_records[0].id
    output_name = path.join(options.outputfoldername, "%s_genecluster_proteins.fa" % basename)
    logging.debug("Writing seq_records to %r" % output_name)

    with open(output_name, 'w+') as handle:
        for seq_record in seq_records:
            clusters = utils.get_cluster_features(seq_record)
            for cluster in clusters:
                clustertype = utils.get_cluster_type(cluster)
                clusternr = utils.get_cluster_number(cluster)
                for feature in utils.get_cluster_cds_features(cluster, seq_record):
                    qual = feature.qualifiers
                    fasta_header = '>%s:%s %s #%s - %s\n' % (qual['locus_tag'][0], qual['protein_id'][0], clustertype, clusternr, qual['product'][0])
                    handle.write( fasta_header )
                    handle.write( '%s\n' % '\n'.join( textwrap.wrap(qual['translation'][0], 60) ) )
