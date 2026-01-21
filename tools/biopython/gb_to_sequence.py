#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main(genbank_file, ofile, feature_name, transl_table = -1, to_protein = False, is_complete_cds = False, stop = True):

    handle = open(genbank_file, "rU")
    recs = _extract_features(SeqIO.parse(handle, 'genbank'), feature_name, transl_table, is_complete_cds, stop, to_protein)
    SeqIO.write(recs, ofile, 'fasta')
    handle.close()


def extract_features(gb_iter, feature_name, transl_table, is_complete_cds, stop, to_protein):

    feature_count = 0
    for gb_record in gb_iter:
        nucleotide = []
        protein = []
        for feature in gb_record.features:
            if feature.type == feature_name:
                # if the translation table is not set explicitly through the user, try to extract it from the genbank file
                # otherwise take the standard value 1
                if transl_table == -1 && qualifiers.key_exists('transl_table'):
                    transl_table = feature.qualifiers['transl_table'][0]
                else transl_table == -1:
                    transl_table = 1

                feature_count += 1
                gene = feature.extract(gb_record.seq)
                desc = feature.qualifiers.get('product', [''])[0]
                # try to extract an id and name for the feature
                if "locus_tag" in feature.qualifiers:
                    feature_id = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers:
                    feature_id = feature.qualifiers["gene"][0]
                else:
                    feature_id =  "%s_%s" % (feature_name, feature_count)
                if not to_protein:
                    seq = SeqRecord(gene, 
                        id = feature_id, 
                        name = feature_id, 
                        description = desc)
                else:
                    seq = SeqRecord( gene.translate( table = transl_table, to_stop=True, cds=is_complete_cds), 
                        id = feature_id, 
                        name = feature_id, 
                        description = desc)
                yield seq


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = "gbk_to_orf.py", description = 'Extract features from GenBank files.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--output', '-o', dest="ofile", 
        required = True,
        help='Path to the FASTA output file.')
    parser.add_argument('--input', '-i', dest="ifile", 
        required = True,
        help='Path to the GenBank file to read.')
    parser.add_argument('--feature', '-f', 
        dest="feature", 
        help="Feature to extract, for exampple 'gene' or 'cds'")
    parser.add_argument('--to-protein', '-p', 
        dest="to_protein", 
        action='store_true', 
        default=False, 
        help='Convert the output to amino acid sequences.')
    translation = parser.add_argument_group('Traslation parameters', 
        description='If you want to translate your extracted sequences to amino acids, you can specify some additional parameters.')
    translation.add_argument('--translation-table', '-t', 
        dest="translation_table", 
        type=int, 
        help='Number of the ncbi codon table, see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.')
    translation.add_argument('--complete-cds', '-c', 
        dest="is_complete_cds", 
        action='store_true', 
        default=False, 
        help='Indicate that the CDS is complete.')
    translation.add_argument('--stop', '-s', 
        dest="stop", 
        default = True,
        action='store_true', 
        help='Quit the translation of the CDS sequence on any stop codon.')
    options = parser.parse_args()

    main(options.ifile, options.ofile, options.feature, options.translation_table, options.to_protein, options.is_complete_cds, options.stop)



