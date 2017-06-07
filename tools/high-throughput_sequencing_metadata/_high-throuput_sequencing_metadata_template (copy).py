#!/usr/bin/env python
# -*- coding: UTF-8 -*-

__author__ = "Bjoern Gruening"
__version__ = "0.1"
__date__ = ""
__license__ = 'GLP3+'

"""
    These script converts a Blast XML Output to the NCBI Feature Table format,
    required for submission of Serquence data.
    http://www.ncbi.nlm.nih.gov/WebSub/index.cgi?tool=genbank
    TODO: include t-RNA results
"""
import os, sys
import csv
from Bio.Blast import NCBIXML
from Bio import SeqIO
import re
from collections import defaultdict
from glimmerhmm_gff_to_sequence import main as gff2seq
import subprocess
from RemoveVectorContamination import *
import zipfile
from utils import change_according_reviewer


def augustus_prediction(fasta_dir, trainingset = "aspergillus_nidulans"):
    for root, dirs, files in os.walk( fasta_dir ):
        for filename in files:
            if os.path.splitext(filename)[-1] in ['.fasta','.fa']:
                path = os.path.join(root, filename)
                com = """augustus \
                            --gff3=on \
                            --strand=both \
                            --stopCodonExcludedFromCDS=false \
                            --genemodel=complete \
                            --noInFrameStop=true \
                            --species='%s' \
                            --outfile=%s \
                            %s \
                            """ % (trainingset, os.path.splitext( path )[0] + '.gff3', path,)
                subprocess.call( com, shell=True, stdout=subprocess.PIPE )



def gff2sequence(fasta_dir, traslation_table):
    for root, dirs, files in os.walk( fasta_dir ):
        for filename in files:
            if os.path.splitext(filename)[-1] in ['.gff3']:
                path = os.path.join(root, filename)
                base_path_name = os.path.splitext( path )[0]
                gff2seq(path, base_path_name + '.fasta', base_path_name + '_predicted.fasta', False, ncbi_traslation_table = traslation_table)


class Segment():
    def __init__(self, start, stop, strand, seq_type, frame, note, iid):
        if strand == '-':
            # we need to swap start and stop
            self.start = stop
            self.stop = start
        else:
            self.start = start
            self.stop = stop
        self.strand = strand
        self.seq_type = seq_type
        self.note = note
        self.iid = iid
        self.frame = frame

    def __repr__(self):
        return "Segment (%s:%s): [start: %s, end: %s, strand: %s, frame: %s, note: %s]" % (
            self.seq_type, self.iid, self.start, self.stop,
            self.strand, self.frame, self.note)


class AugustusEntry():
    def __init__(self, name, seq_id = None):
        self.name = name
        self.seq_id = seq_id
        self.index = 0
        #self.has_initial_exon = False
        #self.has_final_exon = False
        self.seq_type = 'gene' # can be gene or CDS
        self.exons = []
        self.mRNA = None

    def add_segment(self, start, stop, strand, seq_type, frame, note = ''):
        if seq_type == 'CDS':
            self.seq_type = seq_type
            self.add_exon( int(start), int(stop), strand, self.name,
                           int(frame), note)
        elif seq_type == 'gene':
            self. seq_type = seq_type
            self.add_mRNA( int(start), int(stop), strand, self.name, note )
        else:
            print seq_type
            raise ValueError
        """
        if note == 'final-exon':
            self.has_final_exon = True
        elif note == 'initial-exon':
            self.has_initial_exon = True
        elif note == 'single-exon':
            # if it is only one exon it has per definition and end and a start
            self.has_final_exon = True
            self.has_initial_exon = True
        """
    def add_mRNA(self, mRNA_start, mRNA_stop, strand, iid, note = ''):
        #mRNA has no frame, default to 0
        self.mRNA = Segment( mRNA_start, mRNA_stop, strand,
                        self.seq_type, '0', note, iid )

    def add_exon(self, exon_start, exon_stop, strand, iid, frame, note = ''):
        self.exon_start = exon_start
        self.exon_stop = exon_stop
        self.strand = strand # +/-
        s = Segment( exon_start, exon_stop, strand, self.seq_type, frame, note, iid)
        if strand == '-':
            self.exons.insert( 0, s )
        else:
            self.exons.append( s )
    
    def __getsegments__(self):
        if self.mRNA == None:
            segments = self.exons
        else:
            segments = [self.mRNA] + self.exons
        return segments

    def __repr__(self):
        segments = self.__getsegments__()
        return '\n'.join( str(s) for s in segments )

    def __getitem__(self, key):
        segments = self.__getsegments__()
        return segments.__getitem__( key )

    def __setitem__(self, key):
        segments = self.__getsegments__()
        return segments.__setitem__( key )

    def __getslice__(self, i, j):
        segments = self.__getsegments__()
        return segments.__getslice__( i, j )

    def __iter__(self):
        return self

    def next(self):
        segments = self.__getsegments__()
        if self.index == len(segments):
            raise StopIteration
        self.index += 1
        return segments[self.index]



def get_augustus_mapping(path):
    """
        ToDo: change the datastructure, some nice class
    """
    mapping = dict()

    with open(path) as augustus_handle:
        for line in augustus_handle:
            line = line.strip()
            if not line.startswith('#') and line:
                columns = line.split('\t')
                seq_id = columns[0]
                seq_type = columns[2]
                start = columns[3]
                stop = columns[4]
                strand = columns[6]
                frame = columns[7]
                metadata = columns[-1]
                """
                Seq4	AUGUSTUS	gene	213143	215159	.	-	.	ID=g1582
                Seq4	AUGUSTUS	transcript	213143	215159	0.31	-	.	ID=g1582.t1;Parent=g1582
                """
                iid = metadata.split('ID=')[-1].split(';')[0].split('.')[0]

                if seq_type in ['CDS', 'gene']:
                    if mapping.has_key( iid ):
                        mapping[iid].add_segment( start, stop, strand, seq_type, frame)
                    else:
                        ge = AugustusEntry( iid, seq_id )
                        ge.add_segment( start, stop, strand, seq_type, frame )
                        mapping[iid] = ge
    return mapping


def get_organism(accession):
    return accession.split('OS=')[-1].split('GN=')[0]


def run_blast(data_dir, blastdb = '/media/mydata/uniprot/blastdb/uniprot_sprot.blastdb', threads = 8):
    for root, dirs, files in os.walk( data_dir ):
        for filename in files:
            if filename.endswith('_predicted.fasta'):
                predicted_proteins = os.path.join(root, filename)
                blastxml_file = predicted_proteins.replace('_predicted.fasta', '.blastxml')
                if open(predicted_proteins).read().strip() == '':
                    continue
                #com = "blastp -query %s -db %s -task blastp -evalue 0.001 -out %s -outfmt 5 -num_threads 10 -seg no -max_target_seqs 1" % (predicted_proteins, blastdb, blastxml_file)
                com = "blastx -query %s -db %s -evalue 0.001 -out %s -outfmt 5 -num_threads %s" % (predicted_proteins, blastdb, blastxml_file, threads)
                subprocess.call( com, shell=True, stdout=subprocess.PIPE )


def run(data_dir, feature_table_path, locus_tag, min_coverage, min_ident):
    feature_table = open(feature_table_path, 'w+')
    gene_counter = 0
    annotation_count_with_putative_function = 0

    # To get a sorted Feature table, go through the whole directory save it in a list and sort it
    filepaths = {}
    for root, dirs, files in os.walk( data_dir ):
        for filename in files:
            if filename.endswith('.blastxml'):
                blastxml_file = os.path.join(root, filename)
                if filename.find( 'Seq' ) != -1:
                    seq_number = filename.split('Seq')[-1].split('.blastxml')[0]
                    filepaths[int(seq_number)] = blastxml_file
                else:
                    # if the name is not from the kind of 'Seq34' use the original name
                    filepaths[blastxml_file] = blastxml_file

    sorted_iids = sorted(filepaths.keys())
    for iids in sorted_iids:
        blastxml_file = filepaths[iids]
        if open(blastxml_file).read().strip() == '':
            continue
        augustus_mapping = get_augustus_mapping( blastxml_file.replace('.blastxml', '.gff3') )
        gene_counter, annotation_count_with_putative_function = parse_blastxml(blastxml_file, augustus_mapping, feature_table, annotation_count_with_putative_function, gene_counter, locus_tag, min_coverage, min_ident)

    print 'From %s total genes %s are annotated.' % (gene_counter, annotation_count_with_putative_function)

def iter_islast(iterable):
    """
        iter_islast(iterable) -> generates (item, islast) pairs

        Generates pairs where the first element is an item from the iterable
        source and the second element is a boolean flag indicating if it is the
        last item in the sequence.
    """
    it = iter(iterable)
    prev = it.next()
    for item in it:
        yield prev, False
        prev = item
    yield prev, True


def check_short_introns(cds):
    """ 
        Check for short introns!
        'We prefer not to annotate every possible open reading frame.  Instead, only annotate features when you think they represent real genes.
        Do you believe these are real genes with frameshifts?  If so, we prefer not to include translations for CDS features when the translation is known to be incorrect.
        Please remove the CDS and mRNA features and annotate this with a single gene feature across the entire span.  Include the pseudo qualifier and a note on the gene with a brief description.  
        For example:
        1       200     gene
                gene    phoA
                gene_desc       alkaline phosphatase
                locus_tag     POR_0001
                pseudo
                note    nonfunctional due to frameshift'
    """
    old_end = 0
    old_start = 0
    for region in cds:
        end = int(region.stop)
        start = int(region.start)

        """
            [2] is the strand
            If the Intron-Gap is smaller than 18 and its not the first exon.
            First exon check, because if the start is 0 (or <18), that means
            the CDS begins in first position of the sequence,
            than it would be a short intron (0-0) < 18.
        """
        if abs(start - old_end) < 18 and old_start != 0:
            return True
        old_end = end
        old_start = start
    return False


def parse_blastxml(input_path, augustus_mapping, feature_table, annotation_count_with_putative_function, gene_counter, locus_tag, min_coverage, min_ident):

    #print input_path, augustus_mapping, feature_table, annotation_count_with_putative_function, gene_counter, locus_tag
    # extract the sequence name
    seq_name = os.path.splitext(os.path.split(input_path)[-1])[0]
    #input_path.split('Seq')[-1].split('.blastxml')[0]
    #print seq_name
    feature_table.write('>Feature %s\n' % seq_name)
    #locus_tag = 'M7I_'
    with open(input_path) as blast_handle:
        for entry in NCBIXML.parse(blast_handle):
            if entry.application == "BLASTX":
                query_length = entry.query_length / 3
                if type(query_length) == type(1.7):
                    print "Query length is not a multiple of three"
                    break
                query_id = entry.query.split()[0]
                query_info = augustus_mapping[ query_id ]

                assert query_info.mRNA.seq_type == 'gene'
                gene_start = query_info.mRNA.start
                gene_end = query_info.mRNA.stop
                cds = query_info.exons
                mRNA = query_info.mRNA
            else:
                break

            gene_counter += 1
            hsp_has_annotation = False
            feature_table_text = dict()
            for alignment in entry.alignments:
                for hsp in alignment.hsps:
                    nident = hsp.identities
                    ident = (100*float(nident)/float(hsp.align_length))
                    """
                        Coverage: 'c8-c7+1 >= 0.5*c23'
                    """
                    coverage = False
                    if int(hsp.query_end) - int(hsp.query_start) + 1 >= min_coverage * query_length:
                        coverage = True
                    # only annotate hits with an identity over 50% and a coverage over 50%
                    if ident > min_ident and coverage:
                        feature_table_text[ hsp.bits ] = ""
                        hsp_has_annotation = True
                        """
                        Hit_def changed: It now looks like: 
                        'RecName: Full=Erythronolide synthase, modules 3 and 4; Short=PKS; AltName: Full=6-deoxyerythronolide B synthase II; AltName: Full=DEBS 2; AltName: Full=ORF 2'
                        """
                        print alignment.hit_def
                        accession = alignment.hit_def.encode('utf8')
                        accession = filter(lambda token: token.startswith('RecName:'), map(str.strip, accession.split(';')))[0].split('Full=')[-1]
                        accession = change_according_reviewer(accession, note_line = False)

                        feature_table_text[ hsp.bits ] += '%i\t%i\tgene\n' % (gene_start, gene_end)
                        feature_table_text[ hsp.bits ] += '\t\t\tlocus_tag\t%s%04d\n' % (locus_tag, gene_counter)

                        short_intron = check_short_introns(cds)
                        if short_intron:
                            feature_table_text[ hsp.bits ] += '\t\t\tpseudo\n'
                            feature_table_text[ hsp.bits ] += '\t\t\tnote\tnonfunctional; similar to %s\n' % accession
                            continue
                        """
                            Write the CDS section for the 'annotation' case and save a string for the mRNA section.
                        """
                        mRNA_annotation = ''
                        mRNA_annotation += '%i\t%i\tmRNA\n' % (cds[0].start, cds[0].stop)

                        feature_table_text[ hsp.bits ] += '%i\t%i\tCDS\n' % (cds[0].start, cds[0].stop)
                        for region in cds[1:]:
                            feature_table_text[ hsp.bits ] += '%i\t%i\n' % (region.start, region.stop)
                            mRNA_annotation += '%i\t%i\n' % (region.start, region.stop)

                        if accession.startswith('hypothetical protein') or \
                                accession.startswith('predicted protein') or \
                                accession == '' or accession == 'protein':
                            feature_table_text[ hsp.bits ] += '\t\t\tproduct\thypothetical protein\n'
                            mRNA_annotation += '\t\t\tproduct\thypothetical protein\n'
                        else:
                            feature_table_text[ hsp.bits ] += '\t\t\tproduct\tputative %s\n' % (accession)
                            mRNA_annotation += '\t\t\tproduct\tputative %s\n' % (accession)

                        feature_table_text[ hsp.bits ] += '\t\t\tprotein_id\tgnl|PBUF|%s%04d\n' % (locus_tag, gene_counter)
                        mRNA_annotation += '\t\t\tprotein_id\tgnl|PBUF|%s%04d\n' % (locus_tag, gene_counter)
                        feature_table_text[ hsp.bits ] += '\t\t\ttranscript_id\tgnl|PBUF|%smrna%04d\n' % (locus_tag, gene_counter)
                        mRNA_annotation += '\t\t\ttranscript_id\tgnl|PBUF|%smrna%04d\n' % (locus_tag, gene_counter)

                        # Write mRNA section
                        feature_table_text[ hsp.bits ] += mRNA_annotation

                        if str(hsp.expect).find('e') != -1:
                            """ Der evalue ist eine lange Zahl und muss gekuertzt werden. Z.B. 4.787347812347e-124"""
                            evalue_first, evalue_last = str(hsp.expect).split('e')
                            evalue = str(round(float(evalue_first), 1)) + 'e' + evalue_last
                        else:
                            evalue = round(hsp.expect, 1)
                        """
                        hit_def = change_according_reviewer(alignment.hit_def, note_line = True)
                        if hit_def.split('|')[:-1] != []:
                            hit_def = hit_def.split('|')[-1].split()[0]
                        else:
                            hit_def = accession
                        """
                        hit_def = accession
                        """
                        try:
                            protein_accession_gb = hit_def.split('gb|')[1].split('|')[0] #try to extract the genbank accession number >gi|302432474|gb|EFL04290.1|; -> EFL04290.1
                            inference = "similar to AA sequence:INSD: %s" % protein_accession_gb
                            feature_table_text[ hsp.bits ] += '\t\t\tinference\t%s\n' % (inference)

                            protein_accession_ref = hit_def.split('ref|')[1].split('|')[0] #try to extract the genbank accession number >gi|302432474|gb|EFL04290.1|; -> EFL04290.1
                            inference = "similar to AA sequence:RefSeq: %s" % protein_accession_ref
                            feature_table_text[ hsp.bits ] += '\t\t\tinference\t%s\n' % (inference)
                        except:
                            pass
                        """
                        inference = """ab initio prediction:Augustus:2.5.5"""
                        feature_table_text[ hsp.bits ] += '\t\t\tinference\t%s\n' % (inference)

                        note = """similar to UniProtKB/Swiss-Prot Entry: %(hit_accession)s""" % {'gene_counter': gene_counter, 
                                'accession':accession,
                                'alignment_hit_def': hit_def,
                                'hit_accession': alignment.accession,
                                'len': query_length,
                                'evalue': evalue,
                                'bit_score': round(hsp.bits, 2),
                                'locus_tag': locus_tag,
                                }

                        feature_table_text[ hsp.bits ] += '\t\t\tnote\t%s\n' % (note)


                #for region in cds[1:]:
                #    mRNA_annotation += '%i\t%i\n' % (region.start, region.stop)
                #    feature_table.write('%i\t%i\n' % (region.start, region.stop))


            if hsp_has_annotation == False:

                """
                    If hsp has no annotation, insert a hypothetical protein
                """
                feature_table.write('%i\t%i\tgene\n' % (gene_start, gene_end))
                feature_table.write('\t\t\tlocus_tag\t%s%04d\n' % (locus_tag, gene_counter))
                assert cds[0].seq_type == 'CDS'
                short_intron = check_short_introns(cds)
                if short_intron:
                    feature_table.write('\t\t\tpseudo\n')
                    feature_table.write('\t\t\tnote\tnonfunctional\n')
                """
                    Write the CDS section for the 'no-annotation' case.
                """
                feature_table.write('%i\t%i\tCDS\n' % (cds[0].start, cds[0].stop))

                for region in cds[1:]:
                    feature_table.write('%i\t%i\n' % (region.start, region.stop))

                feature_table.write('\t\t\tproduct\thypothetical protein\n')
                feature_table.write('\t\t\tprotein_id\tgnl|PBUF|%s%04d\n' % (locus_tag, gene_counter))
                feature_table.write('\t\t\ttranscript_id\tgnl|PBUF|%smrna%04d\n' % (locus_tag, gene_counter))
                feature_table.write('\t\t\tnote\tpredicted with Augustus 2.5.5\n')

                """
                    Write the mRNA section for the 'no-annotation' case.
                """
                feature_table.write('%i\t%i\tmRNA\n' % (cds[0].start, cds[0].stop))
                for region in cds[1:]:
                    feature_table.write('%i\t%i\n' % (region.start, region.stop))
                feature_table.write('\t\t\tproduct\thypothetical protein\n')
                feature_table.write('\t\t\tprotein_id\tgnl|PBUF|%s%04d\n' % (locus_tag, gene_counter))
                feature_table.write('\t\t\ttranscript_id\tgnl|PBUF|%smrna%04d\n' % (locus_tag, gene_counter))

            else:
                bitscores = feature_table_text.keys()
                bitscores.sort(reverse=True)
                feature_table.write(feature_table_text[ bitscores[0] ])
                if feature_table_text[ bitscores[0] ].find('\t\t\tproduct\thypothetical protein\n') == -1:
                    annotation_count_with_putative_function += 1

    return (gene_counter, annotation_count_with_putative_function)


def splitFastaFile(infile, informat, outdir):
    for record in SeqIO.parse(open(infile), informat):
        iid = record.id
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        f_out = os.path.join(outdir,iid+'.fasta')
        SeqIO.write([record],open(f_out,'w'),"fasta")


def zipper(dir, zip_file):
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


if __name__ == '__main__':

    import argparse
    import tempfile
    import Scaffold2Fasta
    import shutil
    parser = argparse.ArgumentParser(description='Creates a NCBI FeatureTable.')

    parser.add_argument("-u", "--univec", dest="univec_path",
                      default="/media/mydata/univec/tue6071_stupl",
                      help="Path to the UniVex BLAST Database. It is used for contamination checking.")

    parser.add_argument("--num-threads", dest="num_threads",
                      default=8, type=int,
                      help="Number of threads to use for similarity searching.")

    parser.add_argument("--blastdb", dest="blastdb_path",
                      default="/media/mydata/uniprot/blastdb/uniprot_sprot.blastdb",
                      help="Path to the SwissProt BLAST Database. It is used for annotating proteins.")

    parser.add_argument("--trainingset", dest="trainingset",
                      default="aspergillus_nidulans",
                      help="Name of the trainingset, that should be used.")

    parser.add_argument("-d", "--datadir", dest="data_dir",
                      default = tempfile.mkdtemp(),
                      help="the data directoy to store all temporary files")
    parser.add_argument("--scaffold", dest="scaffold",
                      required=False,
                      help="FASTA scaffold file, can contain gaps in the sequence")
    parser.add_argument("-f", "--feature-table",
                      dest="feature_table", default=False,
                      help="output path of the NCBI Feature-Table")
    parser.add_argument("--cleaned-sequence", dest="cleaned_sequence",
                      default = False,
                      help="""The output path where the processed sequence,
                            without gaps larger than 1 and no contamination,
                            should be stored""")
    parser.add_argument("--agp-file", dest="agp_file",
                      default = False,
                      help="""The output path where the AGP file
                            should be stored""")
    parser.add_argument("--genbank-file", dest="gbf_file",
                      default = False,
                      help="""The output path where the GenBank file
                            should be stored""")
    parser.add_argument("--sequin-file", dest="sqn_file",
                      default = False,
                      help="""The output path where the Seqin file
                            should be stored""")
    parser.add_argument("--validation-file", dest="val_file",
                      default = False,
                      help="""The output path where the validation file
                            should be stored""")
    parser.add_argument("--sbt-file", dest="sbt_file",
                      default = False,
                      help="""The input path from the author file. It can be created on the ncbi page.\n
                      http://www.ncbi.nlm.nih.gov/WebSub/template.cgi""")

    parser.add_argument("--sequence-description", dest="seq_description",
                      default = "",
                      help="""The sequence description will be inserted in each FASTA header and will be 
                      included in the genbank file. The NCBI reviever suggested something like [organism=Glarea lozoyensis 74030] [strain=74030]""")
    parser.add_argument("--compressed-archive", dest="compressed_results",
                      default = False,
                      help="""Path to an archive of all results. It should contain all relevant files for a NCBI submission.""")
    parser.add_argument("--translation-table", dest="translation_table", type=int,
                      default = 1,
                      help="""The ncbi translation table number. See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.""")

    parser.add_argument("--minimal-coverage", dest="min_coverage",
                      default = 0.5, type=float,
                      help="""Minimal coverage of a BLAST hit to include the annotation in the output.""")
    parser.add_argument("--minimal-identity", dest="min_ident",
                      default = 0.5, type=float,
                      help="""Minimal identity of a BLAST hit to include the annotation in the output.""")

    parser.add_argument("--discrepancy_report", dest="discrepancy_report",
                      default = False,
                      help="""Path to the Discrepancy Report Output File.""")


    """
        If you want to choose a different locus_tag prefix, it must be at least 3 characters.
        See http://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf.
        You can check to see if a prefix is available here:  http://www.ncbi.nlm.nih.gov/genomes/lltp.cgi
    """
    parser.add_argument("--locus-tag", dest="locus_tag",
                      default = False, required=True,
                      help="""The locus tag each scaffold and gene gets. It should be at least 3 characters.""")

    options = parser.parse_args()

    if not os.path.exists( options.data_dir ):
        os.mkdir( options.data_dir )

    print options.data_dir

    nocontamination_path = os.path.join( options.data_dir, 'nocontamination.fsa' )
    if options.cleaned_sequence == False:
        # its not set throught the user
        options.cleaned_sequence = os.path.join( options.data_dir, 'ncbi_submission_sequence.fsa' )

    if options.agp_file == False:
        # its not set throught the user
        options.agp_file = os.path.join( options.data_dir, 'scaffold.agp' )

    if options.feature_table == False:
        # its not set throught the user
        options.feature_table = os.path.join( options.data_dir, 'feature_table.tbl' )

    # If no cleand sequence is provided, clean it in two steps.
    # The galaxy wrapper provides a cleaned sequence, because the two steps running extern, we we don't need that steps.
    if options.cleaned_sequence == False:
        remove_vector_contamination( options.scaffold, options.univec_path, nocontamination_path )
        # convert the scaffold file in smaller contigs, remove the N-runs
        Scaffold2Fasta.run( nocontamination_path, options.cleaned_sequence, options.agp_file, options.locus_tag, options.seq_description)

    # Split a multiple FASTA file in separate FASTA files
    splitFastaFile( options.cleaned_sequence, 'fasta', options.data_dir )

    augustus_prediction( options.data_dir, options.trainingset)
    gff2sequence( options.data_dir, options.translation_table)

    run_blast( options.data_dir, blastdb = options.blastdb_path, threads = options.num_threads)
    run( options.data_dir, options.feature_table, options.locus_tag + '_', options.min_coverage, options.min_ident )

    #./BlastXML_to_NCBIFeatureTable.py -d ./split --scaffold Glarea-losoyensis_scaffold.fasta -f submission/Youssar_Glarea-lozoyensis.tbl
    # use tbl2asn to create a full submission out of the feature table, the sequence and the 'submission template'
    # to create a submission template you can use: http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
    # tbl2asn -t ./submission/Youssar_Glarea-lozoyensis.sbt -p ./submission/ -V vb -a s  -M n -Z discrep.report

    """
        Copy the generated file in a new temporary folder to run the tbl2asn tool against it.
        tbl2asn will produce a genbank file with a validation of the generated data.
    """

    tmp_path = tempfile.mkdtemp()
    shutil.copyfile(options.feature_table, os.path.join( tmp_path, options.locus_tag + '.tbl'))
    shutil.copyfile(options.cleaned_sequence, os.path.join( tmp_path, options.locus_tag + '.fsa'))
    shutil.copyfile(options.agp_file, os.path.join( tmp_path, options.locus_tag + '.agp'))

    if options.sbt_file:
        com = "tbl2asn -t %(sbt_file)s -p %(tmp_path)s -V vb -a s -M n -Z %(tmp_path)s/discrep.report" % {'sbt_file': options.sbt_file, 'tmp_path': tmp_path}
        shutil.copyfile(options.sbt_file, os.path.join( tmp_path, options.locus_tag + '.sbt'))
    else:
        com = "tbl2asn -p %(tmp_path)s -V vb -a s -M n -Z %(tmp_path)s/discrep.report" % {'tmp_path': tmp_path}
    subprocess.call( com, shell=True, stdout=subprocess.PIPE )

    if options.val_file:
        shutil.copyfile(os.path.join( tmp_path, options.locus_tag + '.val'), options.val_file)
    if options.gbf_file:
        shutil.copyfile(os.path.join( tmp_path, options.locus_tag + '.gbf'), options.gbf_file)
    if options.sqn_file:
        shutil.copyfile(os.path.join( tmp_path, options.locus_tag + '.sqn'), options.sqn_file)
    if options.discrepancy_report:
        shutil.copyfile(os.path.join( tmp_path, 'discrep.report'), options.discrepancy_report)

    if options.compressed_results:
        zipper(tmp_path, options.compressed_results)

    # clean temp data files
    shutil.rmtree(options.data_dir)
    shutil.rmtree(tmp_path)


