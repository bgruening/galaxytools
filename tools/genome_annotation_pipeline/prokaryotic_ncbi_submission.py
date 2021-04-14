#!/usr/bin/env python
# -*- coding: UTF-8 -*-

__author__ = "Bjoern Gruening"
__version__ = "0.1"
__date__ = ""
__license__ = "GLP3+"

"""
These script converts a Blast XML Output to the NCBI Feature Table format,
required for submission of Serquence data.
http://www.ncbi.nlm.nih.gov/WebSub/index.cgi?tool=genbank
TODO: include t-RNA results
"""
import os
import shutil
import subprocess
import zipfile

from Bio import SeqIO
from Bio.Blast import NCBIXML
from glimmer_orf_to_seq import glimmer2sequence as g2seq
from RemoveVectorContamination import remove_vector_contamination
from utils import change_according_reviewer, change_glimmer3_prediction_output


def glimmer_prediction(
    fasta_dir, glimmer_trainingset="./glimmer3_strepto.icm", translation_table=11
):
    for root, dirs, files in os.walk(fasta_dir):
        for filename in files:
            if os.path.splitext(filename)[-1] in [".fasta", ".fa"]:
                path = os.path.join(root, filename)
                # com = "glimmerhmm %s %s -o %s -g -f" % (path, glimmer_trainingset, os.path.splitext( path )[0] + '.glimmerhmm')
                com = (
                    "glimmer3 -o50 -g90 -t30 -l --trans_table %s %s %s %s 2> /dev/null"
                    % (
                        translation_table,
                        path,
                        glimmer_trainingset,
                        os.path.splitext(path)[0] + ".glimmer",
                    )
                )
                # com = "tigr-glimmer glimmer3 -o0 -g90 -t30 -l --trans_table %s %s %s %s 2> /dev/null" % (translation_table, path, glimmer_trainingset, os.path.splitext( path )[0] + '.glimmer')
                subprocess.call(com, shell=True, stdout=subprocess.PIPE)
                # Glimmer3 will generate two file one is called _base_.predict the other _base_.detail
                # detail - will be deleted, and predict will be renamed to _base_ without any ending
                change_glimmer3_prediction_output(
                    os.path.splitext(path)[0] + ".glimmer.predict",
                    os.path.splitext(path)[0] + ".glimmer",
                )
                os.remove(os.path.splitext(path)[0] + ".glimmer.predict")
                os.remove(os.path.splitext(path)[0] + ".glimmer.detail")


def glimmer2sequence(fasta_dir, translation_table=11):
    for root, dirs, files in os.walk(fasta_dir):
        for filename in files:
            if os.path.splitext(filename)[-1] in [".glimmer"]:
                path = os.path.join(root, filename)
                base_path_name = os.path.splitext(path)[0]
                g2seq(
                    base_path_name + ".fasta",
                    path,
                    base_path_name + "_predicted.fasta",
                    to_protein=False,
                    translation_table=translation_table,
                )


def get_glimmer_mapping(path):
    mapping = dict()
    with open(path) as glimmer_handle:
        for line in glimmer_handle:
            if not line.startswith(">"):
                orf_name, start, end = line.split()[:3]
                mapping[orf_name] = (start, end)
    return mapping


def run_blast(data_dir, blastdb, threads=8):
    for root, dirs, files in os.walk(data_dir):
        for filename in files:
            if filename.endswith("_predicted.fasta"):
                predicted_proteins = os.path.join(root, filename)
                blastxml_file = predicted_proteins.replace(
                    "_predicted.fasta", ".blastxml"
                )
                if open(predicted_proteins).read().strip() == "":
                    continue
                # com = "blastp -query %s -db %s -task blastp -evalue 0.001 -out %s -outfmt 5 -num_threads 3 -seg no" % (predicted_proteins, blastdb, blastxml_file)
                com = (
                    "blastx -query %s -db %s -evalue 0.001 -out %s -outfmt 5 -num_threads %s"
                    % (predicted_proteins, blastdb, blastxml_file, threads)
                )
                subprocess.call(com, shell=True, stdout=subprocess.PIPE)


def run(data_dir, feature_table_path, locus_tag, min_coverage, min_ident):
    feature_table = open(feature_table_path, "w+")
    gene_counter = 0
    annotation_count_with_putative_function = 0

    # To get a sorted Feature table, go through the whole directory save it in a list and sort it
    filepaths = {}
    for root, dirs, files in os.walk(data_dir):
        for filename in files:
            if filename.endswith(".blastxml"):
                blastxml_file = os.path.join(root, filename)
                seq_number = filename.split("Seq")[-1].split(".blastxml")[0]
                filepaths[int(seq_number)] = blastxml_file
    sorted_iids = sorted(filepaths.keys())
    for iids in sorted_iids:
        blastxml_file = filepaths[iids]
        if open(blastxml_file).read().strip() == "":
            continue
        glimmer_mapping = get_glimmer_mapping(
            blastxml_file.replace(".blastxml", ".glimmer")
        )
        gene_counter, annotation_count_with_putative_function = parse_blastxml(
            blastxml_file,
            glimmer_mapping,
            feature_table,
            annotation_count_with_putative_function,
            gene_counter,
            locus_tag,
            min_coverage,
            min_ident,
        )
    print(
        "From %s total genes %s are annotated."
        % (
            gene_counter,
            annotation_count_with_putative_function,
        )
    )


def parse_blastxml(
    input_path,
    glimmer_mapping,
    feature_table,
    annotation_count_with_putative_function,
    gene_counter,
    locus_tag,
    min_coverage,
    min_ident,
):
    # extract the sequence number
    seq_number = input_path.split("Seq")[-1].split(".blastxml")[0]
    feature_table.write(">Feature Seq%s\n" % seq_number)

    with open(input_path) as blast_handle:
        for entry in NCBIXML.parse(blast_handle):
            if entry.application == "BLASTX":
                query_length = entry.query_length
                if type(query_length) == type(1.7):
                    print("Query length is not a multiple of three")
                    break
                query_id = entry.query.split()[0]
                query_info = glimmer_mapping[query_id]
                query_start = int(query_info[0])
                query_end = int(query_info[1])
            else:
                break

            gene_counter += 1
            """
            if not entry.alignments:
                feature_table.write('%i\t%i\tgene\n' % (query_start, query_end))
                feature_table.write('\t\t\tlocus_tag\t%s%04d\n' % (locus_tag, gene_counter))
                feature_table.write('%i\t%i\tCDS\n' % (query_start, query_end))
                feature_table.write('\t\t\tproduct\thypothetical protein\n')
                feature_table.write('\t\t\tprotein_id\tgnl|PBUF|%s%04d\n' % (locus_tag, gene_counter))
                feature_table.write('\t\t\tnote\tpredicted with glimmer3\n')
                break
            """
            hsp_has_annotation = False
            feature_table_text = dict()

            for alignment in entry.alignments:
                for hsp in alignment.hsps:
                    nident = hsp.identities
                    ident = 100 * float(nident) / float(hsp.align_length)
                    coverage = False
                    if (
                        int(hsp.query_end) - int(hsp.query_start) + 1
                        >= min_coverage * query_length
                    ):
                        coverage = True

                    # only annotate hits with an identity over 50% and a coverage over 50%
                    if ident > min_ident and coverage:
                        feature_table_text[hsp.bits] = ""
                        hsp_has_annotation = True

                        """
                        Hit_def changed: It now looks like: 
                        'RecName: Full=Erythronolide synthase, modules 3 and 4; Short=PKS; AltName: Full=6-deoxyerythronolide B synthase II; AltName: Full=DEBS 2; AltName: Full=ORF 2'
                        """
                        print(alignment.hit_def)
                        accession = alignment.hit_def.encode("utf8")
                        accession = filter(
                            lambda token: token.startswith("RecName:"),
                            map(str.strip, accession.split(";")),
                        )[0].split("Full=")[-1]

                        assert (
                            change_according_reviewer(
                                "Pimelyl-[acyl-carrier protein] methyl ester esterase",
                                note_line=False,
                            )
                            == "Pimelyl-[acyl-carrier protein] methyl ester esterase"
                        )
                        assert (
                            change_according_reviewer(
                                "putative D-malate dehydrogenase [decarboxylating] [gnl|PBUF|STVIR_0046:1-352] [gnl|PBUF|STVIR_0046: raw, aa len= 352]",
                                note_line=False,
                            )
                            == "D-malate dehydrogenase"
                        )
                        accession = change_according_reviewer(
                            accession, note_line=False
                        )

                        feature_table_text[hsp.bits] += "%i\t%i\tgene\n" % (
                            query_start,
                            query_end,
                        )
                        feature_table_text[hsp.bits] += "\t\t\tlocus_tag\t%s%04d\n" % (
                            locus_tag,
                            gene_counter,
                        )
                        feature_table_text[hsp.bits] += "%i\t%i\tCDS\n" % (
                            query_start,
                            query_end,
                        )

                        if (
                            accession.startswith("hypothetical protein")
                            or accession.startswith("predicted protein")
                            or accession == ""
                            or accession == "protein"
                        ):
                            feature_table_text[
                                hsp.bits
                            ] += "\t\t\tproduct\thypothetical protein\n"
                        else:
                            feature_table_text[
                                hsp.bits
                            ] += "\t\t\tproduct\tputative %s\n" % (accession)

                        feature_table_text[
                            hsp.bits
                        ] += "\t\t\tprotein_id\tgnl|PBUF|%s%04d\n" % (
                            locus_tag,
                            gene_counter,
                        )

                        if str(hsp.expect).find("e") != -1:
                            """ Der evalue ist eine lange Zahl und muss gekuertzt werden. Z.B. 4.787347812347e-124"""
                            evalue_first, evalue_last = str(hsp.expect).split("e")
                            evalue = (
                                str(round(float(evalue_first), 1)) + "e" + evalue_last
                            )
                        else:
                            evalue = round(hsp.expect, 1)

                        """
                        hit_def = change_according_reviewer(alignment.hit_def, note_line = True)
                        if hit_def.split('|')[:-1] != []:
                            hit_def = hit_def.split('|')[-1].split()[0]
                        else:
                            hit_def = accession
                        """

                        """"
                        try:
                            protein_accession_gb = hit_def.split('gb|')[1].split('|')[0] #try to extract the genbank accession number >gi|302432474|gb|EFL04290.1|; -> EFL04290.1
                            inference = "similar to AA sequence:INSD: %s" % protein_accession_gb
                            feature_table_text[ hsp.bits ] += '\t\t\tinference\t%s\n' % (inference)

                            protein_accession_ref = hit_def.split('ref|')[1].split('|')[0] #try to extract the genbank accession number >gi|302432474|gb|EFL04290.1|; -> EFL04290.1
                            inference = "similar to AA sequence:RefSeq: %s" % protein_accession_ref
                            feature_table_text[ hsp.bits ] += '\t\t\tinference\t%s\n' % (inference)
                        except Exception:
                            pass
                        """
                        inference = """ab initio prediction:Glimmer:3"""
                        feature_table_text[hsp.bits] += "\t\t\tinference\t%s\n" % (
                            inference
                        )

                        note = """similar to UniProtKB/Swiss-Prot Entry: %(hit_accession)s""" % {
                            "gene_counter": gene_counter,
                            "accession": accession,
                            "alignment_hit_def": accession,
                            "hit_accession": alignment.accession,
                            "len": query_length,
                            "evalue": evalue,
                            "bit_score": round(hsp.bits, 2),
                            "locus_tag": locus_tag,
                        }
                        feature_table_text[hsp.bits] += "\t\t\tnote\t%s\n" % (note)

            if not hsp_has_annotation:
                """
                If hsp has no annotation with the specified identity and coverage, insert a hypothetical protein
                """
                feature_table.write("%i\t%i\tgene\n" % (query_start, query_end))
                feature_table.write(
                    "\t\t\tlocus_tag\t%s%04d\n" % (locus_tag, gene_counter)
                )
                feature_table.write("%i\t%i\tCDS\n" % (query_start, query_end))
                feature_table.write("\t\t\tproduct\thypothetical protein\n")
                feature_table.write(
                    "\t\t\tprotein_id\tgnl|PBUF|%s%04d\n" % (locus_tag, gene_counter)
                )
                feature_table.write("\t\t\tnote\tab initio prediction:Glimmer3\n")
            else:
                bitscores = feature_table_text.keys()
                bitscores.sort(reverse=True)
                feature_table.write(feature_table_text[bitscores[0]])
                if (
                    feature_table_text[bitscores[0]].find(
                        "\t\t\tproduct\thypothetical protein\n"
                    )
                    == -1
                ):
                    annotation_count_with_putative_function += 1

    return (gene_counter, annotation_count_with_putative_function)


def splitFastaFile(infile, informat, outdir):
    for record in SeqIO.parse(open(infile), informat):
        iid = record.id
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        f_out = os.path.join(outdir, iid + ".fasta")
        SeqIO.write([record], open(f_out, "w"), "fasta")


def zipper(dir, zip_file):
    zip = zipfile.ZipFile(zip_file, "w", compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


if __name__ == "__main__":
    import argparse
    import shutil
    import tempfile

    import Scaffold2Fasta
    from change_fasta_header import change_fasta_header

    parser = argparse.ArgumentParser(description="Creates a NCBI FeatureTable.")

    parser.add_argument(
        "-u",
        "--univec",
        dest="univec_path",
        default="/media/mydata/univec/tue6071_stupl",
        help="Path to the UniVex BLAST Database. It is used for contamination checking.",
    )

    parser.add_argument(
        "--blastdb",
        dest="blastdb_path",
        required=True,
        help="Path to the SwissProt BLAST Database. It is used for annotating proteins.",
    )

    parser.add_argument(
        "--num-threads",
        dest="num_threads",
        default=8,
        type=int,
        help="Number of threads to use for similarity searching.",
    )

    parser.add_argument(
        "--glimmer-trainingset",
        dest="glimmer_trainingset",
        default="./glimmer3_strepto.icm",
        help="Path to th glimmer trainingset.",
    )

    parser.add_argument(
        "-d",
        "--datadir",
        dest="data_dir",
        default=tempfile.mkdtemp(),
        help="the data directoy to store all temporary files",
    )
    parser.add_argument(
        "--scaffold",
        dest="scaffold",
        required=False,
        help="FASTA scaffold file, can contain gaps in the sequence",
    )
    parser.add_argument(
        "-f",
        "--feature-table",
        dest="feature_table",
        default=False,
        help="output path of the NCBI Feature-Table",
    )
    parser.add_argument(
        "--cleaned-sequence",
        dest="cleaned_sequence",
        default=False,
        help="""The output path where the processed sequence,
                            without gaps larger than 1 and no contamination,
                            should be stored""",
    )
    parser.add_argument(
        "--agp-file",
        dest="agp_file",
        default=False,
        help="""The output path where the AGP file
                            should be stored""",
    )
    parser.add_argument(
        "--genbank-file",
        dest="gbf_file",
        default=False,
        help="""The output path where the GenBank file
                            should be stored""",
    )
    parser.add_argument(
        "--sequin-file",
        dest="sqn_file",
        default=False,
        help="""The output path where the Seqin file
                            should be stored""",
    )
    parser.add_argument(
        "--validation-file",
        dest="val_file",
        default=False,
        help="""The output path where the validation file
                            should be stored""",
    )
    parser.add_argument(
        "--sbt-file",
        dest="sbt_file",
        default=False,
        help="""The input path from the author file. It can be created on the ncbi page.\n
                      http://www.ncbi.nlm.nih.gov/WebSub/template.cgi""",
    )
    parser.add_argument(
        "--compressed-archive",
        dest="compressed_results",
        default=False,
        help="""Path to an archive of all results. It should contain all relevant files for a NCBI submission.""",
    )
    parser.add_argument(
        "--discrepancy_report",
        dest="discrepancy_report",
        default=False,
        help="""Path to the Discrepancy Report Output File.""",
    )
    parser.add_argument(
        "--translation-table",
        dest="translation_table",
        type=int,
        default=11,
        help="""The ncbi translation table number. See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.""",
    )

    parser.add_argument(
        "--minimal-coverage",
        dest="min_coverage",
        default=0.5,
        type=float,
        help="""Minimal coverage of a BLAST hit to include the annotation in the output.""",
    )
    parser.add_argument(
        "--minimal-identity",
        dest="min_ident",
        default=0.5,
        type=float,
        help="""Minimal identity of a BLAST hit to include the annotation in the output.""",
    )

    parser.add_argument(
        "--sequence-description",
        dest="seq_description",
        default="",
        help="""The sequence description will be inserted in each FASTA header and will be 
                      included in the genbank file. The NCBI reviever suggested something like [organism=Glarea lozoyensis 74030] [strain=74030]""",
    )

    """
        If you want to choose a different locus_tag prefix, it must be at least 3 characters.
        See http://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf.
        You can check to see if a prefix is available here:  http://www.ncbi.nlm.nih.gov/genomes/lltp.cgi
    """
    parser.add_argument(
        "--locus-tag",
        dest="locus_tag",
        default=False,
        required=True,
        help="""The locus tag each scaffold and gene gets. It should be at least 3 characters.""",
    )

    options = parser.parse_args()

    if not os.path.exists(options.data_dir):
        os.mkdir(options.data_dir)

    # options.data_dir = '/tmp/tmpomltKT'

    print(options.data_dir)

    nocontamination_path = os.path.join(options.data_dir, "nocontamination.fsa")

    if not options.agp_file:
        # its not set throught the user
        options.agp_file = os.path.join(options.data_dir, "scaffold.agp")

    if not options.feature_table:
        # its not set throught the user
        options.feature_table = os.path.join(options.data_dir, "feature_table.tbl")

    # If no cleand sequence is provided, clean it in two steps.
    # The galaxy wrapper provides a cleaned sequence, because the two steps running extern, we we don't need that steps.
    if not options.cleaned_sequence:
        # its not set throught the user
        options.cleaned_sequence = os.path.join(
            options.data_dir, "ncbi_submission_sequence.fsa"
        )
        # change_title_path = os.path.join( options.data_dir, 'changed_title.fsa' )
        # change_fasta_header(options.scaffold, change_title_path, 'id_only')
        remove_vector_contamination(
            options.scaffold, options.univec_path, nocontamination_path
        )
        # convert the scaffold file in smaller contigs, remove the N-runs
        Scaffold2Fasta.run(
            nocontamination_path,
            options.cleaned_sequence,
            options.agp_file,
            options.locus_tag,
            options.seq_description,
        )

    # Split a multiple FASTA file in separate FASTA files
    splitFastaFile(options.cleaned_sequence, "fasta", options.data_dir)

    glimmer_prediction(
        options.data_dir, options.glimmer_trainingset, options.translation_table
    )
    glimmer2sequence(options.data_dir, options.translation_table)

    run_blast(
        options.data_dir, blastdb=options.blastdb_path, threads=options.num_threads
    )
    run(
        options.data_dir,
        options.feature_table,
        options.locus_tag + "_",
        options.min_coverage,
        options.min_ident,
    )

    # ./BlastXML_to_NCBIFeatureTable.py -d ./split --scaffold Glarea-losoyensis_scaffold.fasta -f submission/Youssar_Glarea-lozoyensis.tbl
    # use tbl2asn to create a full submission out of the feature table, the sequence and the 'submission template'
    # to create a submission template you can use: http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
    # tbl2asn -t ./submission/Youssar_Glarea-lozoyensis.sbt -p ./submission/ -V vb -a s -M n -Z discrep.report

    """
    Copy the generated file in a new temporary folder to run the tbl2asn tool against it.
    tbl2asn will produce a genbank file with a validation of the generated data.
    """

    tmp_path = tempfile.mkdtemp()

    shutil.copyfile(
        options.feature_table, os.path.join(tmp_path, options.locus_tag + ".tbl")
    )
    shutil.copyfile(
        options.cleaned_sequence, os.path.join(tmp_path, options.locus_tag + ".fsa")
    )
    shutil.copyfile(
        options.agp_file, os.path.join(tmp_path, options.locus_tag + ".agp")
    )

    if options.sbt_file:
        com = (
            "tbl2asn -t %(sbt_file)s -p %(tmp_path)s -V vb -a s -M n -Z %(tmp_path)s/discrep.report"
            % {"sbt_file": options.sbt_file, "tmp_path": tmp_path}
        )
        shutil.copyfile(
            options.sbt_file, os.path.join(tmp_path, options.locus_tag + ".sbt")
        )
    else:
        com = (
            "tbl2asn -p %(tmp_path)s -V vb -a s -M n -Z %(tmp_path)s/discrep.report"
            % {"tmp_path": tmp_path}
        )
    subprocess.call(com, shell=True, stdout=subprocess.PIPE)

    if options.val_file:
        shutil.copyfile(
            os.path.join(tmp_path, options.locus_tag + ".val"), options.val_file
        )
    if options.gbf_file:
        shutil.copyfile(
            os.path.join(tmp_path, options.locus_tag + ".gbf"), options.gbf_file
        )
    if options.sqn_file:
        shutil.copyfile(
            os.path.join(tmp_path, options.locus_tag + ".sqn"), options.sqn_file
        )
    if options.discrepancy_report:
        shutil.copyfile(
            os.path.join(tmp_path, "discrep.report"), options.discrepancy_report
        )

    if options.compressed_results:
        zipper(tmp_path, options.compressed_results)

    # clean temp data files
    shutil.rmtree(options.data_dir)
    shutil.rmtree(tmp_path)

    """
    ./NCBI_annotation_prokaryotes.py -d ./temp/ --scaffold viridochromogens_part.fasta --locus-tag VIR --sequence-description "[organism=Streptomyces viridochromogenes Tue494] [strain=Tue494] [gcode=11]"
    """
