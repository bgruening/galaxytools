#!/usr/bin/env python

"""
Extract sequences from a GFF3 file.

This works for examlpe with the GlimmerHMM GFF3 output format:

#gff-version 3
#sequence-region Glarea 1 38755942
Glarea	GlimmerHMM	mRNA	1017	1644	.	-	.	ID=Glarea.path1.gene1;Name=Glarea.path1.gene1
Glarea	GlimmerHMM	CDS	1017	1126	.	-	2	ID=Glarea.cds1.1;Parent=Glarea.path1.gene1;Name=Glarea.path1.gene1;Note=final-exon
Glarea	GlimmerHMM	CDS	1258	1317	.	-	2	ID=Glarea.cds1.2;Parent=Glarea.path1.gene1;Name=Glarea.path1.gene1;Note=internal-exon
Glarea	GlimmerHMM	CDS	1638	1644	.	-	0	ID=Glarea.cds1.3;Parent=Glarea.path1.gene1;Name=Glarea.path1.gene1;Note=initial-exon
Glarea	GlimmerHMM	mRNA	2755	5399	.	-	.	ID=Glarea.path1.gene2;Name=Glarea.path1.gene2
Glarea	GlimmerHMM	CDS	2755	4826	.	-	2	ID=Glarea.cds2.1;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=final-exon
Glarea	GlimmerHMM	CDS	5071	5145	.	-	2	ID=Glarea.cds2.2;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=internal-exon
Glarea	GlimmerHMM	CDS	5202	5214	.	-	0	ID=Glarea.cds2.3;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=internal-exon
Glarea	GlimmerHMM	CDS	5337	5399	.	-	0	ID=Glarea.cds2.4;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=initial-exon
Glarea	GlimmerHMM	mRNA	6580	7583	.	+	.	ID=Glarea.path1.gene3;Name=Glarea.path1.gene3

"""

import argparse
import os
import shutil
import tempfile
import zipfile

from BCBio import GFF
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm
from svg_stack import (AlignCenter, Document, VBoxLayout, convert_to_pixels,
                       get_unit_attr)


def basic_parsing(gffFile):
    """ GFF3 parse, extract information """

    with open(gffFile) as in_handle:
        for rec in GFF.parse(in_handle):
            # iterate features
            for feature in rec.features:
                # iterate sub features
                sub_features_temp = []
                for sub_feature in feature.sub_features:
                    sub_features_temp.append(sub_feature)
                yield SeqFeature(
                    type=feature.type,
                    location=feature.location,
                    strand=feature.strand,
                    qualifiers=feature.qualifiers,
                    sub_features=sub_features_temp,
                )


def visualize_gff3(gffFile, oformat, circular, outfile, compressed_archive=False):
    """ """
    index = 0
    outname = os.path.splitext(os.path.basename(gffFile))[0]
    temp_dir = tempfile.mkdtemp()

    if oformat == "svg" and not compressed_archive:
        # Use svg_stack to bundle the different plots together in one file
        margin = "100px"
        margin_gx = convert_to_pixels(*get_unit_attr(margin))
        layout = VBoxLayout()
        doc = Document()
        svg_handle = open(outfile, mode="w")

    # iterate the features in gff3 file
    for feature in basic_parsing(gffFile):
        index += 1
        temp_plot = os.path.join(temp_dir, "%s_%d" % (outname, index))

        seq_starts = []
        seq_ends = []
        gd_diagram = GenomeDiagram.Diagram()
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
        gd_feature_set = gd_track_for_features.new_set()

        # iterate sub-features
        for sub_feature in feature.sub_features:
            if len(gd_feature_set) % 2 == 0:
                color = colors.blue
            else:
                color = colors.lightblue
            gd_feature_set.add_feature(
                sub_feature,
                sigil="ARROW",
                arrowshaft_height=0.2,
                arrowhead_length=0.25,
                label_position="middle",
                color=color,
                label=True,
                label_size=14,
                label_angle=20,
                name=sub_feature.qualifiers.get("Name", ["no name"])[0],
            )
            seq_ends.append(sub_feature.location.end.position)
            seq_starts.append(sub_feature.location.start.position)

        # draw linear or circular
        if not circular:
            gd_diagram.draw(
                format="linear",
                pagesize="A4",
                fragments=4,
                start=min(seq_starts),
                end=max(seq_ends),
            )
        else:
            gd_diagram.draw(
                format="circular",
                circular=True,
                pagesize=(20 * cm, 20 * cm),
                start=min(seq_starts),
                end=max(seq_ends),
            )

        # specify the output format
        if oformat == "pdf":
            gd_diagram.write(temp_plot + ".pdf", "PDF")
        elif oformat == "eps":
            gd_diagram.write(temp_plot + ".eps", "EPS")
        elif oformat == "svg":
            gd_diagram.write(temp_plot + ".svg", "SVG")
            if not compressed_archive:
                layout.addSVG(temp_plot + ".svg", alignment=AlignCenter)

    if compressed_archive:
        outfile = zipper(temp_dir, outfile)
    else:
        # merge files and clean files
        if oformat == "pdf":
            os.system(
                "gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=%s -dBATCH %s/*pdf"
                % (outfile, temp_dir)
            )
        elif oformat == "eps":
            os.system(
                "gs -dNOPAUSE -sDEVICE=epswrite -dSAFER -sOutputFile=%s.eps -dBATCH  %s/*eps"
                % (outfile, temp_dir)
            )
        elif oformat == "svg":
            layout.setSpacing(margin_gx)
            doc.setLayout(layout)
            doc.save(svg_handle)

    # cairosvg.svg2png(url=svg_filename, write_to=png_filename, dpi=72)
    # cairosvg.svg2pdf(url=outfile, write_to=outfile+'.png')
    # remove temp dir with all single images
    shutil.rmtree(temp_dir)


def zipper(dir, zip_file):
    """ """
    zip = zipfile.ZipFile(zip_file, "w", compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for filename in files:
            fullpath = os.path.join(root, filename)
            archive_name = os.path.join(archive_root, filename)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract sequence information from gff file."
    )
    parser.add_argument("-o", "--outfile", help="path to the image")
    parser.add_argument(
        "-f",
        "--oformat",
        dest="oformat",
        default="svg",
        choices=["pdf", "eps", "svg"],
        help="Output format.",
    )
    parser.add_argument(
        "--circular", action="store_true", default=False, help="Circulat plots."
    )
    parser.add_argument(
        "--compressed-archive",
        dest="compressed_archive",
        action="store_true",
        default=False,
        help="create one single file with all genes, rather than multiple files",
    )
    parser.add_argument(
        "-g", "--gff", dest="gff_path", required=True, help="GFF input file."
    )

    options = parser.parse_args()

    visualize_gff3(
        options.gff_path,
        options.oformat,
        options.circular,
        options.outfile,
        options.compressed_archive,
    )
