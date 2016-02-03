#!/bin/bash
#echo "The script you are running has basename `basename $0`, dirname `dirname $0`"
#echo "The present working directory is `pwd`"

SCRIPT_DIR=`dirname $0`
report_script="${SCRIPT_DIR}/rcas.Rmd"
css="${SCRIPT_DIR}/base/custom.css"
outdir=`pwd`
header="${SCRIPT_DIR}/base/header.html"

function usage()
{
    	echo "RCAS v.0.1"
	echo "Authors: AA, DY, BU, RW"
    	echo "./generate_report.sh"
    	echo "	-h --help					print this help message"
	echo "	-o --output_filename=<out_file.html>		output filename for the report.html"
	echo "	-a --annot=<annot_file> 			annotation file as output by parse_anot.py"
	echo "	-p --peaks=<peaks.bed>				input BED file which was used as input to the pipeline"
	echo "	-g --gff3=<input.gff3>				genome annotations in gff3 format"
	echo "	-b --go_bp=<go_bp.tsv>				file containing GO term analysis results for Biological Processes"
	echo "	-m --go_mf=<go_mf.tsv>				file containing GO term analysis results for Molecular Functions"
	echo "	-c --go_cc=<go_cc.tsv>				file containing GO term analysis results for Cellular Compartments"
	echo "	-s --msigdb=<msigdb_file.tsv>			file containing MSIGDB gene set enrichment analysis results as output from rcas.msigdb.R"
	echo "	-e --meme_out=</path/to/meme_out_dir>		path to the folder that contains the output of meme results (e.g. /path/to/meme_out/)"
	echo "	-t --motif_annot=<motif_anot_file.tsv>		file containing the output of top_motifs.py scripts, which contains the annotations related to the detected motifs in the experiment" 
	echo "	-i --coverage_profile_option=<option>		option to run or not run coverage profile calculations: choose 'NOT_RUN' in order to turn it off. choose 'RUN' to keep it on."
	echo ""
}

if [ "$1" == "" ]; then
	usage
	exit 1
fi

while [ "$1" != "" ]; do
	PARAM=`echo $1 | awk -F= '{print $1}'`
	VALUE=`echo $1 | awk -F= '{print $2}'`
    	case $PARAM in
        	-h | --help)
	        	usage
        	    	exit
            		;;
		-o | --output_filename)
	    		output_filename=$VALUE
	    		;;
		-a | --annot)
			annotation_file=$VALUE
			;;
		-p | --peaks)
			peaks_file=$VALUE
			;;
		-g | --gff3)
			gff3_file=$VALUE
			;;
		-b | --go_bp)
			go_bp_results=$VALUE
			;;
		-m | --go_mf)
			go_mf_results=$VALUE
			;;
		-c | --go_cc)
			go_cc_results=$VALUE
			;;
		-s | --msigdb)
			msigdb_results=$VALUE
			;;
		-e | --meme_out)
			meme_outdir=$VALUE
			;;	
		-t | --motif_annot)
			motif_annot_file=$VALUE
			;;
		-i | --coverage_profile_option)
			coverage_profile_option=$VALUE
			;;
        	*)
            		echo "ERROR: unknown parameter \"$PARAM\""
            		usage
            		exit 1
            		;;
	esac
	shift
done


if [ -z "${output_filename}" ]; then
	echo "Error: missing argument, provide output filename"
	usage
	exit 1 
fi

if [ -z "${annotation_file}" ]; then
	echo "Error: missing annotation file (output of parse_anot.py)"
	usage
	exit 1
fi

if [ -z "${peaks_file}" ]; then
	echo "Error: missing peaks file (BED format)"
	usage
	exit 1
fi

if [ -z "${gff3_file}" ]; then
	echo "Error: missing gff3 file"
	usage
	exit 1
fi

if [ -z "${go_bp_results}" ]; then
	echo "Error: missing GO term results for Biological Processes"
	usage
	exit 1
fi

if [ -z "${go_mf_results}" ]; then
	echo "Error: missing GO term results for Molecular Functions"
	usage
	exit 1
fi

if [ -z "${go_cc_results}" ]; then
	echo "Error: missing GO term results for Cellular Components"
	usage
	exit 1
fi

if [ -z "${msigdb_results}" ]; then
	echo "Error: missing MSIGDB results"
	usage
	exit 1
fi

if [ -z "${meme_outdir}" ]; then 
	echo "Error: missing MEME output directory"
	usage
	exit 1
fi

if [ -z "${motif_annot_file}" ]; then
	echo "Error: missing motif annotation results"
	usage
	exit 1 
fi
if [ -z "${coverage_profile_option}" ]; then
	echo "Error: missing option for coverage profile calculations"
	usage
	exit 1
fi

Rscript -e "library('rmarkdown'); rmarkdown::render('${report_script}', output_file = '${output_filename}', output_dir='${outdir}', html_document(toc=TRUE, theme='cerulean', number_sections=TRUE, css='${css}', includes=includes(before_body='${header}')))" ${outdir} ${annotation_file} ${peaks_file} ${gff3_file} ${go_bp_results} ${go_mf_results} ${go_cc_results} ${msigdb_results} ${meme_outdir} ${motif_annot_file} ${coverage_profile_option}
