configfile: "config.json"

def report_arguments(infile, TRACK_gff, run_motif, run_PATHrich, run_GOrich, run_coverage):
	input_options = {
		"peaks": "NOT_RUN",
		"gff3": "NOT_RUN",
		"go_bp": "NOT_RUN",
		"go_mf": "NOT_RUN",
		"go_cc": "NOT_RUN",
		"msigdb": "NOT_RUN",
		"meme_out": "NOT_RUN",
		"motif_annot": "NOT_RUN"
	}

	coverage_profile = "NOT_RUN"

	if run_coverage:
		input_options["peaks"] = lambda wildcards: infile[wildcards.sample]
		input_options["gff3"] = TRACK_gff
		coverage_profile = "RUN"

	if run_GOrich:
		input_options["go_bp"] = "{sample}-GO-term/BP.GO.results.tsv"
		input_options["go_mf"] = "{sample}-GO-term/MF.GO.results.tsv"
		input_options["go_cc"] = "{sample}-GO-term/CC.GO.results.tsv"

	if run_PATHrich:
		input_options["msigdb"] = "{sample}.msigdb.results.tsv"

	if run_motif:
		input_options["meme_out"] = "{sample}_memechip_output/meme_out/"
		input_options["motif_annot"] = "{sample}.anot-motif.tsv"

	return input_options, coverage_profile

RCAS_path = config["RCAS_path"]
anot = RCAS_path  + "/src/RCAS.anot"
motif = RCAS_path  + "/src/RCAS.motif"
GOrich = RCAS_path  + "/src/RCAS.GOrich"
PATHrich = RCAS_path  + "/src/RCAS.PATHrich"

TRACK_gff = config["gff3"]
genome_reference = config["genome"]
species = config["species"]

if species == "human":
	gmt = "c2.cp.v5.0.entrez.gmt"
elif species == "fly":
	gmt = "c2.cp.v5.0.entrez.dm3.gmt"
elif species == "worm":
	gmt = "c2.cp.v5.0.entrez.ce10.gmt"
elif species == "mouse":
	gmt = "c2.cp.v5.0.entrez.mm9.gmt"

infile = config["infile"]

run_motif = eval(config["switch"]["run_motif"])
run_PATHrich = eval(config["switch"]["run_PATHrich"])
run_GOrich = eval(config["switch"]["run_GOrich"])
run_coverage = eval(config["switch"]["run_coverage"])

input_options, coverage_profile = report_arguments(infile, TRACK_gff,
								run_motif, run_PATHrich,
								run_GOrich, run_coverage)

rule target:
	 input:
		   expand("{sample}.rcas.html", sample=infile)

#default step
include: anot

#optionl step
if run_motif:
	include: motif

#optional step
if run_PATHrich:
	include: PATHrich

#optional step
if run_GOrich:
	include: GOrich

#a temporary file used to switch on/off modules
rule NOT_RUN:
	output:
		temp("NOT_RUN")
	shell:
		"touch {output}"

#default step
rule html_report:
	 input:
			annot="{sample}.anot.tsv",
			peaks=input_options["peaks"],
			gff3=input_options["gff3"],
			go_bp=input_options["go_bp"],
			go_mf=input_options["go_mf"],
			go_cc=input_options["go_cc"],
			msigdb=input_options["msigdb"],
			meme_out=input_options["meme_out"],
			motif_annot=input_options["motif_annot"]
	 output:
			"{sample}.rcas.html"
	 shell:
			"bash {RCAS_path}/src/generate_report.sh"
			" --output_filename={output} --annot={input.annot}"
			" --peaks={input.peaks}"
			" --gff3={input.gff3}"
			" --go_bp={input.go_bp}"
			" --go_mf={input.go_mf}"
			" --go_cc={input.go_cc}"
			" --msigdb={input.msigdb}"
			" --meme_out={input.meme_out} --motif_annot={input.motif_annot}"
			" --coverage_profile_option=%s" % coverage_profile
