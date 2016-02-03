#! /usr/bin/env python

def check_argv(argv):

	usage = """
	top_motifs.py

	version: 0.0.1
	author: Dilmurat Yusuf

	Usage:

	$ top_motifs.py -m centrimo.html -e centrimo.txt -c cor_flank -a annotation_table -n number_of_top_motifs > annotation_motifs.tsv

	Default number of top motifs according to E-value: 10

	"""

	try:
		opts, args = getopt.getopt(argv,"m:e:c:a:h:n:",[""])

	except getopt.GetoptError:
		print usage
		sys.exit(1)

	if set([opt for opt, arg in opts]) not in [set(["-m", "-e", "-a", "-c"]), set(["-m", "-e", "-a", "-c", '-n'])]:
		print usage
		sys.exit(1)

	number_of_top_motifs = 10

	for opt, arg in opts:
		if opt == '-h' :
			print  usage
			sys.exit()

		if opt == '-m':
			centrimo_html = arg

		if opt == '-e':
			centrimo_txt = arg

		if opt == "-c":
			coordinates_original_flank = arg

		if opt == '-a':
			annotation = arg

		if opt == '-n':
			number_of_top_motifs = int(arg)

	return 	centrimo_html, centrimo_txt, coordinates_original_flank, annotation, number_of_top_motifs

def extract_motif_evalue_id_sites(line, motif_evalue_id_sites):
	if line[0] == " ":
		line = line.strip().split()

		motif_id, evalue, sites = line[1], line[3], line[9]

		evalue_id_sites  =  float(evalue), motif_id, sites

		motif_evalue_id_sites.add(evalue_id_sites)

	return motif_evalue_id_sites

def extract_motif_data(line, on, motif_data):
	if "//@JSON_VAR data" in line:
		on = True

	if  line.endswith("};\n"):
		on = False

	if on:
		motif_data.append(line)

	return on, motif_data

def show_nucleotides(pwm):
		# order of PWM: A C G T

		def assess(site_pwm):
			nucleotide_composition = []
			if site_pwm[0] >= 0.1:
				nucleotide_composition.append("A")
			if site_pwm[1] >= 0.1:
				nucleotide_composition.append("C")
			if site_pwm[2] >= 0.1:
				nucleotide_composition.append("G")
			if site_pwm[3] >= 0.1:
				nucleotide_composition.append("T")
			return nucleotide_composition

		return [tuple(sorted(assess(site_pwm))) for site_pwm in pwm]

def pwm2IUPAC(pwm):

	IUPAC_nucleotide = {
		('A',): 'A',
		('C',): 'C',
		('G',): 'G',
		('T',): 'T',
		('A', 'G') : 'R',
		('C', 'T') : 'Y',
		('C', 'G') : 'S',
		('A', 'T') : 'W',
		('G', 'T') : 'K',
		('A', 'C') : 'M',
		('C', 'G', 'T') : 'B',
		('A', 'G', 'T') : 'D',
		('A', 'C', 'T') : 'H',
		('A', 'C', 'G') : 'V',
		('A', 'C', 'G', 'T') : 'N'
	}

	compositions = show_nucleotides(pwm)

	consensus = [IUPAC_nucleotide[site_composition] for site_composition in compositions]

	return "".join(consensus)

def update_motif_info(motif, motif_info):
	motif_id =  motif['id']
	IUPAC_consensus = pwm2IUPAC(motif['pwm'])
	matched_seq_ids = motif['seqs']

	#there can be duplicated motifs
	#use set to remove duplications
	if motif_id not in motif_info:
		motif_info[motif_id] = IUPAC_consensus, matched_seq_ids

	return motif_info

def update_mapped_coordinates(line, mapped_coordinates):
		line = line[:-1]
		line = line.split()

		seq_id = line[0]
		start_flank = line[1]
		end_flank = line[2]
		start_peak = line[-2]
		end_peak = line[-1]
		strand = line[5]

		coordinates_flank = "%s:%s-%s(%s)" % (seq_id, start_flank, end_flank, strand)
		coordinates_peak = "%s:%s-%s(%s)" % (seq_id, start_peak, end_peak, strand)

		mapped_coordinates[coordinates_flank] = coordinates_peak

		return mapped_coordinates

def update_anot_info(line, anot_info):
		line = line[:-1]
		line = line.split()

		seq_id = line[0]
		start_peak = line[1]
		end_peak = line[2]
		strand = line[5]

		coordinates_peak = "%s:%s-%s(%s)" % (seq_id, start_peak, end_peak, strand)
		annotation = line[8:]

		try:
			anot_info[coordinates_peak].append(annotation)
		except:
			anot_info[coordinates_peak] = [annotation]

		return anot_info

def update_table(evalue_id_sites, motif_info, sequneces, mapped_coordinates, anot_info, table):
	evalue, motif_id, sites = evalue_id_sites
	IUPAC_consensus, matched_seq_ids = motif_info[motif_id]

	for seq_id in matched_seq_ids:
		coordinates_flank = sequneces[seq_id]
		coordinates_peak = mapped_coordinates[coordinates_flank]

		anot_features = anot_info[coordinates_peak]

		chromosome_id, rest = coordinates_peak.split(":")
		coordinate, strand = rest[:-1].split("(")
		start, end = coordinate.split("-")

		for feature in anot_features:
			feature = "\t".join(feature)

			table.append("\t".join([chromosome_id, start, end, strand, feature, IUPAC_consensus, motif_id, sites]))

	return table

if __name__ == '__main__':
	import sys, getopt, json

	argv = sys.argv[1:]

	#check commandline options
	centrimo_html, centrimo_txt, coordinates_peak_flank, annotation, number_of_top_motifs = check_argv(argv)

	###########process for centrimo.txt, start
	#obtain motif id, E-value, site_in_bin
	#sort motif according to E-value and retian motifs according to number_of_top_motifs

	with open(centrimo_txt) as handle:
		motif_evalue_id_sites = set([])

		for line in handle:
			motif_evalue_id_sites = extract_motif_evalue_id_sites(line, motif_evalue_id_sites)

	motif_evalue_id_sites = sorted(motif_evalue_id_sites)

	###########process for centrimo.txt, end

	###########process for centrimo.html, start
	#obtain motif id, pwm, seqs
	#build consensus from pwm

	with open(centrimo_html) as handle:
		on = False
		motif_data = []

		for line in handle:
			on, motif_data = extract_motif_data(line, on, motif_data)

	motif_data = motif_data[2:]
	motif_data[0] = "{\n"
	motif_data.append("}\n")
	motif_data = "".join(motif_data)
	motif_data = json.loads(motif_data)
	#data stucture of motif_data is following:
	#[u'motif_dbs', u'cmd', u'seqlen', u'tested', u'program', u'motifs', u'sequences', u'release', u'sequence_db', u'options', u'revision']
	#u'motifs': [u'peaks', u'score_threshold', u'total_sites', u'db', u'sites', u'len', u'motif_nsites', u'n_tested', u'seqs', u'alt', u'pwm', u'motif_evalue', u'id']

	sequneces = motif_data['sequences']

	motif_info ={}
	for motif in motif_data['motifs']:
		motif_info = update_motif_info(motif, motif_info)
	###########process for centrimo.html, end

	###########coordinate mapping, start
	#map original coordinates and 100-nt coordinates
	mapped_coordinates = {}

	with open(coordinates_peak_flank) as handle:
		for line in handle:
			mapped_coordinates = update_mapped_coordinates(line, mapped_coordinates)
	###########coordinate mapping, end

	###########extract annotation information, start
	#extract anot info: start, end, strand, feature, feature_id
	anot_info = {}

	with open(annotation) as handle:
		handle.next()

		for line in handle:
			anot_info = update_anot_info(line, anot_info)
	###########extract annotation information, end

	#generate table
	header = 'chromosome_id\tstart_position\tend_position\tstrand\tfeature\tfeature_id\tgene_id\ttranscript_id\tgene_type\tgene_name\tgene_strand\tmotif\tmotif_id\tsites'
	table = [header]

	for evalue_id_sites in motif_evalue_id_sites[:number_of_top_motifs]:
		table = update_table(evalue_id_sites, motif_info, sequneces, mapped_coordinates, anot_info, table)

	print "\n".join(table) #_tmp
