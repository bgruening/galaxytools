#! /usr/bin/env python
import math
import optparse
import os
import string
import sys
import tempfile

import fasta


def format(seq, N=60):
    nseg = int(math.ceil(len(seq) / (N + 0.0)))
    return "\n".join([seq[i * N : (i + 1) * N] for i in range(nseg)])


# write into fasta format file

parser = optparse.OptionParser(version="%prog ")

parser.add_option(
    "-i",
    "--input",
    dest="input_fasta",
    action="store",
    help="name of input file in fasta format",
)
parser.add_option(
    "-L",
    "--LibHmm",
    dest="hmm_path",
    action="store",
    default="HMM3",
    help="path of hmm database",
)
parser.add_option(
    "--gff", dest="out_gff", action="store", help="name of output gff file"
)
parser.add_option(
    "--seq", dest="out_seq", action="store", help="name of output sequence file"
)
parser.add_option(
    "--mask", dest="out_mask", action="store", help="name of output mask file"
)
parser.add_option(
    "-k",
    "--kingdoms",
    dest="kingdoms",
    action="store",
    default="arc,bac,euk",
    help="kingdom used",
)
parser.add_option(
    "-m",
    "--moltypes",
    dest="moltypes",
    action="store",
    default="lsu,ssu,tsu",
    help="molecule type detected",
)
parser.add_option(
    "-e",
    "--Evalue",
    dest="evalue",
    action="store",
    type="float",
    default=0.01,
    help="evalue cut-off for hmmsearch",
)
parser.add_option(
    "--cpu",
    dest="cpu",
    action="store",
    type="int",
    default=4,
    help="number of cpus hmmsearch should use",
)


try:
    (options, args) = parser.parse_args()
except Exception:
    parser.print_help()
    sys.exit(1)

if options.input_fasta is None or options.hmm_path is None:
    parser.print_help()
    sys.exit(1)

# print "%s"% os.path.abspath(options.hmm_path)
# os.environ["HMMERDB"] += ":"+os.path.abspath(options.hmm_path)
# print os.environ["HMMERDB"]
fname = os.path.abspath(options.input_fasta)

tr = string.maketrans(
    "gatcryswkmbdhvnGATCRYSWKMBDHVN", "ctagyrswmkvhdbnCTAGYRSWMKVHDBN"
)


def rev_record(record):
    return ">" + record.header + "|rev\n" + format(record.sequence[::-1].translate(tr))


records = [rec for rec in fasta.fasta_itr(fname)]
headers = [[rec.header, len(rec.sequence)] for rec in records]


temp_fasta = tempfile.NamedTemporaryFile(delete=False)
ff = open(temp_fasta.name, "w")
for (i, rec) in enumerate(records):
    ff.write(">s" + str(i) + "\n" + format(rec.sequence) + "\n")
    ff.write(">s" + str(i) + "|rev\n" + format(rec.sequence[::-1].translate(tr)) + "\n")
ff.close()
# sys.exit(1)
# a temporary fasta file, use s(int) to easy the parsing


def parse_hmmsearch(kingdom, moltype, src):
    # function to parse hmmsearch output
    resu = []
    data = open(src).readlines()
    inds = [-1] + [i for (i, x) in enumerate(data[2]) if x == " "]
    inds = [(inds[j] + 1, inds[j + 1]) for j in range(len(inds) - 1)]
    data = [line for line in data if line[0] != "#"]
    for line in data:
        if not len(line.strip()):
            continue
        [
            read,
            acc,
            tlen,
            qname,
            qaccr,
            qlen,
            seq_evalue,
            seq_score,
            seq_bias,
            seq_num,
            seq_of,
            dom_cEvalue,
            dom_iEvalue,
            dom_score,
            dom_bias,
            hmm_start,
            hmm_end,
            dom_start,
            dom_end,
            env_start,
            env_end,
        ] = line.split()[:21]
        #            [line[x[0]:x[1]].strip() for x in inds[:21]]
        if string.atof(dom_iEvalue) < options.evalue:
            #            resu.append("\t".join([read, acc, tlen, qname, qaccr, \
            #                    qlen, seq_evalue, seq_score, seq_bias, seq_num, seq_of, \
            #                    dom_cEvalue, dom_iEvalue, dom_score, dom_bias, hmm_start, \
            #                    hmm_end, dom_start, dom_end, env_start, env_end]))
            resu.append("\t".join([qname, dom_start, dom_end, read, dom_iEvalue]))

    #    print resu[0]
    #    print resu[-1]
    return resu


hmm_resu = []
for kingdom in options.kingdoms.split(","):
    for moltype in options.moltypes.split(","):
        # print kingdom, moltype
        # hmm_out_fname = "%s.%s_%s.out" % (out_fname, kingdom, moltype)
        temp_hmm_out = tempfile.NamedTemporaryFile(delete=False)
        # dom_out_fname = "%s.%s_%s.dom" % (out_fname, kingdom, moltype)
        temp_dom_out = tempfile.NamedTemporaryFile(delete=False)
        # cmd = '/home/thumper6/camera-annotation/sitao/hmmer3.0/hmmer-3.0-linux-intel-x86_64/binaries/hmmsearch --cpu 1 -o %s --domtblout %s -E %g %s/%s_%s.hmm %s' % \
        cmd = "hmmsearch --cpu %s -o %s --domtblout %s -E %g %s/%s_%s.hmm %s" % (
            options.cpu,
            temp_hmm_out.name,
            temp_dom_out.name,
            options.evalue,
            os.path.abspath(options.hmm_path),
            kingdom,
            moltype,
            temp_fasta.name,
        )
        #        hmm_resu += parse_hmmsearch(os.popen(cmd))
        os.system(cmd)
        hmm_resu += parse_hmmsearch(kingdom, moltype, temp_dom_out.name)
        os.remove(temp_hmm_out.name)
        os.remove(temp_dom_out.name)

dict_read2kingdom = {}
for line in hmm_resu:
    [feature_type, r_start, r_end, read, evalue] = line.strip().split("\t")
    read = read.split("|")[0]
    evalue = string.atof(evalue)
    kingdom = feature_type.split("_")[0]
    if read in dict_read2kingdom:
        if evalue < dict_read2kingdom[read][1]:
            dict_read2kingdom[read] = [kingdom, evalue]
    else:
        dict_read2kingdom[read] = [kingdom, evalue]

header = ["##seq_name", "method", "feature", "start", "end", "evalue", "strand", "gene"]
ff = open(options.out_gff, "w")
dict_rRNA = {
    "arc_lsu": "Archaeal:23S_rRNA",
    "arc_ssu": "Archaeal:16S_rRNA",
    "arc_tsu": "Archaeal:5S_rRNA",
    "bac_lsu": "Bacterial:23S_rRNA",
    "bac_ssu": "Bacterial:16S_rRNA",
    "bac_tsu": "Bacterial:5S_rRNA",
    "euk_lsu": "Eukaryotic:28S_rRNA",
    "euk_ssu": "Eukaryotic18S_rRNA",
    "euk_tsu": "Eukaryotic:8S_rRNA",
}

ff.write("\t".join(header) + "\n")
for line in hmm_resu:
    #    [kingdom, moltype, read, acc, tlen, qname, qaccr, \
    #         qlen, seq_evalue, seq_score, seq_bias, seq_num, seq_of, \
    #         dom_cEvalue, dom_iEvalue, dom_score, dom_bias, hmm_start, \
    #         hmm_end, dom_start, dom_end, env_start, env_end] = line.strip().split('\t')
    [feature_type, r_start, r_end, read, evalue] = line.strip().split("\t")
    if dict_read2kingdom[read.split("|")[0]][0] != feature_type.split("_")[0]:
        continue
    feature_type = dict_rRNA[feature_type]
    if read.endswith("|rev"):
        strand = "-"
        tmp = map(string.atoi, [r_start, r_end])
        pos = string.atoi(read[1:-4])
        header = headers[pos][0]
        L = headers[pos][1]
        [r_end, r_start] = [str(L + 1 - x) for x in tmp]
    else:
        strand = "+"
        pos = string.atoi(read[1:])
        header = headers[pos][0]
    ff.write(
        "\t".join(
            [
                header.split()[0],
                "rna_hmm3",
                "rRNA",
                r_start,
                r_end,
                evalue,
                strand,
                feature_type,
            ]
        )
        + "\n"
    )
ff.close()


os.remove(temp_fasta.name)

# postprocessing
rRNA_dict = {}

if options.out_gff:
    for line in open(options.out_gff).readlines()[1:]:
        line = line.strip().split()
        # print '%s %s'% (line[0],line[1])
        if rRNA_dict.get(line[0], None) is None:
            rRNA_dict[line[0]] = []
        # else:
        rRNA_dict[line[0]].append(
            (string.atoi(line[3]), string.atoi(line[4]), line[6], line[7])
        )

if options.out_mask:
    f_mask = open(options.out_mask, "w")

if options.out_seq:
    f_seq = open(options.out_seq, "w")
for rec in fasta.fasta_itr(options.input_fasta):
    header = rec.header.split()[0]
    seq = rec.sequence[:]
    tno = 1
    for (start, end, strand, rRNA_type) in rRNA_dict.get(header, []):
        seq = seq[: (start - 1)] + "N" * (end - start + 1) + seq[end:]

        # print "%s %d"% (header,tno)
        if options.out_seq:
            f_seq.write(
                ">%s.%d /start=%d /end=%d /strand=%s %s\n"
                % (header, tno, start, end, strand, rRNA_type)
            )
            if strand == "+":  # forward strand
                f_seq.write(format(rec.sequence[(start - 1) : end]) + "\n")
            else:
                f_seq.write(
                    format(rec.sequence[(start - 1) : end][::-1].translate(tr)) + "\n"
                )
        tno = tno + 1

    # f_mask.write(">"+header+"\n")
    # next line by liwz
    if options.out_mask:
        f_mask.write(">" + rec.header + "\n")
        f_mask.write(format(seq) + "\n")

if options.out_mask:
    f_mask.close()
if options.out_seq:
    f_seq.close()
