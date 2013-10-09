#!/usr/bin/env python
## Copyright (c) 2010 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

##Functions necessary for this script

import linecache, cPickle

DEBUG = True


def invalidoptions(argument):
  if len(argument) > 0:
    print >> sys.stderr, "Invalid options input:"
    print >> sys.stderr, argument
  print "From the command line, input antismash --help for more information."
  logfile.write("Invalid options input: " + argument + "\n")
  logfile.close()
  sys.exit(1)

def sortdictkeysbyvalues(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    return [key for value, key in items]

def sortdictkeysbyvaluesrev(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]

def sortdictkeysbyvaluesrevv(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    items.reverse()
    return [value for value, key in items]

def get_sequence(fasta):
    """get the description and trimmed dna sequence"""
    #in_file = open(fasta, 'r')
    #content = in_file.readlines()
    #in_file.close()
    #content2 = []
    #for i in content:
      #if i != "":
      #  content2.append(i)
    content = []
    [content.append(line) for line in open(fasta, 'r') if line]
    #content = content2
    while content[0] == "" or content[0] == "\n":
      content = content[1:]
    header = content[0]
    content = content[1:]
    content = [x.rstrip() for x in content]
    seq = "".join(content)
    if ">" not in header or ">" in seq:
      print >> sys.stderr, "FASTA file not properly formatted; should be single sequence starting with '>' and sequence name."
      logfile.write("FASTA file not properly formatted; should started with '>' and sequence name on first line.\n")
      logfile.close()
      sys.exit(1)
    return seq

def complement(seq):  
  complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}  
  complseq = []
  for base in seq:
    if base in complement.keys():
      complbase = complement[str(base)]
      complseq.append(complbase)
    else:
      complbase = 'n'
      complseq.append(complbase)
  return complseq 

def reverse_complement(seq):  
  seq = list(seq)  
  seq.reverse()  
  revcompl = complement(seq)
  revcomplstr = str()
  for i in revcompl:
    revcomplstr = revcomplstr + str(i)
  return  revcomplstr

def fastaseqlengths(proteins):
  names = proteins[0]
  seqs = proteins[1]
  seqlengths = {}
  a = 0
  for i in names:
    #seq = seqs[a]
    #seqlength = len(seq)
    #seqlengths[i] = seqlength
    seqlengths[i] = len(seqs[a])
    a += 1
  return seqlengths

# Function that reads the fasta file into a dictionary
def fastadict(fasta):
  file = open(fasta,"r")
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  filetext = filetext.strip()
  #Replaces all spaces with "_" to avoid problems
  filetext = filetext.replace(' ','_')
  filetext = filetext.split()
  dictseq = {}
  for a in filetext:
    if ">" in a[0]:
      f = str()
      d = a[1:68]
    else:
      e = a
      f += e
      dictseq[d] = f
  return dictseq

# Function that extracts all sequence names from the fasta dictionary
def lnames(fastadict):
  items = fastadict.items()
  items.sort()
  return [names for names, seqs in items]

# Function that extracts all sequences from the fasta dictionary
def lseqs(fastadict):
  items = fastadict.items()
  items.sort()
  return [seqs for names, seqs in items]

def extractpositions(refmusclefile,newmusclefile,positions,refsequencename,querysequencename):
  dict = fastadict(refmusclefile)
  seqs = lseqs(dict)
  names = lnames(dict)
  #startpos = 2
  residues = []
  #Count residues in ref sequence and put positions in list
  muscle_dict = fastadict(newmusclefile)
  muscle_seqs = lseqs(muscle_dict)
  muscle_names = lnames(muscle_dict)
  refseqnr = muscle_names.index(refsequencename)
  #Extract activity signature
  refseq = muscle_seqs[refseqnr]
  poslist = []
  b = 0
  c = 0
  while refseq != "":
    i = refseq[0]
    if c in positions and i != "-":
      poslist.append(b)
    if i != "-":
      c += 1
    b += 1
    refseq = refseq[1:]
  #Extract positions from query sequence
  query_seqnr = muscle_names.index(querysequencename)
  query_seq = muscle_seqs[query_seqnr]
  for j in poslist:
    residues.append(query_seq[j])
  return residues

def parsegenes(genes):
  genedict = {}
  genelist = []
  joinlist = []
  joindict = {}
  accessiondict = {}
  error = "n"
  errorlocations = []
  genenr = 0
  for i in genes:
    if "     gene            " in i:
      i = i.split("     gene            ")[0]
    elif "FT   gene            " in i:
      i = i.split("FT   gene            ")[0]
    join = "no"
    genenr += 1
    #Find gene location info for each gene
    if "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("   /")[0]
      while ")" not in location.replace(" ","")[-3:]:
        location = location.rpartition("\n")[0]
      location = location.replace("\n","")
      location = location.replace(" ","")
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("/")[0]
      while ")" not in location.replace(" ","")[-3:]:
        location = location.rpartition("\n")[0]
      location = location.replace("\n","")
      location = location.replace(" ","")
    else:
      location = i.split("\n")[0]
    original_location = location
    #location info found in gbk/embl file, now extract start and end positions
    if location.count("(") != location.count(")"):
      error = "y"
      errorlocations.append(original_location)
      continue
    if "join(complement" in location.lower():
      location = location.lower()
      join = "yes"
      location2 = location.partition("join(")[2][:-1].replace("<","").replace(">","")
      if ("complement(" in location2[0:12] and location2[-1] != ")") or ")," in location2:
        error = "y"
        errorlocations.append(original_location)
        continue
      elif ("complement(" in location2[0:12] and location2[-1] == ")" and location2[12:-2].count(")") == 0 and location2[12:-2].count("(") == 0):
        location2 = location2.partition("complement(")[2][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
        strand = "-"
      else:
        error = "y"
        errorlocations.append(original_location)
        continue
    elif "complement" in location.lower():
      location = location.lower()
      location = location.partition("complement(")[2][:-1]
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.partition("join(")[2][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","")
        if ".." in location:
          end = location.split("..")[1]
        else:
          end = location
        end = end.replace(">","")
      strand = "-"
    else:
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.partition("join(")[2][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","")
        if ".." in location:
          end = location.split("..")[1]
        else:
          end = location
        end = end.replace(">","")
      strand = "+"
    try:
      if int(start) > int(end):
        start2 = end
        end2 = start
        start = start2
        end = end2
    except ValueError:
        error = "y"
        errorlocations.append(original_location)
        continue
    #Correct for alternative codon start positions
    if "codon_start=" in i.lower():
      temp = i.lower().split("codon_start=")[1].split()[0]
      if '"' in temp:
        # temp ist "1" oder "2", dies kommt aus biopython
        temp = temp[1]
      else:
        # ohne anfuhrungszeichen ... 1 oder 2
        temp = temp[0]
      codonstart = temp
      if strand == "+":
        start = str(int(start) +  (int(codonstart) - 1))
      elif strand == "-":
        end = str(int(end) - (int(codonstart) - 1))
    #Find gene name for each gene, preferably locus_tag, than gene, than protein_ID
    a = 0
    b = 0
    genename = ""
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      if "protein_id=" in line:
        genename = (line.split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif "protein_id=" in line.lower():
        genename = (line.lower().split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        genename = ""
        b += 1
      else:
        a += 1
    if len(genename) > 1:
      accnr = genename
    else:
      accnr = "no_accession_number_found"
    a = 0
    b = 0
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      if "gene=" in line:
        genename = (line.split("gene=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif "gene=" in line.lower():
        genename = (line.lower().split("gene=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        b += 1
      else:
        a += 1
    a = 0
    b = 0
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      if "locus_tag=" in line:
        genename = (line.split("locus_tag=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif "locus_tag=" in line.lower():
        genename = (line.lower().split("locus_tag=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        if genename == "":
          genename = "prot_ID_" + str(genenr)
        b += 1
      else:
        a += 1
    #Find sequence for each gene
    a = 0                                             ###Not all gbks contain protein sequences as translations, therefore sequences from gene clusters are now extracted from the database at a later stage if sequence is not in gbk
    b = 0
    sequence = ""
    while b < 2:
      line = i.split("\n")[a]
      if "translation=" in line:
        sequence = line.split("translation=")[1][1:]
        b += 1
        a += 1
        if line.count('"') > 1:
          sequence = line.split("translation=")[1][1:-1]
          b = 2
      elif "translation=" in line.lower():
        sequence = line.lower().split("translation=")[1][1:]
        b += 1
        a += 1
        if line.count('"') > 1:
          sequence = line.lower().split("translation=")[1][1:-1]
          b = 2
      elif a == (nrlines - 2) or a == (nrlines - 1):
        sequence = ""
        b = 2
      elif b == 1:
        if '"' in line:
          seqline = line.replace(" ","")
          seqline = seqline.split('"')[0]
          sequence = sequence + seqline
          b += 1
        else:
          seqline = line.replace(" ","")
          sequence = sequence + seqline
        a += 1
      else:
        a += 1
    sequence = sequence.upper()
    #Quality-check sequence
    forbiddencharacters = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for z in forbiddencharacters:
      if z in sequence:
        sequence = ""
    #Find annotation for each gene
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      if "product=" in line:
        annotation = line.split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif "product=" in line.lower():
        annotation = line.lower().split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif a == (nrlines - 1):
        annotation = "not_annotated"
        b += 1
      else:
        a += 1
    accessiondict[genename] = accnr
    if join == "yes":
      joinlist.append(genename)
      joindict[genename] = joinedparts2
    #Save data to dictionary
    if len(genename) > 1:
      genedict[genename] = [start,end,strand,annotation,sequence]
    genelist.append(genename)
  if error == "y":
    errorinfo = "\n".join(errorlocations)
    print >> sys.stderr, "Exit: locations in GBK/EMBL file not properly formatted:\n" + errorinfo
    logfile.write("Exit: GBK file not properly formatted, no sequence found or no CDS annotation found.\n")
    logfile.close()
    sys.exit(1)
  return [genelist, genedict, joinlist, joindict, accessiondict]

def cleandnaseq(dnaseq):
  dnaseq = dnaseq.replace(" ","")
  dnaseq = dnaseq.replace("\t","")
  dnaseq = dnaseq.replace("\n","")
  dnaseq = dnaseq.replace("0","")
  dnaseq = dnaseq.replace("1","")
  dnaseq = dnaseq.replace("2","")
  dnaseq = dnaseq.replace("3","")
  dnaseq = dnaseq.replace("4","")
  dnaseq = dnaseq.replace("5","")
  dnaseq = dnaseq.replace("6","")
  dnaseq = dnaseq.replace("7","")
  dnaseq = dnaseq.replace("8","")
  dnaseq = dnaseq.replace("9","")
  dnaseq = dnaseq.replace("/","")
  dnaseq = dnaseq.replace("u","t")
  dnaseq = dnaseq.replace("U","T")
  dnaseq = dnaseq.replace("r","n")
  dnaseq = dnaseq.replace("R","n")
  dnaseq = dnaseq.replace("y","n")
  dnaseq = dnaseq.replace("Y","n")
  dnaseq = dnaseq.replace("w","n")
  dnaseq = dnaseq.replace("W","n")
  dnaseq = dnaseq.replace("s","n")
  dnaseq = dnaseq.replace("S","n")
  dnaseq = dnaseq.replace("m","n")
  dnaseq = dnaseq.replace("M","n")
  dnaseq = dnaseq.replace("k","n")
  dnaseq = dnaseq.replace("K","n")
  dnaseq = dnaseq.replace("h","n")
  dnaseq = dnaseq.replace("H","n")
  dnaseq = dnaseq.replace("b","n")
  dnaseq = dnaseq.replace("B","n")
  dnaseq = dnaseq.replace("v","n")
  dnaseq = dnaseq.replace("V","n")
  dnaseq = dnaseq.replace("d","n")
  dnaseq = dnaseq.replace("D","n")
  return dnaseq

def extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict):
  names = []
  seqs = []
  for i in genelist:
    genename = i
    #If suitable translation found in gbk, use that
    if len(genedict[i][4]) > 5:
      protseq = genedict[i][4]
      i = genedict[i]
    #If no suitable translation found in gbk, extract from DNA sequence
    else:
      i = genedict[i]
      y = int(i[0])
      z = int(i[1])
      if i[2] == "+":
        if genename in joinlist:
          geneseq = ""
          for j in joindict[genename]:
            partstart = int(j.split("..")[0])
            if ".." in j:
              partend = int(j.split("..")[1])
            else:
              partend = int(j)
            geneseqpart = dnaseq[(partstart - 1):partend]
            geneseq = geneseq + geneseqpart
        else:
          geneseq = dnaseq[(y - 1):z]
        protseq = translate(geneseq)
      elif i[2] == "-":
        if genename in joinlist:
          geneseq = ""
          joinlistrev = joindict[genename]
          joinlistrev.reverse()
          for j in joinlistrev:
            partstart = int(j.split("..")[0])
            if ".." in j:
              partend = int(j.split("..")[1])
            else:
              partend = int(j)
            geneseqpart = rc_dnaseq[(len(rc_dnaseq) - partend):(len(rc_dnaseq) - partstart + 1)]
            geneseq = geneseq + geneseqpart
        else:
          geneseq = rc_dnaseq[(len(rc_dnaseq) - z):(len(rc_dnaseq) - y + 1)]
        protseq = translate(geneseq)
    name = "input" + "|" + "c1" + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename + "|" + i[3]
    seqs.append(protseq)
    names.append(name)
  proteins = [names,seqs,genelist,genedict,accessiondict]
  return proteins

def gbk2proteins(gbkfile):
  file = open(gbkfile,"r")
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  if "     CDS             " not in filetext or "\nORIGIN" not in filetext:
    print >> sys.stderr, "Exit: GBK file not properly formatted, no sequence found or no CDS annotation found."
    logfile.write("Exit: GBK file not properly formatted, no sequence found or no CDS annotation found.\n")
    logfile.close()
    sys.exit(1)
  cdspart = filetext.split("\nORIGIN")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = filetext.split("\nORIGIN")[1]
  dnaseq = cleandnaseq(dnaseq)
  sequence = dnaseq
  if (sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
    print >> sys.stderr, "Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file."
    sys.exit(1)
  dnaseqlength = len(dnaseq)
  rc_dnaseq = reverse_complement(dnaseq)
  #Extract genes
  genes = cdspart.split("     CDS             ")
  genes = genes[1:]
  try:
    genesdetails = parsegenes(genes)
  except ValueError, e:
    print >> sys.stderr, "Could not parse genes from GBK/EMBL file. Please check if your GBK/EMBL file is valid."
    raise
    print >> sys.stderr, "Error was: %s" % e
    print len(genes)
    sys.exit(1)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  #Locate all genes on DNA sequence and translate to protein sequence
  proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict)
  textlines = filetext.split("\n//")[0]
  textlines = textlines.split("\n")
  accession = ""
  for i in textlines:
    if accession == "":
      if "LOCUS       " in i:
        j = i.split("LOCUS       ")[1]
        accession = j.split(" ")[0]
        if len(accession) < 4:
          accession = ""
  #Test if accession number is probably real GenBank/RefSeq acc nr
  numbers = range(0,10)
  letters = []
  for i in ascii_letters:
    letters.append(i)
  nrnumbers = 0
  nrletters = 0
  for i in accession:
    if i in letters:
      nrletters += 1
    try:
      j = int(i)
      if j in numbers:
        nrnumbers += 1
    except:
      pass
  if nrnumbers < 3 or nrletters < 1:
    accession = ""
  return [proteins,accession,dnaseqlength]

def embl2proteins(emblfile,sequence):
  file = open(emblfile,"r")
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  file.close()
  if "FT   CDS " not in filetext or ("\nSQ" not in filetext and len(sequence) < 1):
    logfile.write("Exit: EMBL file not properly formatted, no sequence found or no CDS annotation found.\n")
    print >> sys.stderr, "Exit: EMBL file not properly formatted, no sequence found or no CDS annotation found.\n"
    logfile.close()
    sys.exit(1)
  cdspart = filetext.split("\nSQ  ")[0]
  #Extract DNA sequence and calculate reverse complement of it
  seqpart = filetext.split("\nSQ  ")[1]
  seqlines = seqpart.split("\n")[1:]
  dnaseq = ""
  for i in seqlines:
    dnaseq = dnaseq + i
  dnaseq = cleandnaseq(dnaseq)
  sequence = dnaseq
  if (sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
    print >> sys.stderr, "Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file."
    sys.exit(1)
  dnaseqlength = len(dnaseq)
  rc_dnaseq = reverse_complement(dnaseq)
  #Extract genes
  genes = cdspart.split("FT   CDS             ")
  genes = genes[1:]
  try:
    genesdetails = parsegenes(genes)
  except ValueError, e:
    print >> sys.stderr, "Could not parse genes from GBK/EMBL file. Please check if your GBK/EMBL file is valid."
    print >> sys.stderr, "Error was: %s" % e
    sys.exit(1)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  #Locate all genes on DNA sequence and translate to protein sequence
  proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict)
  textlines = filetext.split("SQ   ")[0]
  textlines = textlines.split("\n")
  accession = ""
  for i in textlines:
    if accession == "":
      if "AC   " in i:
        j = i.split("AC   ")[1]
        j = j.replace(" ","")
        accession = j.split(";")[0]
        if len(accession) < 4:
          accession = ""
  #Test if accession number is probably real GenBank/RefSeq acc nr
  numbers = range(0,10)
  letters = []
  for i in ascii_letters:
    letters.append(i)
  nrnumbers = 0
  nrletters = 0
  for i in accession:
    if i in letters:
      nrletters += 1
    try:
      j = int(i)
      if j in numbers:
        nrnumbers += 1
    except:
      pass
  if nrnumbers < 3 or nrletters < 1:
    accession = ""
  return [proteins,accession,dnaseqlength]

def translate(sequence):
  #Translation table standard genetic code; according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
  transldict = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 
                 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', 
                 'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 
                 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 
                 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 
                 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 
                 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 
                 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
                 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 
                 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
                 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
                 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 
                 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 
                 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
                 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
                 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                 'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                 'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
                 'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
                 'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
                 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
                 'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
                 'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
                 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
                 'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
                 'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
                 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
                 'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
                 'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
                 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
                 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
  triplets = []
  triplet = ""
  a = 0
  for i in sequence:
    if a < 2:
      a += 1
      triplet = triplet + i
    elif a == 2:
      triplet = triplet + i
      triplets.append(triplet)
      triplet = ""
      a = 0
  protseq = ""
  aanr = 0
  for i in triplets:
    aanr += 1
    if aanr == 1:
      protseq = protseq + "M"
    else:
      if "n" in i or "N" in i or i not in transldict.keys():
        protseq = protseq + "X"
      else:
        protseq = protseq + transldict[i]
  if  len(protseq) > 0 and protseq[-1] == "*":
    protseq = protseq[:-1]
  return protseq

def writefasta(names,seqs,file):
  e = 0
  f = len(names) - 1
  try:
    out_file = open(file,"w")
    while e <= f:
      out_file.write(">%s\n%s\n" % (names[e], seqs[e]) )
      #out_file.write(">")
      #out_file.write(names[e])
      #out_file.write("\n")
      #out_file.write(seqs[e])
      #out_file.write("\n")
      e += 1
    out_file.close()
  except(IOError,OSError,NotImplementedError):
    print >> sys.stderr, "FASTA file not created."
    logfile.write("FASTA file not created.\n")

def parsehmmoutput(cutoff,file):
  #file = open(file,"r")
  #filetext = file.read()
  #filetext = filetext.replace("\r","\n")
  #lines = filetext.split("\n")
  protlines = []
  #for i in lines:
  #  if len(i) > 1 and i[0] != "#":
  #    protlines.append(i)
  [protlines.append(line.strip()) for line in open(file,"r") if len(line) > 1 and not line.startswith('#')]
  proteins = []
  scores = []
  #measuringline = lines[2]
  measuringline = linecache.getline(file, 3)
  x = 0
  y = 0
  for i in measuringline:
    y += 1
    if "-" in i:
      x += 1
    else:
      if x > 1:
        break
  for i in protlines:
    #accession = ""
    #a = 0
    protname = i[0:y]
    protnameparts = protname.split("|")
    accession = protnameparts[4]
    score = i[(y+76):(y+82)]
    score = float(score.replace(" ",""))
    if score > cutoff and len(accession) > 1:
      proteins.append(accession)
      scores.append(score)
  return [proteins,scores]

def sortonsecondvalueoflist(first,second):
  f = int(first[1])
  s = second[1]
  if f > s:
    value = 1
  elif f < s:
    value = -1
  elif f == s:
    value = 0
  return value

def hmmlengths(hmmfile):
  hmmlengthsdict = {}
  file = open(hmmfile,"r")
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  hmms = filetext.split("//")[:-1]
  for i in hmms:
    namepart = i.split("NAME  ")[1]
    name = namepart.split("\n", 1)[0]
    lengthpart = i.split("LENG  ")[1]
    #print lengthline
    #tabs = lengthline.split(" ")
    #tabs2 = []
    #for j in tabs:
    #  if j != "":
    #    tabs2.append(j)
    #print tabs2
    length = lengthpart.split("\n", 1)[0]
    hmmlengthsdict[name] = int(length)
  return hmmlengthsdict

def hmmscanparse(hmmscanoutputfile,hmmlengthsdict):
  domaindict = {}
  file = open(hmmscanoutputfile,"r")
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  outputs = filetext.split("Query:       ")[1:]
  for i in outputs:
    protname = i.split("\n", 1)[0]
    protname = protname.split(" ", 1)[0]
    domainresults = i.split("Domain annotation for each model:\n")[1]
    domainresults = domainresults.split("\n\nInternal pipeline statistics summary:")[0]
    domains = domainresults.split(">> ")
    domainlist = []
    #Find all domains
    for i in domains:
      tokens = i.split('\n')
      domainname = tokens[0]
      domainname = domainname.split(" ", 1)[0]
      domainresults = tokens[3:-2]
      for i in domainresults:
        tabs = i.split(" ")
        tabs2 = []
        [tabs2.append(tab) for tab in tabs if tab != '']
        #for i in tabs:
        #  if i != "":
        #    tabs2.append(i)
        tabs = tabs2
        start = int(tabs[12])
        end = int(tabs[13])
        evalue = tabs[5]
        score = float(tabs[2])
        domainlist.append([domainname,start,end,evalue,score])
    domainlist.sort(sortonsecondvalueoflist)
    #Purify domain list to remove overlapping domains, only keeping those with the highest scores
    if len(domainlist) > 1:
      domainlist2 = [domainlist[0]]
      for i in domainlist[1:]:
        maxoverlap = 20
        if i[1] < (domainlist2[-1][2] - maxoverlap):
          if i[4] < domainlist2[-1][4]:
            pass
          elif i[4] > domainlist2[-1][4]:
            del domainlist2[-1]
            domainlist2.append(i)
        else:
          domainlist2.append(i)
      domainlist = domainlist2
    #Merge domain fragments which are really one domain
    if len(domainlist) > 1:
      domainlist2 = [domainlist[0]]
      for i in domainlist[1:]:
        alilength1 = int(domainlist2[-1][2]) - int(domainlist2[-1][1])
        alilength2 = int(i[2]) - int(i[1])
        domainlength = hmmlengthsdict[i[0]]
        if i[0] == domainlist2[-1][0] and (alilength1 < (0.75 * domainlength) or alilength2 < (0.75 * domainlength)) and (alilength1 + alilength2) < (1.5 * domainlength):
          name = i[0]
          start = domainlist2[-1][1]
          end = i[2]
          evalue = str(float(domainlist2[-1][3]) * float(i[3]))
          score = str(float(domainlist2[-1][4]) + float(i[4]))
          del domainlist2[-1]
          domainlist2.append([name,start,end,evalue,score])
        else:
          domainlist2.append(i)
      domainlist = domainlist2
    #Remove incomplete domains (covering less than 60% of total domain hmm length)
    if len(domainlist) > 1:
      domainlist2 = []
      for i in domainlist:
        alilength = int(i[2]) - int(i[1])
        domainlength = hmmlengthsdict[i[0]]
        if alilength > (0.6 * domainlength):
          domainlist2.append(i)
      domainlist = domainlist2
    #Save domainlist to domaindict
    domaindict[protname] = domainlist
  return domaindict

def blastparse(blasttext,minseqcoverage,minpercidentity,seqlengths,geneclustergenes):
  blastdict = {}
  querylist = []
  hitclusters = []
  blastlines = blasttext.split("\n")[:-1]
  #Filter for best blast hits (of one query on each subject)
  query_subject_combinations = []
  blastlines2 = []
  for i in blastlines:
    tabs = i.split("\t")
    query = tabs[0]
    subject = tabs[1]
    query_subject_combination = query + "_" + subject
    if query_subject_combination in query_subject_combinations:
      pass
    else:
      query_subject_combinations.append(query_subject_combination)
      blastlines2.append(i)
  blastlines = blastlines2
  #Filters blastlines to get rid of hits that do not meet criteria
  blastlines2 = []
  for i in blastlines:
    tabs = i.split("\t")
    perc_ident = int(tabs[2].split(".",1)[0])
    alignmentlength = float(tabs[3])
    evalue = str(tabs[10])
    blastscore = int(tabs[11].split(".",1)[0])
    if seqlengths.has_key(query):
      perc_coverage = (float(tabs[3]) / seqlengths[query]) * 100
    if perc_ident > minpercidentity and (perc_coverage > minseqcoverage or alignmentlength > 40):
      blastlines2.append(i)
  blastlines = blastlines2
  #Goes through the blastlines. For each query, creates a querydict and hitlist, and adds these to the blastdict when finding the next query
  firstquery = "y"
  for i in blastlines:
    tabs = i.split("\t")
    query = tabs[0]
    
    second_column_split = tabs[1].split("|")
    
    subject = second_column_split[4]
    if subject == "no_locus_tag":
      subject = second_column_split[6]
    if subject in geneclustergenes:
      subject = "h_" + subject
    if len(second_column_split) > 6:
      locustag = second_column_split[6]
    else:
      locustag = ""
    subject_genecluster = second_column_split[0] + "_" + second_column_split[1]
    subject_start = (second_column_split[2]).split("-")[0]
    subject_end = (second_column_split[2]).split("-")[1]
    subject_strand  = second_column_split[3]
    subject_annotation = second_column_split[5]
    perc_ident = int(tabs[2].split(".")[0])
    alignmentlength = float(tabs[3])
    evalue = str(tabs[10])
    blastscore = int(tabs[11].split(".", 1)[0])
    if seqlengths.has_key(query):
      perc_coverage = (float(tabs[3]) / seqlengths[query]) * 100
    else:
      seqlength = len(seqdict[query.split("|")[4]])
      perc_coverage = (float(tabs[3]) / seqlength) * 100
    if firstquery == "y": #Only until the first blastline with good hit
      firstquery = "n"
      querylist.append(query)
      subjectlist = []
      querydict = {}
      subjectlist.append(subject)
      querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
      if subject_genecluster not in hitclusters:
        hitclusters.append(subject_genecluster)
      last_query = query
    elif i == blastlines[-1]: #Only for the last blastline
      if query not in querylist:
        subjectlist.append(subject)
        querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
        blastdict[query] = [subjectlist,querydict]
        querylist.append(query)
        if subject_genecluster not in hitclusters:
          hitclusters.append(subject_genecluster)
      else:
        subjectlist.append(subject)
        querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
        blastdict[query] = [subjectlist,querydict]
    else: #For all but the first and last blastlines
      if query not in querylist:
        blastdict[last_query] = [subjectlist,querydict]
        querylist.append(query)
        subjectlist = []
        querydict = {}
        subjectlist.append(subject)
        querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
        if subject_genecluster not in hitclusters:
          hitclusters.append(subject_genecluster)
        last_query = query
      else:
        subjectlist.append(subject)
        querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
        if subject_genecluster not in hitclusters:
          hitclusters.append(subject_genecluster)
  return [blastdict,querylist,hitclusters]

def getdircontents():
  return os.listdir('.')
  """
  if sys.platform == ('win32'):
    dircontents = os.popen("dir/w")
    dircontents = dircontents.read()
    dircontents = dircontents.replace("\n"," ")
    dircontents = dircontents.split(" ")
  if sys.platform == ('linux2'):
    dircontents = os.popen("ls")
    dircontents = dircontents.read()
    dircontents = dircontents.replace("\n"," ")
    dircontents = dircontents.replace("\r"," ")
    dircontents = dircontents.split(" ")
  
  return dircontents
  """

def _gene_arrow(start,end,strand,color,base,height):
    halfheight = height/2
    if start > end:
      start2 = end
      end2 = start
      start = start2
      end = end2
    dist = 100
    oh = ShapeBuilder()
    if (end - start) < halfheight:
        if (strand == "+"):
            pointsAsTuples=[(start,base),
                            (end,base - halfheight),
                            (start,base - height),
                            (start,base)
    ]
        elif (strand == "-"):
            pointsAsTuples=[(start,base - halfheight),
                            (end,base - height),
                            (end,base),
                            (start,base - halfheight)
                            ]
    else:
        if (strand == "+"):
            arrowstart = end-halfheight
            pointsAsTuples=[(start,base),
                            (arrowstart,base),
                            (end,base-halfheight),
                            (arrowstart,base - height),
                            (start,base - height),
                            (start,base)
                            ]
        elif (strand == "-"):
            arrowstart = start + halfheight
            pointsAsTuples=[(start,base - halfheight),
                            (arrowstart,base - height),
                            (end,base - height),
                            (end,base),
                            (arrowstart,base),
                            (start,base - halfheight)
                            ]
    pg=oh.createPolygon(points=oh.convertTupleArrayToPoints(pointsAsTuples),strokewidth=1, stroke='black', fill=color)
    return pg        

def _gene_label(start,end,name,y,screenwidth):
    #Add gene label
    txt = name
    myStyle = StyleBuilder()
    myStyle.setFontFamily(fontfamily="Verdana")
    #myStyle.setFontWeight(fontweight='bold')
    myStyle.setFontStyle(fontstyle='italic')
    myStyle.setFontSize('10px')
    myStyle.setFilling('#600000')
    x =  ((start + end)/2)
    base = 35
    height = 10
    halfheight = height/2
    y =  base + halfheight
    t1 = text(txt,x,y)
    t1.set_style(myStyle.getStyle())
    return t1

def relativepositions(starts,ends,largestclustersize):
  rel_starts = []
  rel_ends = []
  #Assign relative start and end sites for visualization
  lowest_start = int(starts[0])
  leftboundary = lowest_start
  for i in starts:
    i = float(float(int(i) - int(leftboundary)) / largestclustersize) * screenwidth * 0.75
    i = int(i)
    rel_starts.append(i)
  for i in ends:
    i = float(float(int(i) - int(leftboundary)) / largestclustersize) * screenwidth * 0.75
    i = int(i)
    rel_ends.append(i)
  return [rel_starts,rel_ends]

def startendsitescheck(starts,ends):
  #Check whether start sites are always lower than end sites, reverse if necessary
  starts2 = []
  ends2 = []
  a = 0
  for i in starts:
    if int(i) > int(ends[a]):
      starts2.append(ends[a])
      ends2.append(i)
    else:
      starts2.append(i)
      ends2.append(ends[a])
    a += 1
  ends = ends2
  starts = starts2
  return [starts,ends]

def RadialGradient(startcolor,stopcolor,gradientname):
  d = defs()
  rg = radialGradient()
  rg.set_id(gradientname)
  s = stop(offset="0%")
  s.set_stop_color(startcolor)
  s.set_stop_opacity(1)
  rg.addElement(s)
  s = stop(offset="100%")
  s.set_stop_color(stopcolor)
  s.set_stop_opacity(1)
  rg.addElement(s)
  d.addElement(rg)
  return d
  
def LinearGradient(startcolor,stopcolor,gradientname):
  d = defs()
  lg = linearGradient()
  lg.set_id(gradientname)
  s = stop(offset="0%")
  s.set_stop_color(startcolor)
  s.set_stop_opacity(1)
  lg.addElement(s)
  s = stop(offset="100%")
  s.set_stop_color(stopcolor)
  s.set_stop_opacity(1)
  lg.addElement(s)
  d.addElement(lg)
  return d
  
def generate_rgbscheme(nr):
  usablenumbers = [1,2,4,8,12,18,24,32,48,64,10000]
  lengthsdict = {1:[1,1,1],2:[1,1,2],4:[1,2,2],8:[2,2,2],12:[2,2,3],18:[2,3,3],24:[3,3,3],32:[3,3,4],48:[3,4,4],64:[4,4,4]}
  shortestdistance = 10000
  for i in usablenumbers:
    distance = i - nr
    if distance >= 0:
      if distance < shortestdistance:
        shortestdistance = distance
        closestnr = i
  toohigh = "n"
  if closestnr == 10000:
    toohigh = "y"
    closestnr = 64
  xyznumbers = lengthsdict[closestnr]
  x = xyznumbers[0]
  y = xyznumbers[1]
  z = xyznumbers[2]
  xpoints = []
  xpoint = (255/z)/2
  for i in range(x):
    xpoints.append(xpoint)
    xpoint += (255/x)
  ypoints = []
  ypoint = (255/z)/2
  for i in range(y):
    ypoints.append(ypoint)
    ypoint += (255/y)
  zpoints = []
  zpoint = (255/z)/2
  for i in range(z):
    zpoints.append(zpoint)
    zpoint += (255/z)
  colorlist = []
  for i in xpoints:
    for j in ypoints:
      #for k in zpoints:
      #  rgb = "rgb(%s,%s,%s)" % (i, j, k)
      #  #rgb = "rgb(" + str(i) + "," + str(j) + "," + str(k) + ")"
      #  colorlist.append(rgb)
      [colorlist.append("rgb(%s,%s,%s)" % (i, j, k)) for k in zpoints]
  if toohigh == "y":
    colorlist = colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist
  if closestnr == 24:
    colorlist = colorlist[:15] + colorlist[18:]
  if closestnr == 32:
    colorlist = colorlist[:21] + colorlist[24:]
  colorlist2 = []
  if closestnr == 1:
    colorlist2.append("red")
  if closestnr == 2:
    colorlist2.append("red")
    colorlist2.append("green")
  if closestnr == 4:
    colorlist2.append("red")
    colorlist2.append("green")
    colorlist2.append("blue")
    colorlist2.append("yellow")
  if closestnr == 8:
    neworder=[4,1,2,5,6,7,3,0]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 12:
    neworder=[6,3,5,9,7,2,11,4,8,1,10,0]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 18:
    neworder=[9,6,2,14,15,8,12,10,3,5,7,11,4,1,16,13,0]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 24:
    neworder=[15,12,9,6,5,0,21,1,16,14,8,17,2,23,22,3,13,7,10,4,18,20,19,11]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 32:
    neworder = [21,19,27,6,8,1,14,7,20,13,9,30,4,23,18,12,5,29,24,17,11,31,2,28,22,15,26,3,20,16,10,25]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr > 32:
    random.shuffle(colorlist)
    colorlist2 = colorlist
  colorlist = colorlist2
  return colorlist

def geneclustersvg(genes,rel_starts,rel_ends,strands,geneposdict,pksnrpsprots,pksnrpsdomains,qclusternr):
  nrgenes = len(genes)
  #Define relative start and end positions for plotting
  s = svg(x = 0, y = 0, width = (screenwidth * 0.75), height = (259 + 99 * len(pksnrpsprots)))
  viewbox = "0 -30 " + str(screenwidth * 0.8) + " " + str(185 + 70 * len(pksnrpsprots))
  s.set_viewBox(viewbox)
  s.set_preserveAspectRatio("none")
  
  #Add line behind gene arrows
  oh = ShapeBuilder()
  group = g()
  group.addElement(oh.createLine(10,60,10 + (screenwidth * 0.75),60, strokewidth = 2, stroke = "grey"))
  s.addElement(group)
  #Add gene arrows
  a = 0
  y = 0
  for x in range(nrgenes):
    group = g()       
    #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
    group.addElement(_gene_arrow(10 + rel_starts[a],10 + rel_ends[a],strands[a],colors[a],65,10))
  #Can be used for domains
  #   group.addElement(oh.createRect(rel_starts[a],45,(rel_ends[a]-rel_starts[a]),10, strokewidth = 2, stroke = "black", fill="#237845"))
    group.set_id("a" + str(qclusternr) + "_00%s"%x)
    s.addElement(group)
    if y == 0:
      y = 1
    elif y == 1:
      y = 0
    a += 1
  #Add domain depictions
  oh = ShapeBuilder()
  group = g()
  #Determine longest protein to decide on scaling
  longestprot = 0
  protlengthdict = {}
  for i in pksnrpsprots:
    protlength = (geneposdict[i][1] - geneposdict[i][0]) / 3
    protlengthdict[i] = protlength
    if protlength > longestprot:
      longestprot = protlength
  z = 1
  w = 0
  ksnr = 1
  atnr = 1
  dhnr = 1
  krnr = 1
  ernr = 1
  acpnr = 1
  cnr = 1
  enr = 1
  anr = 1
  pcpnr = 1
  tenr = 1
  othernr = 1 
  for i in pksnrpsprots:
    domains = pksnrpsdomains[i][0]
    domainsdict = pksnrpsdomains[i][1]
    protlength = protlengthdict[i]
    group.addElement(oh.createLine(10,(125 + z * 60 ),10 + ((float(protlength) / float(longestprot)) * (screenwidth * 0.75)),(125 + z * 60 ), strokewidth = 1, stroke = "grey"))
    s.addElement(group)
    try:
        aa2pixelratio = longestprot * 0.75 / screenwidth
    except:
        aa2pixelratio = 0.1
    #print 'logestprot', longestprot
    #print 'scrennwidth', screenwidth
    #print aa2pixelratio
    myStyle = StyleBuilder()
    myStyle.setFontFamily(fontfamily="MS Reference Sans Serif")
    myStyle.setFontWeight(fontweight='bold')
    myStyle.setFontSize('12px')
    for j in domains:
      startpos = domainsdict[j][0]
      endpos = domainsdict[j][1]
      if "PKS_KS" in j:
        c = LinearGradient("#08B208","#81F781","KS_domain"+str(qclusternr) + "_" + str(ksnr))
        d = LinearGradient("#81F781","#08B208","KS_line"+str(qclusternr) + "_" + str(ksnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#KS_line' + str(qclusternr) + "_" + str(ksnr) + ")",fill="url(#KS_domain" + str(qclusternr) + "_" + str(ksnr) + ")")
        f = text("KS",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A2A0A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("KS",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#3B0B0B')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        ksnr += 1
      elif "PKS_AT" in j:
        c = LinearGradient("#DC0404","#F78181","AT_domain"+str(qclusternr) + "_" + str(atnr))
        d = LinearGradient("#F78181","#DC0404","AT_line"+str(qclusternr) + "_" + str(atnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#AT_line' + str(qclusternr) + "_" + str(atnr) + ")",fill="url(#AT_domain" + str(qclusternr) + "_" + str(atnr) + ")")
        f = text("AT",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#2A1B0A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("AT",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#2A1B0A')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        atnr += 1
      elif "PKS_DH" in j:
        c = LinearGradient("#B45F04","#F7BE81","DH_domain"+str(qclusternr) + "_" + str(dhnr))
        d = LinearGradient("#F7BE81","#B45F04","DH_line"+str(qclusternr) + "_" + str(dhnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#DH_line' + str(qclusternr) + "_" + str(dhnr) + ")",fill="url(#DH_domain" + str(qclusternr) + "_" + str(dhnr) + ")")
        f = text("DH",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#3B0B0B')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("DH",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#3B0B0B')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        dhnr += 1
      elif "PKS_KR" in j:
        c = LinearGradient("#089E4B","#81F781","KR_domain"+str(qclusternr) + "_" + str(krnr))
        d = LinearGradient("#81F781","#089E4B","KR_line"+str(qclusternr) + "_" + str(krnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#KR_line' + str(qclusternr) + "_" + str(krnr) + ")",fill="url(#KR_domain" + str(qclusternr) + "_" + str(krnr) + ")")
        f = text("KR",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A2A1B')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("KR",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0A2A1B')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        krnr += 1
      elif "PKS_ER" in j:
        c = LinearGradient("#089E85","#81F7F3","ER_domain"+str(qclusternr) + "_" + str(ernr))
        d = LinearGradient("#81F7F3","#089E85","ER_line"+str(qclusternr) + "_" + str(ernr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#ER_line' + str(qclusternr) + "_" + str(ernr) + ")",fill="url(#ER_domain" + str(qclusternr) + "_" + str(ernr) + ")")
        f = text("ER",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A2A29')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("ER",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0A2A29')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        ernr += 1
      elif "ACP" in j:
        c = LinearGradient("#084BC6","#81BEF7","ACP_domain"+str(qclusternr) + "_" + str(acpnr))
        d = LinearGradient("#81BEF7","#084BC6","ACP_line"+str(qclusternr) + "_" + str(acpnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#ACP_line' + str(qclusternr) + "_" + str(acpnr) + ")",fill="url(#ACP_domain" + str(qclusternr) + "_" + str(acpnr) + ")")
        f = text("ACP",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A1B2A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("ACP",((-2 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0A1B2A')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        acpnr += 1
      elif ("C" in j or "Heterocyclization" in j ) and "ACP" not in j and "PCP" not in j and "NRPS-COM" not in j and "CAL" not in j:
        c = LinearGradient("#393989","#8181F7","C_domain"+str(qclusternr) + "_" + str(cnr))
        d = LinearGradient("#8181F7","#393989","C_line"+str(qclusternr) + "_" + str(cnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#C_line' + str(qclusternr) + "_" + str(cnr) + ")",fill="url(#C_domain" + str(qclusternr) + "_" + str(cnr) + ")")
        f = text("C",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A0A2A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("C",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0A0A2A')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        cnr += 1
      elif "Epimerization" in j and "ER" not in j and "TE" not in j:
        c = LinearGradient("#393989","#8181F7","E_domain"+str(qclusternr) + "_" + str(enr))
        d = LinearGradient("#8181F7","#393989","E_line"+str(qclusternr) + "_" + str(enr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#E_line' + str(qclusternr) + "_" + str(enr) + ")",fill="url(#E_domain" + str(qclusternr) + "_" + str(enr) + ")")
        f = text("E",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A0A2A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("E",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0A0A2A')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        enr += 1
      elif ("AMP" in j or "A-OX" in j):
        c = LinearGradient("#56157F","#BE81F7","A_domain"+str(qclusternr) + "_" + str(anr))
        d = LinearGradient("#BE81F7","#56157F","A_line"+str(qclusternr) + "_" + str(anr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#A_line' + str(qclusternr) + "_" + str(anr) + ")",fill="url(#A_domain" + str(qclusternr) + "_" + str(anr) + ")")
        f = text("A",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#1B0A2A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("A",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#1B0A2A')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        anr += 1
      elif "PCP" in j:
        c = LinearGradient("#084BC6","#81BEF7","PCP_domain"+str(qclusternr) + "_" + str(pcpnr))
        d = LinearGradient("#81BEF7","#084BC6","PCP_line"+str(qclusternr) + "_" + str(pcpnr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#PCP_line' + str(qclusternr) + "_" + str(pcpnr) + ")",fill="url(#PCP_domain" + str(qclusternr) + "_" + str(pcpnr) + ")")
        f = text("PCP",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0A1B2A')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          f = text("PCP",((-2 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0A1B2A')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        pcpnr += 1
      elif "Thioesterase" in j or "TD" in j:
        c = LinearGradient("#750072","#F5A9F2","TE_domain"+str(qclusternr) + "_" + str(tenr))
        d = LinearGradient("#F5A9F2","#750072","TE_line"+str(qclusternr) + "_" + str(tenr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#TE_line' + str(qclusternr) + "_" + str(tenr) + ")",fill="url(#TE_domain" + str(qclusternr) + "_" + str(tenr) + ")")
        if "Thioesterase" in j:
          f = text("TE",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#2A0A29')
        else:
          f = text("TD",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#2A0A29')
        if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
          myStyle.setFontSize('8px')
          if "Thioesterase" in j:
            f = text("TE",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#2A0A29')
          else:
            f = text("TD",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#2A0A29')
        elif ((endpos-startpos) / aa2pixelratio) < 20:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        tenr += 1
      else:
        c = LinearGradient("#929292","#DBDBDB","other_domain"+str(qclusternr) + "_" + str(othernr))
        d = LinearGradient("#DBDBDB","#929292","other_line"+str(qclusternr) + "_" + str(othernr))
        e = oh.createRect(str(10 + startpos / aa2pixelratio),str((125 + z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#other_line' + str(qclusternr) + "_" + str(othernr) + ")",fill="url(#other_domain" + str(qclusternr) + "_" + str(othernr) + ")")
        domname = (((((((((j.replace("0","")).replace("1","")).replace("2","")).replace("3","")).replace("4","")).replace("5","")).replace("6","")).replace("7","")).replace("8","")).replace("9","")
        if len(domname) == 1:
          f = text(domname,((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0B0B0B')
        elif len(domname) == 2:
          f = text(domname,((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0B0B0B')
        elif len(domname) == 3:
          f = text(domname,((-12 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 4),fill='#0B0B0B')
        if len(domname) > 3 or ((endpos-startpos) / aa2pixelratio) < 100:
          myStyle.setFontSize('8px')
          f = text(domname,((-16 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0B0B0B')
        if len(domname) > 4 and ((endpos-startpos) / aa2pixelratio) < 100:
          myStyle.setFontSize('6px')
          f = text(domname,((-16 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((125 + z * 60 ) + 3),fill='#0B0B0B')
        if ((endpos-startpos) / aa2pixelratio) < 60:
          f = "notext"
        if f != "notext":
          f.set_style(myStyle.getStyle())
        myStyle.setFontSize('12px')
        group = g()
        group.addElement(c)
        group.addElement(d)
        group.addElement(e)
        if f != "notext":
          group.addElement(f)
        group.set_id("b" + str(qclusternr) + "_00%s"%w)
        s.addElement(group)
        othernr += 1
      w += 1
    z += 1
  s.addElement(group)
  return s
  
def calculate_colorgroups(queryclusternumber,hitclusternumbers,queryclusterdata,internalhomologygroupsdict):
  #Extract data and generate color scheme
  nrhitclusters = queryclusterdata[queryclusternumber][0]
  hitclusterdata = queryclusterdata[queryclusternumber][1]
  queryclustergenes = hitclusterdata[1][3]
  queryclustergenesdetails = hitclusterdata[1][4]
  colorgroupsdict = {}
  colorgroupslengthlist = []
  colorgroupslist = []
  for hitclusternumber in hitclusternumbers:
    colorgroups = hitclusterdata[hitclusternumber][0][hitclusternumber]
    colorgroupsdict[hitclusternumber] = colorgroups
    colorgroupslengthlist.append(len(colorgroups))
    colorgroupslist.append(colorgroups)
  metacolorgroups = []
  internalgroups = internalhomologygroupsdict[queryclusternumber]
  for i in internalgroups:
    metagroup = []
    for j in i:
      for m in colorgroupslist:
        for l in m:
          if j in l:
            #for k in l:
            #  if k not in metagroup:
            #    metagroup.append(k)
            [metagroup.append(k) for k in l if k not in metagroup]
    if len(metagroup) > 1 and metagroup not in metacolorgroups:
      metacolorgroups.append(metagroup)
  #Generate RGB scheme
  rgbcolorscheme = generate_rgbscheme(len(metacolorgroups))
  rgbcolorscheme.append("#FFFFFF")
  #Create colorschemedict in which all genes that are hits of the same query gene get the same color
  colorschemedict = {}
  z = 0
  for i in queryclustergenes:
    for j in metacolorgroups:
      if i in j:
        for l in j:
          if colorschemedict.has_key(l):
            pass
          else:
            colorschemedict[l] = z
        #[colorschemedict[l] = z for l in j if not coloschemedict.has_key(l)]
    if z in colorschemedict.values():
      z += 1
  return colorschemedict,rgbcolorscheme
 
def clusterblastresults(queryclusternumber,hitclusternumbers,queryclusterdata,colorschemedict,rgbcolorscheme):
  #print "Generating svg for cluster",queryclusternumber
  #Extract data and generate color scheme
  nrhitclusters = queryclusterdata[queryclusternumber][0]
  hitclusterdata = queryclusterdata[queryclusternumber][1]
  queryclustergenes = hitclusterdata[1][3]
  queryclustergenesdetails = hitclusterdata[1][4]
  colorgroupsdict = {}
  colorgroupslengthlist = []
  colorgroupslist = []
  for hitclusternumber in hitclusternumbers:
    colorgroups = hitclusterdata[hitclusternumber][0][hitclusternumber]
    colorgroupsdict[hitclusternumber] = colorgroups
    colorgroupslengthlist.append(len(colorgroups))
    colorgroupslist.append(colorgroups)
  #Find out whether hit gene cluster needs to be inverted compared to query gene cluster
  strandsbalancedict = {}
  for m in hitclusternumbers:
    hitclustergenesdetails = hitclusterdata[m][2]
    strandsbalance = 0
    for i in queryclustergenes:
      refstrand = queryclustergenesdetails[i][2]
      for j in colorgroupsdict[m]:
        if i in j:
          for k in j:
            if k in hitclusterdata[m][1] and hitclustergenesdetails[k][2] == refstrand:
              strandsbalance += 1
            elif k in hitclusterdata[m][1] and hitclusterdata[m][2][k][2] != refstrand:
              strandsbalance = strandsbalance - 1
    strandsbalancedict[m] = strandsbalance
  #Generate coordinates for SVG figure
  qnrgenes = len(queryclustergenes)
  qstarts =[]
  qends = []
  qstrands =[]
  qcolors = []
  for i in queryclustergenes:
    qgenedata = queryclustergenesdetails[i]
    if qgenedata[0] > qgenedata[1]:
      qstarts.append(qgenedata[0])
      qends.append(qgenedata[1])
    else:
      qstarts.append(qgenedata[1])
      qends.append(qgenedata[0])
    qstrands.append(qgenedata[2])
    if colorschemedict.has_key(i):
      qcolors.append(colorschemedict[i])
    else:
      qcolors.append("white")
  qstarts_ends = startendsitescheck(qstarts,qends)
  qstarts = qstarts_ends[0]
  qends = qstarts_ends[1]
  hdata = {}
  for m in hitclusternumbers:
    hitclustergenes = hitclusterdata[m][1]
    hitclustergenesdetails = hitclusterdata[m][2]
    hnrgenes = len(hitclustergenes)
    hstarts =[]
    hends = []
    hstrands =[]
    hcolors = []
    for i in hitclustergenes:
      hgenedata = hitclustergenesdetails[i]
      if hgenedata[0] > hgenedata[1]:
        hstarts.append(hgenedata[0])
        hends.append(hgenedata[1])
      else:
        hstarts.append(hgenedata[1])
        hends.append(hgenedata[0])
      hstrands.append(hgenedata[2])
      if colorschemedict.has_key(i):
        hcolors.append(colorschemedict[i])
      else:
        hcolors.append("white")
    #Invert gene cluster if needed
    if strandsbalancedict[m] < 0:
      hstarts2 = []
      hends2 = []
      hstrands2 = []
      for i in hstarts:
        hstarts2.append(str(100000000 - int(i)))
      hstarts = hstarts2
      hstarts.reverse()
      for i in hends:
        hends2.append(str(100000000 - int(i)))
      hends = hends2
      hends.reverse()
      for i in hstrands:
        if i == "+":
          hstrands2.append("-")
        elif i == "-":
          hstrands2.append("+")
      hstrands = hstrands2
      hstrands.reverse()
      hcolors.reverse()
    hstarts_ends = startendsitescheck(hstarts,hends)
    hstarts = hstarts_ends[0]
    hends = hstarts_ends[1]
    hdata[m] = [hstarts,hends,hstrands,hcolors]
  #Find cluster size of largest cluster of query & all hit clusters assessed
  clustersizes = []
  for m in hitclusternumbers:
    hclustersize = int(hdata[m][1][-1]) - int(hdata[m][0][0])
    clustersizes.append(hclustersize)
  qclustersize = int(qends[-1]) - int(qstarts[0])
  clustersizes.append(qclustersize)
  largestclustersize = max(clustersizes) 
  smallestclustersize = min(clustersizes)
  #Find relative positions
  qrelpositions = relativepositions(qstarts,qends,largestclustersize)
  qrel_starts = qrelpositions[0]
  qrel_ends = qrelpositions[1]
  qdata = [qrel_starts,qrel_ends,qstrands,qcolors]
  hdata2 = {}
  qdata2 = []
  for m in hitclusternumbers:
    hclustersize = int(hdata[m][1][-1]) - int(hdata[m][0][0])
    hrelpositions = relativepositions(hdata[m][0],hdata[m][1],largestclustersize)
    hrel_starts = hrelpositions[0]
    hrel_ends = hrelpositions[1]
    #Center-align smallest gene cluster
    if largestclustersize == hclustersize:
      qrel_ends2 = []
      qrel_starts2 = []
      for i in qrel_starts:
        qrel_starts2.append(int(i) + int(float(float((largestclustersize - qclustersize) / 2) / largestclustersize) * screenwidth * 0.75))
      for i in qrel_ends:
        qrel_ends2.append(int(i) + int(float(float((largestclustersize - qclustersize) / 2) / largestclustersize) * screenwidth * 0.75))
      qrel_ends = qrel_ends2
      qrel_starts = qrel_starts2
    else:
      hrel_ends2 = []
      hrel_starts2 = []
      for i in hrel_starts:
        hrel_starts2.append(int(i) + int(float(float((largestclustersize - hclustersize) / 2) / largestclustersize) * screenwidth * 0.75))
      for i in hrel_ends:
        hrel_ends2.append(int(i) + int(float(float((largestclustersize - hclustersize) / 2) / largestclustersize) * screenwidth * 0.75))
      hrel_ends = hrel_ends2
      hrel_starts = hrel_starts2
    hdata2[m] = [hrel_starts,hrel_ends,hdata[m][2],hdata[m][3]]
    qdata2 = [qrel_starts,qrel_ends,qdata[2],qdata[3]]
  hdata = hdata2
  qdata = qdata2
  s = svg(x = 0, y = 0, width = (screenwidth * 0.75), height = (270 + len(hitclusternumbers) * 50))
  viewbox = "0 0 " + str(screenwidth * 0.8) + " " + str(180 + len(hitclusternumbers) * 50)
  s.set_viewBox(viewbox)
  s.set_preserveAspectRatio("none")
  #Add line behind query gene cluster gene arrows
  oh = ShapeBuilder()
  group = g()
  group.addElement(oh.createLine(10,35,10 + (screenwidth * 0.75),35, strokewidth = 1, stroke = "grey"))
  s.addElement(group)
  #Add query gene cluster gene arrows
  a = 0
  y = 0
  for x in range(qnrgenes):
    group = g()
    #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
    if qcolors[a] == "white":
      group.addElement(_gene_arrow(10 + qrel_starts[a],10 + qrel_ends[a],qstrands[a],rgbcolorscheme[-1],40,10))
    else:
      group.addElement(_gene_arrow(10 + qrel_starts[a],10 + qrel_ends[a],qstrands[a],rgbcolorscheme[qcolors[a]],40,10))
    #Can be used for domains
    #group.addElement(oh.createRect(rel_starts[a],45,(rel_ends[a]-rel_starts[a]),10, strokewidth = 2, stroke = "black", fill="#237845"))
    if len(hitclusternumbers) == 1:
      group.set_id("q" + str(queryclusternumber) + "_" + str(hitclusternumbers[0]) + "_" + "%s"%x)
    else:
      group.set_id("all_" + str(queryclusternumber) + "_0_" + "%s"%x)
    s.addElement(group)
    if y == 0:
      y = 1
    elif y == 1:
      y = 0
    a += 1
  for m in hitclusternumbers:
    #Add line behind hit gene cluster gene arrows
    group.addElement(oh.createLine(10,35 + 50 * (hitclusternumbers.index(m) + 1),10 + (screenwidth * 0.75),35 + 50 * (hitclusternumbers.index(m) + 1), strokewidth = 1, stroke = "grey"))
    s.addElement(group)
    #Add hit gene cluster gene arrows
    hitclustergenes = hitclusterdata[m][1]
    hnrgenes = len(hitclustergenes)
    hrel_starts = hdata[m][0]
    hrel_ends = hdata[m][1]
    hstrands = hdata[m][2]
    hcolors = hdata[m][3]
    a = 0
    y = 0
    for x in range(hnrgenes):
      group = g()       
      #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
      if hcolors[a] == "white":
        group.addElement(_gene_arrow(10 + hrel_starts[a],10 + hrel_ends[a],hstrands[a],rgbcolorscheme[-1],40 + 50 * (hitclusternumbers.index(m) + 1),10))
      else:
        group.addElement(_gene_arrow(10 + hrel_starts[a],10 + hrel_ends[a],hstrands[a],rgbcolorscheme[hcolors[a]],40 + 50 * (hitclusternumbers.index(m) + 1),10))
      #Can be used for domains
      #   group.addElement(oh.createRect(rel_starts[a],45,(rel_ends[a]-rel_starts[a]),10, strokewidth = 2, stroke = "black", fill="#237845"))
      if len(hitclusternumbers) == 1:
        group.set_id("h" + str(queryclusternumber) + "_" + str(m) + "_" + "%s"%x)
      else:
        group.set_id("all_" + str(queryclusternumber) + "_" + str(m) + "_" + "%s"%x)
      s.addElement(group)
      if y == 0:
        y = 1
      elif y == 1:
        y = 0
      a += 1
  return [s,[qdata,hdata,strandsbalancedict]]

def runblast(query):
  blastsearch = "blastp  -db "+antismash_path+"clusterblast/geneclusterprots.fasta -query " + query + " -outfmt 6 -max_target_seqs 1000 -evalue 1e-05 -out " + query.split(".")[0] + ".out"
  os.system(blastsearch)

def smcog_analysis(inputgenes,inputnr,accessiondict,seqdict,smcogdict,smcogsoutputfolder):
  #create input.fasta file with single query sequence to be used as input for MSA
  for k in inputgenes:
    gene = accessiondict[k]
    tag = k
    seq = seqdict[k]
    writefasta([tag],[seq],"input" + str(inputnr) + ".fasta")
    if len(smcogdict[k]) > 0:
      smcog = (smcogdict[k][0][0]).split(":")[0]
      #Align to multiple sequence alignment, output as fasta file
      fastafile = "input" + str(inputnr) + ".fasta"
      musclecommand = "muscle -quiet -profile -in1 " + str(smcog).lower() + "_muscle.fasta -in2 input" + str(inputnr) + ".fasta -out muscle" + str(inputnr) + ".fasta"
      os.system(musclecommand)
      #Trim alignment 
      #edit muscle fasta file: remove all positions before the first and after the last position shared by >33% of all sequences
      file = open("muscle" + str(inputnr) + ".fasta","r")
      filetext = file.read()
      filetext = filetext.replace("\r","\n")
      lines = filetext.split("\n")
      ##Combine all sequence lines into single lines
      lines2 = []
      seq = ""
      nrlines = len(lines)
      a = 0
      lines = lines[:-1]
      for i in lines:
        if a == (nrlines - 2):
          seq = seq + i
          lines2.append(seq)
        if i[0] == ">":
          lines2.append(seq)
          seq = ""
          lines2.append(i)
        else:
          seq = seq + i
        a += 1
      lines = lines2[1:]
      #Retrieve names and seqs from muscle fasta lines
      seqs = []
      names = []
      for i in lines:
        if len(i) > 0 and i[0] == ">":
          name = i[1:]
          names.append(name)
        else:
          seq = i
          seqs.append(seq)
      #Find first and last amino acids shared conserved >33%
      #Create list system to store conservation of residues
      conservationlist = []
      lenseqs = len(seqs[0])
      nrseqs = len(seqs)
      for i in range(lenseqs):
        conservationlist.append({"A":0,"B":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"J":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"U":0,"V":0,"W":0,"X":0,"Y":0,"Z":0,"-":0})
      a = 0
      for i in seqs:
        aa = list(i)
        for i in aa:
          conservationlist[a][i] += 1
          a += 1
        a = 0
      firstsharedaa = 0
      lastsharedaa = lenseqs
      #Find first amino acid shared
      first = "yes"
      nr = 0
      for i in conservationlist:
        aa = sortdictkeysbyvaluesrev(i)
        if aa[0] != "-" and i[aa[1]] > (nrseqs / 3) and first == "yes":
          firstsharedaa = nr
          first = "no"
        nr += 1
      #Find last amino acid shared
      conservationlist.reverse()
      first = "yes"
      nr = 0
      for i in conservationlist:
        aa = sortdictkeysbyvaluesrev(i)
        if aa[0] != "-" and i[aa[1]] > (nrseqs / 3) and first == "yes":
          lastsharedaa = lenseqs - nr
          first = "no"
        nr += 1
      #Shorten sequences to detected conserved regions
      seqs2 = []
      for i in seqs:
        seq = i[firstsharedaa:lastsharedaa]
        seqs2.append(seq)
      seqs = seqs2
      seedfastaname = "trimmed_alignment" + str(inputnr) + ".fasta"
      writefasta(names,seqs,seedfastaname)
      #Draw phylogenetic tree with fasttree 2.1.1
      nwkfile = "tree" + str(inputnr) + ".nwk"
      if sys.platform == ('win32'):
        fasttreecommand = "fasttree -quiet -fastest -noml trimmed_alignment" + str(inputnr) + ".fasta > " + nwkfile
      elif sys.platform == ('linux2'):
        fasttreecommand = "./FastTree -quiet -fastest -noml trimmed_alignment" + str(inputnr) + ".fasta > " + nwkfile
      os.system(fasttreecommand)
      #Convert tree to XTG and draw PNG image using TreeGraph
      p = subprocess.Popen("java -Djava.awt.headless=true -jar TreeGraph.jar -convert tree" + str(inputnr) + ".nwk -xtg tree" + str(inputnr) + ".xtg", shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
      processes_starttime = time.time()
      while True:
        if (time.time() - processes_starttime) > 300:
          if sys.platform == ('linux2'):
            os.kill(p.pid,signal.SIGKILL)
            break
          if sys.platform == ('win32'):
            subprocess.Popen("taskkill /F /T /PID %i"%p.pid , shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            break
        if p.poll() == 0:
          break
        time.sleep(2)
      out, err = p.communicate()
      output = out
      if "exception" not in output and "Exception" not in output:
        p = subprocess.Popen("java -Djava.awt.headless=true -jar TreeGraph.jar -image tree" + str(inputnr) + ".xtg " + tag.split(".")[0] + ".png", shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        processes_starttime = time.time()
        while True:
          if (time.time() - processes_starttime) > 300:
            if sys.platform == ('linux2'):
              os.kill(p.pid,signal.SIGKILL)
              break
            if sys.platform == ('win32'):
              subprocess.Popen("taskkill /F /T /PID %i"%p.pid , shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
              break
          if p.poll() == 0:
            break
          time.sleep(2)
        out, err = p.communicate()
        output = out
        if "exception" not in output and "Exception" not in output:
          if sys.platform == ('win32'):
            copycommand = 'copy/y ' + tag.split(".")[0] + '.png "..\\' + smcogsoutputfolder + '" > nul'
          elif sys.platform == ('linux2'):
            copycommand = 'cp ' + tag.split(".")[0] + '.png "../' + smcogsoutputfolder + '" > /dev/null'
          os.system(copycommand)
          if sys.platform == ('win32'):
            os.system("del " + tag.split(".")[0] + ".png")
            os.system("del tree" + str(inputnr) + ".xtg")
            os.system("del trimmed_alignment" + str(inputnr) + ".fasta")
          elif sys.platform == ('linux2'):
            os.system("rm " + tag.split(".")[0] + ".png")
            os.system("rm tree" + str(inputnr) + ".xtg")
            os.system("rm trimmed_alignment" + str(inputnr) + ".fasta")

def depict_smile(genecluster,structuresfolder):
  if sys.platform == ('win32'):
    indigo_depict_command1 = "indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + "_icon.png -query -w 200 -h 150"
    indigo_depict_command2 = "indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + ".png -query"
  elif sys.platform == ('linux2'):
    indigo_depict_command1 = "./indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + "_icon.png -query -w 200 -h 150"
    indigo_depict_command2 = "./indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + ".png -query"
  os.system(indigo_depict_command1)
  os.system(indigo_depict_command2)
  dircontents = getdircontents()
  geneclusterstring = "genecluster" + str(genecluster) + ".png"
  if geneclusterstring in dircontents:
    if sys.platform == ('win32'):
      structuresfolder = structuresfolder.replace("/","\\")
      copycommand1 = "copy/y genecluster" + str(genecluster) + ".png ..\\" + structuresfolder + ' > nul'
      copycommand2 = "copy/y genecluster" + str(genecluster) + "_icon.png ..\\" + structuresfolder + ' > nul'
      delcommand1 = "del genecluster" + str(genecluster) + ".png"
      delcommand2 = "del genecluster" + str(genecluster) + "_icon.png"
      delcommand3 = "del genecluster" + str(genecluster) + ".smi"
      os.system(copycommand1)
      os.system(copycommand2)
      os.system(delcommand1)
      os.system(delcommand2)
      os.system(delcommand3)
    if sys.platform == ('linux2'):
      copycommand1 = "cp genecluster" + str(genecluster) + ".png ../" + structuresfolder
      copycommand2 = "cp genecluster" + str(genecluster) + "_icon.png ../" + structuresfolder
      delcommand1 = "rm genecluster" + str(genecluster) + ".png"
      delcommand2 = "rm genecluster" + str(genecluster) + "_icon.png"
      delcommand3 = "rm genecluster" + str(genecluster) + ".smi"
      os.system(copycommand1)
      os.system(copycommand2)
      os.system(delcommand1)
      os.system(delcommand2)
    return "success"
  else:
    return "failed"

##Core script
import os
from os import system
import sys
import multiprocessing
import time
from multiprocessing import Process, freeze_support
import random
import string
import itertools
from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *
from string import ascii_letters
from pyExcelerator import *
from pyExcelerator.Workbook import *
import signal
import subprocess
starttime = time.time()

os.environ['NRPS2BASEDIR'] = os.path.join(os.getcwd(), 'NRPSPredictor2')

#Fix sys.argv input
options = []
for i in sys.argv:
  if i.count('"') > 1:
    j = i.split(' ')
    for k in j:
      if k[0] == '"':
        k = k + '"'
      elif k[-1] == '"':
        k = '"' + k
      options.append(k)
  else:
    options.append(i)
sys.argv = options
#Redirect stdout and stderr if GUI-executed
if "--gui" in sys.argv and len(sys.argv) < (sys.argv.index("--gui") + 2):
  print >> sys.stderr, "Invalid options input: --gui without n or y"
  print "From the command line, input antismash --help for more information."
  logfile = open("antismash.log","w")
  logfile.write("Invalid options input: --gui without n or y\n")
  logfile.close()
  sys.exit(1)
if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
  stdoutfile = open("stdout.txt","w")
  sys.stdout = stdoutfile
  sys.stderr = stdoutfile

if __name__ == '__main__':
  import shutil
  hmmsearch_path = 'hmmsearch'
  hmmscan_path = 'hmmscan'
  antismash_path = '/home/galaxy/bin/antismash-1.1.0/'
  hmms_path = antismash_path + '/hmms/'
  shutil.copytree(antismash_path + '/NRPSPredictor2/', './NRPSPredictor2/')
  shutil.copytree(antismash_path + '/Minowa/', './Minowa/')
  shutil.copytree(antismash_path + '/pkssignatures/', './pkssignatures/')
  shutil.copytree(antismash_path + '/kr_analysis/', './kr_analysis/')
  shutil.copytree(antismash_path + '/docking_analysis/', './docking_analysis/')
  shutil.copytree(antismash_path + '/NRPeditor/', './NRPeditor/')
  shutil.copy(antismash_path + '/search_form.html', './')
  shutil.copy(antismash_path + '/empty.xhtml', './')
  shutil.copytree(antismash_path + '/vis/', './vis/')
  shutil.copytree(antismash_path + '/smcogtree/', './smcogtree/')

  # add freeze support
  freeze_support()
  
  #Open logfile
  logfile = open("antismash.log","w")

  #Identify screen width
  if sys.platform == ('win32'):
    import ctypes
    user32 = ctypes.windll.user32
    screenwidth = user32.GetSystemMetrics(0)
  if sys.platform == ('linux2'):
    screenwidth = 1024
  #  res = os.popen("xrandr | grep \* | cut -d' ' -f4")  ###FOR SERVER USE###
  #  res = res.read()                                    ###FOR SERVER USE###
  #  screenwidth = int(res.split("x")[0])                ###FOR SERVER USE###
  if screenwidth < 1024:
    screenwidth = 1024
  #temporary for testing
  screenwidth = 1024


  #Reads input
  inputinstructions = "antiSMASH 1.1.0 arguments:\n\nUsage: antismash <query fasta/embl/gbk file>  [options]\n\nOptions (x is an integer number, list x,y,z is a list of integer numbers separated by commas):\n\n--gtransl <x>  : GenBank translation table used for Glimmer (only for FASTA inputs, default: 1)\n1.  The Standard Code\n2.  The Vertebrate Mitochondrial Code\n3.  The Yeast Mitochondrial Code\n4.  The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n5.  The Invertebrate Mitochondrial Code\n6.  The Ciliate, Dasycladacean and Hexamita Nuclear Code\n9.  The Echinoderm and Flatworm Mitochondrial Code\n10. The Euplotid Nuclear Code\n11. The Bacterial, Archaeal and Plant Plastid Code\n12. The Alternative Yeast Nuclear Code\n13. The Ascidian Mitochondrial Code\n14. The Alternative Flatworm Mitochondrial Code\n15. Blepharisma Nuclear Code\n16. Chlorophycean Mitochondrial Code\n21. Trematode Mitochondrial Code\n22. Scenedesmus Obliquus Mitochondrial Code\n23. Thraustochytrium Mitochondrial Code\n--genomeconf <l/c>  : Genome configuration used for Glimmer: linear / circular (only for FASTA inputs, default: l)\n--minglength <x>  : Glimmer minimal gene length (range 30-120, only for FASTA inputs, default: 90)\n--taxon <p/e>  : Taxonomy: prokaryotic / eukaryotic (default: p)\n--cores <x>  : Number of parallel CPUs to use for threading (default: all)\n--clusterblast <y/n> : Include ClusterBlast gene cluster comparison analysis (default:y)\n--smcogs <y/n> : Include smCOG analysis for functional prediction of genes (default:y)\n--fullblast <y/n> : Include genome-wide BLAST analysis (default:n)\n--fullhmm <y/n> : Include genome-wide PFAM HMM analysis (default:n)\n--blastdbpath <path> : Specify folder containing CLUSEAN blast database (default:clusean/db)\n--pfamdbpath <path> : Specify folder containing PFAM database (default:clusean/db)\n--geneclustertypes <x,y,z> : Gene cluster types to scan for (default:1):\n1 = all\n2 = type I polyketide synthases\n3 = type II polyketide synthases\n4 = type III polyketide synthases\n5 = nonribosomal peptide synthetases\n6 = terpene synthases\n7 = lantibiotics\n8 = bacteriocins\n9 = beta-lactams\n10 = aminoglycosides / aminocyclitols\n11 = aminocoumarins\n12 = siderophores\n13 = ectoines\n14 = butyrolactones\n15 = indoles\n16 = nucleosides\n17 = phosphoglycolipids\n18 = melanins\n19 = others\n--help  : this help screen\n"
  #Check input file format
  if len(sys.argv) < 2 or len(sys.argv[1]) < 1:
    print >> sys.stderr, "Please supply valid name for input file."
    print "Usage: antismash <query fasta/embl/gbk file>  [options]"
    print "From the command line, input antismash --help for more information."
    logfile.write("Input format error. Please supply valid name for infile.\n")
    logfile.write("Usage: antismash <query fasta/embl/gbk file>  [options]\n")
    logfile.write("From the command line, input antismash --help for more information.\n")
    logfile.close()
    sys.exit(1)
  if sys.argv[1] != "--help":
    if len(sys.argv[1].split(".")) < 2 or (sys.argv[1].split(".")[-1] != "embl" and sys.argv[1].split(".")[-1] != "EMBL" and sys.argv[1].split(".")[-1] != "emb" and sys.argv[1].split(".")[-1] != "EMB" and sys.argv[1].split(".")[-1] != "genbank" and sys.argv[1].split(".")[-1] != "GENBANK" and sys.argv[1].split(".")[-1] != "gbk" and sys.argv[1].split(".")[-1] != "GBK" and sys.argv[1].split(".")[-1] != "gb" and sys.argv[1].split(".")[-1] != "GB" and sys.argv[1].split(".")[-1] != "fasta" and sys.argv[1].split(".")[-1] != "FASTA" and sys.argv[1].split(".")[-1] != "fas" and sys.argv[1].split(".")[-1] != "FAS" and sys.argv[1].split(".")[-1] != "fa" and sys.argv[1].split(".")[-1] != "FA"):
      print >> sys.stderr, "No EMBL/GBK/FASTA file submitted as input. Please supply a valid file with .embl / .gbk / .fasta extension. "
      print "Usage: antismash <query fasta/embl/gbk file>  [options]"
      print "From the command line, input antismash --help for more information."
      logfile.write("Input format error. Please supply a valid file with .embl / .gbk / .fasta extension.\n")
      logfile.write("Usage: antismash <query fasta/embl/gbk file>  [options]\n")
      logfile.write("From the command line, input antismash --help for more information.\n")
      logfile.close()
      sys.exit(1)
    #Define input filename and make fixes if necessary
    infile = sys.argv[1]
    try:
      testfile = open(infile,"r").read()
    except(IOError):
      print >> sys.stderr, "Please supply valid name for input file."
      print "Usage: antismash <query fasta/embl/gbk file>  [options]"
      print "From the command line, input antismash --help for more information."
      logfile = open("antismash.log","w")
      logfile.write("Input format error. Please supply valid name for infile.\n")
      logfile.write("Usage: antismash <query fasta/embl/gbk file>  [options]\n")
      logfile.write("From the command line, input antismash --help for more information.\n")
      logfile.close()
      sys.exit(1)
    #Parse absolute paths if found
    absolutepath = "n"
    if "/" in infile or "\\" in infile:
      absolutepath = "y"
      lastpos1 = infile.rfind("\\")
      lastpos2 = infile.rfind("/")
      lastpos = max([lastpos1,lastpos2])
      originpath = infile[:(lastpos + 1)]
      infile = infile[(lastpos + 1):]
      if sys.platform == ('win32'):
        copycommand = 'copy/y "' + originpath + infile + '" ' + infile + ' > nul'
        os.system(copycommand)
      if sys.platform == ('linux2'):
        copycommand = 'cp ' + originpath + infile + " . > /dev/null"
        os.system(copycommand)
    #genomename = ".".join(infile.split(".")[:-1])
    #for i in genomename:
    #  if i in '!"#$%&()*+,./:;=>?@[]^`{|}' or i in "'":
    #    genomename = genomename.replace(i,"")
    #    if "/" in genomename:
    #      genomename = genomename.rpartition("/")[2]
    #    if "\\" in genomename:
    #      genomename = genomename.rpartition("\\")[2]
    genomename = os.path.splitext(os.path.basename(infile))[0]
    if sys.platform == ('linux2'):
      if genomename !=  infile.split(".")[-2]:
        oldinfile = infile.replace("(","\\(").replace(")","\\)").replace("*","\\*").replace("&","\\&").replace("!","\\!").replace("$","\\$").replace("{","\\{").replace("}","\\}").replace("|","\\|").replace("`","\\`").replace("'","\\'").replace('"','\\"').replace('?','\\?')
        infile = genomename + "." + infile.split(".")[-1]
        if "/" in genomename:
          genomename = genomename.rpartition("/")[2]
        if "\\" in genomename:
          genomename = genomename.rpartition("\\")[2]
        os.system("cp " + oldinfile + " " + infile)
    #Define outputfolder
    if absolutepath == "y":
      if sys.platform == ('win32'):
        dir1 = os.popen("dir/w/ad " + originpath)
        dir2 = os.popen("dir/w/ad")
        dir1 = dir1.read()
        dir2 = dir2.read()
      if sys.platform == ('linux2'):
        dir1 = os.popen("ls")
        dir2 = os.popen("ls " + originpath)
        dir1 = dir1.read()
        dir2 = dir2.read()
      parts = dir1.split(" ") + dir2.split(" ")
    else:
      if sys.platform == ('win32'):
        dir = os.popen("dir/w/ad")
        dir = dir.read()
      if sys.platform == ('linux2'):
        dir = os.popen("ls")
        dir = dir.read()
      parts = dir.split(" ")
    parts2 = []
    for i in parts:
      partparts = i.split("\n")
      for i in partparts:
        i = i.replace("[","")
        i = i.replace("]","")
        parts2.append(i)
    parts = parts2
    oldgenomename = genomename
    if genomename in parts:
      genomename = genomename + "_" + str(0)
      while genomename in parts:
        finalpart = genomename.split("_")[-1]
        allnumbers = "y"
        for i in finalpart:
          if i not in ["0","1","2","3","4","5","6","7","8","9"]:
            allnumbers = "n"
        if allnumbers == "y" and int(finalpart) in range(0,1000):
          newgenomename = ""
          for i in genomename.split("_")[:-1]:
            newgenomename = newgenomename + "_" + i
          newgenomename = newgenomename + "_" + str(int(finalpart) + 1)
        genomename = newgenomename[1:]
      genomename = genomename.replace("__","_")
    #Output results folder name for output checking by GUI
    resultslocfile = open("resultsfolder.txt","w")
    resultslocfile.write(os.getcwd() + os.sep + genomename)
    resultslocfile.close()
  #Implement defaults
  glimmertransl_table = str(1)
  genomeconf = "l"
  minglength = str(90)
  cores = "all"
  taxon = "p"
  clusterblast = "y"
  smcogs = "y"
  fullblast = "n"
  fullhmm = "n"
  if sys.platform == ('win32'):
    blastdbpath = '"' + os.getcwd() + "/clusean/db" + '"'
  if sys.platform == ('linux2'):
    blastdbpath = os.getcwd() + "/clusean/db"
  if sys.platform == ('win32'):
    pfamdbpath = '"' + os.getcwd() + "/clusean/db/" + '"'
  if sys.platform == ('linux2'):
    pfamdbpath = os.getcwd() + "/clusean/db/"
  geneclustertypes = [1]
  #Read user-specified options which may override defaults
  if len(sys.argv) > 2 or sys.argv[1] == "--help":
    options = sys.argv
    if "--" in options[-1] and sys.argv[1] != "--help":
      invalidoptions(options[-1])
    #identify option identifiers
    identifiers = []
    for i in options:
      if "--" in i:
        if i not in identifiers:
          identifiers.append(i)
        else:
          invalidoptions("No '--' in given options or option given twice.")
    for i in identifiers:
      if i != "--help":
        value = options[options.index(i) + 1].strip()
      if i == "--gtransl":
        for k in value:
          if k not in ["0","1","2","3","4","5","6","7","8","9"]:
            invalidoptions(i + "input is no number")
        if int(value) in range(1,24) and int(value) != 7 and int(value) != 8 and int(value) != 17 and int(value) != 18 and int(value) != 19 and int(value) != 20:
          glimmertransl_table = value
        else:
          invalidoptions(i)
      elif i == "--genomeconf":
        if value == "l" or value == "c":
          genomeconf = value
        else:
          invalidoptions(i)
      elif i == "--minglength":
        for k in value:
          if k not in ["0","1","2","3","4","5","6","7","8","9"]:
            invalidoptions(i)
        if int(value) in range(30,91):
          minglength = value
        else:
          print >> sys.stderr, "Invalid options input: minimal gene length should be a number between 30-90."
          logfile = open("antismash.log","w")
          logfile.write("Invalid options input: minimal gene length should be a number between 30-90.\n")
          logfile.close()
          sys.exit(1)
      elif i == "--cores":
        for k in value:
          if k not in ["0","1","2","3","4","5","6","7","8","9"]:
            invalidoptions(i)
        if int(value) in range(1,1000):
          cores = int(value)
        else:
          invalidoptions(i)
      elif i == "--taxon":
        if value == "p" or value == "e":
          taxon = value
        else:
          invalidoptions(i)
      elif i == "--clusterblast":
        if value == "y" or value == "n":
          clusterblast = value
        else:
          invalidoptions(i)
      elif i == "--smcogs":
        if value == "y" or value == "n":
          smcogs = value
        else:
          invalidoptions(i)
      elif i == "--fullblast":
        if value == "y" or value == "n":
          fullblast = value
        else:
          invalidoptions(i)
      elif i == "--fullhmm":
        if value == "y" or value == "n":
          fullhmm = value
        else:
          invalidoptions(i)
      elif i == "--glimmer_prediction":
        glimmer_prediction_path = value
      elif i == "--blastdbpath":
        if sys.platform == ('win32'):
          if options[options.index(i) + 1][0] != '"':
            value = '"' + options[options.index(i) + 1] + '"'
          else:
            value = options[options.index(i) + 1]
          if ":\\" in value:
            blastdbpath = value
          elif "\\" in value or "/" in value:
            if value[0] == "\\" or value[0] == "/":
              blastdbpath = os.getcwd() + value
            else:
              blastdbpath = os.getcwd() + "\\" + value
          else:
            blastdbpath = os.getcwd() + "\\" + value
        if sys.platform == ('linux2'):
          value = options[options.index(i) + 1]
          if "\\" in value or "/" in value:
            value = value.replace("\\","/")
            if value[0] == "/":
              blastdbpath = value
            else:
              blastdbpath = os.getcwd() + "/" + value
          else:
            blastdbpath = os.getcwd() + "/" + value
      elif i == "--pfamdbpath":
        if sys.platform == ('win32'):
          if options[options.index(i) + 1][0] != '"':
            value = '"' + options[options.index(i) + 1] + '"'
          else:
            value = options[options.index(i) + 1]
          if ":\\" in value:
            pfamdbpath = value
          elif "\\" in value or "/" in value:
            if value[0] == "\\" or value[0] == "/":
              pfamdbpath = os.getcwd() + value
            else:
              pfamdbpath = os.getcwd() + "\\" + value
          else:
            pfamdbpath = os.getcwd() + "\\" + value
        if sys.platform == ('linux2'):
          value = options[options.index(i) + 1]
          if "\\" in value or "/" in value:
            value = value.replace("\\","/")
            if value[0] == "/":
              pfamdbpath = value
            else:
              pfamdbpath = os.getcwd() + "/" + value
          else:
            pfamdbpath = os.getcwd() + "/" + value
      elif i == "--geneclustertypes":
        if "," not in value and value not in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"]:
          invalidoptions(i)
        else:
          types = value.split(",")
          types2 = []
          if "1" in types:
            types2 = [1]
          for j in types:
            if int(j) not in range(1,20):
              invalidoptions(i)
            else:
              types2.append(int(j))
          geneclustertypes = types2
      elif i == "--help":
        print inputinstructions
        sys.exit()
      elif i == "--gui":
        pass
      else:
        invalidoptions(i)

  #Determine number of CPUs used
  if cores == "all":
    try:
      nrcpus = multiprocessing.cpu_count()
    except(IOError,OSError,NotImplementedError):
      nrcpus = 1
  else:
    try:
      nrcpus = multiprocessing.cpu_count()
    except(IOError,OSError,NotImplementedError):
      nrcpus = 1
    if cores < nrcpus:
      nrcpus = cores

  #Create directory structure needed for file storage
  try:
    os.mkdir(genomename)
  except(IOError,OSError):
    pass
  hmmoutputfolder = genomename + "/hmmoutput/"
  try:
    os.mkdir(hmmoutputfolder)
  except(IOError,OSError):
    pass
  nrpspksoutputfolder = genomename + "/nrpspks/"
  try:
    os.mkdir(nrpspksoutputfolder)
  except(IOError,OSError):
    pass
  nrpspredictoroutputfolder = nrpspksoutputfolder + "nrpspredictor/"
  try:
    os.mkdir(nrpspredictoroutputfolder)
  except(IOError,OSError):
    pass
  minowanrpsoutputfolder = nrpspksoutputfolder + "minowanrpspred/"
  try:
    os.mkdir(minowanrpsoutputfolder)
  except(IOError,OSError):
    pass
  minowapksoutputfolder = nrpspksoutputfolder + "minowapkspred/"
  try:
    os.mkdir(minowapksoutputfolder)
  except(IOError,OSError):
    pass
  minowacaloutputfolder = nrpspksoutputfolder + "minowacalpred/"
  try:
    os.mkdir(minowacaloutputfolder)
  except(IOError,OSError):
    pass
  pkssignatureoutputfolder = nrpspksoutputfolder + "pkssignatures/"
  try:
    os.mkdir(pkssignatureoutputfolder)
  except(IOError,OSError):
    pass
  kranalysisoutputfolder = nrpspksoutputfolder + "kr_analysis/"
  try:
    os.mkdir(kranalysisoutputfolder)
  except(IOError,OSError):
    pass
  clusterblastoutputfolder = genomename + "/clusterblast/"
  try:
    os.mkdir(clusterblastoutputfolder)
  except(IOError,OSError):
    pass
  smcogsoutputfolder = genomename + "/smcogs/"
  try:
    os.mkdir(smcogsoutputfolder)
  except(IOError,OSError):
    pass
  substrspecsfolder = genomename + "/substrspecs/"
  try:
    os.mkdir(substrspecsfolder)
  except(IOError,OSError):
    pass
  structuresfolder = genomename + "/structures/"
  try:
    os.mkdir(structuresfolder)
  except(IOError,OSError):
    pass
  svgfolder = genomename + "/svg/"
  try:
    os.mkdir(svgfolder)
  except(IOError,OSError):
    pass
  searchgtrfolder = genomename + "/searchgtr/"
  try:
    os.mkdir(searchgtrfolder)
  except(IOError,OSError):
    pass
  htmlfolder = genomename + "/html/"
  try:
    os.mkdir(htmlfolder)
  except(IOError,OSError):
    pass
  imagesfolder = genomename + "/images/"
  try:
    os.mkdir(imagesfolder)
  except(IOError,OSError):
    pass

  #If input is unannotated GBK/EMBL file, convert to FASTA and use that as input
  if "     CDS             " not in open(infile,"r").read() and "FT   CDS " not in open(infile,"r").read():
    if infile.split(".")[-1] == "embl" or infile.split(".")[-1] == "EMBL" or infile.split(".")[-1] == "emb" or infile.split(".")[-1] == "EMB":
      filetext = open(infile,"r").read()
      if "\nSQ" not in filetext:
        print >> sys.stderr, "Exit: EMBL file not properly formatted, no sequence found."
        logfile = open("antismash.log","w")
        logfile.write("Exit: EMBL file not properly formatted, no sequence found.\n")
        logfile.close()
        sys.exit(1)
      dnaseq = filetext.split("\nSQ")[1]
      dnaseq = cleandnaseq(dnaseq)
      sequence = dnaseq
      if (sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
        print >> sys.stderr, "Protein EMBL file provided. Please provide nucleotide EMBL file."
        sys.exit(1)
      fastafile = open(infile.rpartition(".")[0] + ".fasta","w")
      fastafile.write(">" + infile.rpartition(".")[0] + "|\n")
      fastafile.write(sequence)
      fastafile.close()
      infile = fastafile
    elif infile.split(".")[-1] == "gbk" or infile.split(".")[-1] == "GBK" or infile.split(".")[-1] == "gb" or infile.split(".")[-1] == "GB" or infile.split(".")[-1] == "genbank" or infile.split(".")[-1] == "GENBANK":
      filetext = open(infile,"r").read()
      if "\nORIGIN" not in filetext:
        print >> sys.stderr, "Exit: GBK file not properly formatted, no sequence found."
        logfile = open("antismash.log","w")
        logfile.write("Exit: GBK file not properly formatted, no sequence found.\n")
        logfile.close()
        sys.exit(1)
      dnaseq = filetext.split("\nORIGIN")[1]
      dnaseq = cleandnaseq(dnaseq)
      sequence = dnaseq
      if (sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
        print >> sys.stderr, "Protein GBK file provided. Please provide nucleotide GBK file."
        sys.exit(1)
      fastafile = open(infile.rpartition(".")[0] + ".fasta","w")
      fastafile.write(">" + infile.rpartition(".")[0] + "\n")
      fastafile.write(sequence)
      fastafile.close()
      infile = infile.rpartition(".")[0] + ".fasta"
  #If input is unannotated fasta file, predict genes with Glimmer and create EMBL file. If input is EMBL or GBK file, read input embl/gbk and create input fasta file, read input protein info into memory
  annotated = "y"
  if infile.split(".")[-1] == "fasta" or infile.split(".")[-1] == "FASTA" or infile.split(".")[-1] == "FAS" or infile.split(".")[-1] == "fas" or infile.split(".")[-1] == "FA" or infile.split(".")[-1] == "fa":
    annotated = "n"
    #Check input file formatting
    sequence = get_sequence(infile)
    if (sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
      print >> sys.stderr, "Protein FASTA file provided. Please provide nucleotide FASTA file."
      sys.exit(1)
    nucleotides = ["A","a","C","c","G","g","T","t","N","n"]
    badsequence = "n"
    sequence_name = open(infile,"r").read().partition(">")[2].partition("\n")[0]
    for i in sequence:
      if i not in nucleotides:
        badsequence = "y"
    if badsequence == "y":
      cleaned_sequence = cleandnaseq(sequence)
      badsequence = "n"
      for i in cleaned_sequence:
        if i not in nucleotides:
          badsequence = "y"
      if badsequence == "n":
        writefasta([sequence_name],[cleaned_sequence],infile.rpartition(".")[0] + "_f.fasta")
        infile = infile.rpartition(".")[0] + "_f.fasta"
      else:
        print >>sys.stderr, "Incorrect file formatting. Please submit a properly formatted single-sequence FASTA file."
        logfile = open("antismash.log","w")
        logfile.write("Incorrect file formatting. Please submit a properly formatted single-sequence FASTA file.\n")
        logfile.close()
        sys.exit(1)
    revseq = reverse_complement(sequence)
    seqlength = len(sequence)

    #Print Glimmer notification
    #if taxon == "p":
    #  print "Running Glimmer 3.02 to predict genes in unannotated prokaryotic genome..."
    #elif taxon == "e":
    #  print "Running GlimmerHMM 3.0.1 to predict genes in unannotated eukaryotic genome..."
    logfile = open("antismash.log","w")
    if taxon == "p":
      logfile.write("Running Glimmer 3.02 to predict genes in unannotated prokaryotic genome...\n")
    elif taxon == "e":
      logfile.write("Running GlimmerHMM 3.0.1 to predict genes in unannotated eukaryotic genome...\n")
    #logfile.close()
    loginfo = open("antismash.log","r").read()
    #logfile.close()
    #Copying file and changing to folder to prepare for Glimmer3 prediction
    os.mkdir( os.path.join(os.getcwd(), genomename, "geneprediction"))
    if sys.platform == ('win32'):
      os.system("copy/y " + infile + " geneprediction > nul")
    if sys.platform == ('linux2'):
      os.system("cp " + infile + " geneprediction > /dev/null")

    os.chdir( os.path.join(os.getcwd(), genomename, "geneprediction"))
    fastafile = '../../'+infile

    #Find DNA sequence length
    seq = get_sequence(fastafile)
    dnaseqlength = len(seq)
    #Run Glimmer for prokaryotic sequences, GlimmerHMM for eukaryotic sequences
    if taxon == "p":
      """
      GlimmerPrediction, not needed since we can predict it in galaxy on our own
      if genomeconf == "l":
        if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
          os.popen("tigr-glimmer long-orfs -l -n -t 1.15 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".longorfs")
        else:
          os.system("tigr-glimmer long-orfs -l -n -t 1.15 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".longorfs")
      else:
        if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
          os.popen("tigr-glimmer long-orfs -n -t 1.15 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".longorfs")
        else:
          os.system("tigr-glimmer long-orfs -n -t 1.15 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".longorfs")
      if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
          os.popen("tigr-glimmer extract -t " + fastafile + " " + fastafile.rpartition(".")[0] + ".longorfs > " + fastafile.rpartition(".")[0] + ".train")
      else:
        os.system("tigr-glimmer extract -t " + fastafile + " " + fastafile.rpartition(".")[0] + ".longorfs > " + fastafile.rpartition(".")[0] + ".train")
      if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
        os.popen("tigr-glimmer build-icm -r " + fastafile.rpartition(".")[0] + ".icm < " + fastafile.rpartition(".")[0] + ".train")
      else:
        os.system("tigr-glimmer build-icm -r " + fastafile.rpartition(".")[0] + ".icm < " + fastafile.rpartition(".")[0] + ".train")
      if genomeconf == "l":
        if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
          os.popen("tigr-glimmer glimmer3 -l -o50 -g" + minglength + " -q3000 -t30 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".icm " + fastafile.rpartition(".")[0])
        else:
          os.system("tigr-glimmer glimmer3 -l -o50 -g" + minglength + " -q3000 -t30 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".icm " + fastafile.rpartition(".")[0])
      else:
        if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
          os.popen("tigr-glimmer glimmer3 -o50 -g" + minglength + " -q3000 -t30 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".icm " + fastafile.rpartition(".")[0])        
        else:
          os.system("tigr-glimmer glimmer3 -o50 -g" + minglength + " -q3000 -t30 --trans_table " + glimmertransl_table + " " + fastafile + " " + fastafile.rpartition(".")[0] + ".icm " + fastafile.rpartition(".")[0])
      #Convert glimmer predictions into EMBL with sequence
      glfile = fastafile.rpartition(".")[0] + ".predict"
      
      Ende der Glimmer-Prediction
      """
      glfile = glimmer_prediction_path
      emblfile = fastafile.rpartition(".")[0] + ".embl"
      try:
        file = open(glfile,"r")
        filetext = file.read()
      except:
        print >> sys.stderr, "Glimmer gene prediction failed. Please check the format of your input FASTA file. Error 11."
        logfile = open("antismash.log","w")
        logfile.write("Glimmer gene prediction failed. Please check the format of your input FASTA file. Error 11.\n")
        logfile.close()
        sys.exit(1)
      if "orf" not in filetext:
        print >> sys.stderr, "Glimmer gene prediction failed: no genes found."
        logfile = open("antismash.log","w")
        logfile.write("Glimmer gene prediction failed: no genes found.\n")
        logfile.close()
        sys.exit(1)
      filetext = filetext.replace("\r","\n")
      lines = filetext.split("\n")
      lines = lines[1:-1]
      orfnames = []
      starts = []
      ends = []
      strands = []
      starts2 = []
      ends2 = []
      firstline = "y"
      for i in lines:
        columns = i.split(" ")
        columns2 = []
        for i in columns:
          if i != "":
            columns2.append(i)
        columns = columns2
        if len(columns) > 3:
          frame = columns[3][0]
          strands.append(frame)
        else:
          frame = ""
        if firstline == "y" and frame == "+" and len(columns) > 3:
          orfname = str(columns[0])
          orfnames.append(orfname)
          if genomeconf == "c" and (int(columns[1]) > int(columns[2])) and (int(columns[1]) > (0.5 * dnaseqlength)):
            gstart = (int(columns[2]) % 3) + 1
            if gstart == 3:
              gstart = 0
            starts.append(str(gstart))
            ends.append(columns[2])
            starts.append(columns[1])
            ends.append(str(dnaseqlength))
          else:
            starts.append(columns[1])
            ends.append(columns[2])
          firstline = "n"
        elif firstline == "y" and frame == "-" and len(columns) > 3:
          orfname = str(columns[0])
          orfnames.append(orfname)
          if genomeconf == "c" and (int(columns[1]) > int(columns[2])) and (int(columns[1]) > (0.5 * dnaseqlength)):
            gstart = (int(columns[2]) % 3) + 1
            if gstart == 3:
              gstart = 0
            starts.append("complement(" + str(gstart))
            ends.append(columns[2] + ")")
            starts.append("complement(" + columns[1])
            ends.append(str(dnaseqlength) + ")")
          else:
            complstart = "complement(" + str(columns[1])   
            starts.append(complstart)
            complend = str(columns[2]) + ")"
            ends.append(str(complend))
          firstline = "n"
        elif frame == "+" and len(columns) > 3:
          orfname = str(columns[0])
          orfnames.append(orfname)
          starts.append(columns[1])
          ends.append(columns[2])
        elif frame == "-" and len(columns) > 3:
          orfname = str(columns[0])
          orfnames.append(orfname)
          complstart = "complement(" + str(columns[1])   
          starts.append(complstart)
          complend = str(columns[2]) + ")"
          ends.append(str(complend))
      if len(orfnames) == 0:
        print >> sys.stderr, "Glimmer gene prediction failed. Please check the format of your input FASTA file. Error 10."
        logfile = open("antismash.log","w")
        logfile.write("Glimmer gene prediction failed. Please check the format of your input FASTA file. Error 10.\n")
        logfile.close()
        sys.exit(1)
      out_file = open(emblfile,"w")
      a = 0
      #print "Writing EMBL file with Glimmer-predicted genes..."
      logfile = open("antismash.log","w")
      logfile.write(loginfo)
      logfile.write("Writing EMBL file with Glimmer-predicted genes...\n")
      #logfile.close()
      loginfo = open("antismash.log","r").read()
      #logfile.close()
      if taxon == "p":
        out_file.write("ID   A01; SV 1; linear; DNA; STD; PRO; " + str(dnaseqlength) + " BP.\nXX\n")
      else:
        out_file.write("ID   A01; SV 1; linear; DNA; STD; FUN; " + str(dnaseqlength) + " BP.\nXX\n")
      out_file.write("AC   A01;\nXX\n")
      out_file.write("DE   " + genomename + ";\nXX\n")
      out_file.write("KW   none;\nXX\n")
      out_file.write("OS   unknown;\n")
      if taxon == "p":
        out_file.write("OC   Eubacteria;\nXX\n")
      else:
        out_file.write("OC   Fungi;\nXX\n")
      out_file.write("RN   [1]\n")
      out_file.write("RT   ;\n")
      out_file.write("RL   Unknown.\nXX\n")
      out_file.write("FH   Key             Location/Qualifiers\nFH\n")
      out_file.write("FT   source          1.." + str(dnaseqlength) + "\n")
      for i in orfnames:
        out_file.write("FT   gene            ")
        out_file.write(starts[a])
        out_file.write("..")
        out_file.write(ends[a])
        out_file.write("\n")
        out_file.write('FT                   /gene="' + i + '"\n')
        out_file.write("FT   CDS             ")
        out_file.write(starts[a])
        out_file.write("..")
        out_file.write(ends[a])
        out_file.write("\n")
        out_file.write('FT                   /gene="' + i + '"\n')
        a += 1
    elif taxon == "e":
      """
      GlimmerHMM is executed extern ... in galaxy and will be provided through glimmer_prediction_path
      
      if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
        os.popen("glimmerhmm " + fastafile + " train_crypto -o " + fastafile.rpartition(".")[0] + ".predict -g")
      else:
        os.system("glimmerhmm " + fastafile + " train_crypto -o " + fastafile.rpartition(".")[0] + ".predict -g")
      """
      #Convert glimmerhmm predictions into EMBL with sequence
      #glfile = fastafile.rpartition(".")[0] + ".predict"
      glfile = glimmer_prediction_path
      emblfile = fastafile.rpartition(".")[0] + ".embl"
      try:
        file = open(glfile,"r")
        filetext = file.read().replace("\r","")
      except:
        print >> sys.stderr, "GlimmerHMM gene prediction failed. Please check the format of your input FASTA file. Error 9."
        logfile = open("antismash.log","w")
        logfile.write("GlimmerHMM gene prediction failed. Please check the format of your input FASTA file. Error 9.\n")
        logfile.close()
        sys.exit(1)
      if "CDS" not in filetext:
        print >> sys.stderr, "GlimmerHMM gene prediction failed: no genes found."
        logfile = open("antismash.log","w")
        logfile.write("GlimmerHMM gene prediction failed: no genes found.\n")
        logfile.close()
        sys.exit(1)
      filetext = filetext.replace("\r","\n")
      lines = filetext.split("\n")
      lines = lines[2:-1]
      orfnames = []
      positions = []
      firstline = "y"
      x = 0
      orfnr = 0
      starts = []
      ends = []
      for i in lines:
        columns = i.split("\t")
        if len(columns) > 1:
          if x == 0:
            strand = columns[6]
            if "mRNA" not in i:
              starts.append(columns[3])
              ends.append(columns[4])
          elif x == (len(lines) - 1) or "mRNA" in lines[x + 1]:
            strand = columns[6]
            starts.append(columns[3])
            ends.append(columns[4])
            orfnr += 1
            if len(str(orfnr)) == 1:
              orfname = "orf0000" + str(orfnr)
            elif len(str(orfnr)) == 2:
              orfname = "orf000" + str(orfnr)
            elif len(str(orfnr)) == 3:
              orfname = "orf00" + str(orfnr)
            elif len(str(orfnr)) == 4:
              orfname = "orf0" + str(orfnr)
            elif len(str(orfnr)) == 5:
              orfname = "orf" + str(orfnr)
            orfnames.append(orfname)
            if strand == "+":
              if len(starts) == 1:
                pos = starts[0] + ".." + ends[0]
                positions.append(pos)
              else:
                pos = "join("
                y = 0
                for i in starts:
                  pos = pos + i + ".." + ends[y]
                  if i != starts[-1]:
                    pos = pos + ","
                  y += 1
                pos = pos + ")"
                positions.append(pos)
            elif strand == "-":
              if len(starts) == 1:
                pos = "complement(" + starts[0] + ".." + ends[0] + ")"
                positions.append(pos)
              else:
                pos = "complement(join("
                y = 0
                for i in starts:
                  pos = pos + i + ".." + ends[y]
                  if i != starts[-1]:
                    pos = pos + ","
                  y += 1
                pos = pos + "))"
                positions.append(pos)
            starts = []
            ends = []
          elif "mRNA" not in i:
            starts.append(columns[3])
            ends.append(columns[4])
        x += 1
      if len(orfnames) == 0:
        print >> sys.stderr, "GlimmerHMM gene prediction failed. Please check the format of your input FASTA file. Error: 12"
        logfile = open("antismash.log","w")
        logfile.write("GlimmerHMM gene prediction failed. Please check the format of your input FASTA file. Error 12\n")
        logfile.close()
        sys.exit(1)
      out_file = open(emblfile,"w")
      a = 0
      #print "Writing EMBL file with GlimmerHMM-predicted genes..."
      logfile = open("antismash.log","w")
      logfile.write(loginfo)
      logfile.write("Writing EMBL file with GlimmerHMM-predicted genes...\n")
      #logfile.close()
      loginfo = open("antismash.log","r").read()
      #logfile.close()
      out_file.write("ID   A01; SV 1; linear; DNA; STD; FUN; " + str(dnaseqlength) + " BP.\nXX\n")
      out_file.write("AC   A01;\nXX\n")
      out_file.write("DE   " + genomename + ";\nXX\n")
      out_file.write("KW   none;\nXX\n")
      out_file.write("OS   unknown;\n")
      out_file.write("OC   Fungi;\nXX\n")
      out_file.write("RN   [1]\n")
      out_file.write("RT   ;\n")
      out_file.write("RL   Unknown.\nXX\n")
      out_file.write("FH   Key             Location/Qualifiers\nFH\n")
      out_file.write("FT   source          1.." + str(dnaseqlength) + "\n")
      for i in orfnames:
        out_file.write("FT   gene            ")
        out_file.write(positions[a])
        out_file.write("\n")
        out_file.write('FT                   /gene="' + i + '"\n')
        out_file.write("FT   CDS             ")
        out_file.write(positions[a])
        out_file.write("\n")
        out_file.write('FT                   /gene="' + i + '"\n')
        a += 1
    out_file.write("XX\nSQ   Sequence " + str(dnaseqlength) + " BP; " + str(seq.count("a") + seq.count("A")) + " A; " + str(seq.count("c") + seq.count("C")) + " C; " + str(seq.count("g") + seq.count("G")) + " G; " + str(seq.count("t") + seq.count("T")) + " T; " + str(dnaseqlength - (seq.count("a") + seq.count("A") + seq.count("c") + seq.count("C") + seq.count("g") + seq.count("G") + seq.count("t") + seq.count("T"))) + " other;\n")
    seq2 = seq
    out_file.write("     ")
    grouplen=10
    textlen = len(seq)
    end = textlen - (textlen % grouplen)
    repeated_iterator = [iter(itertools.islice(seq, 0, end))] * grouplen
    parts = list(itertools.imap(lambda *chars: ''.join(chars),*repeated_iterator))
    if dnaseqlength%grouplen != 0:
      parts.append(seq[-1 * (dnaseqlength%grouplen):])
    w = 1
    for l in parts:
      out_file.write(l + " ")
      if w == len(parts):
        if w%6 == 0 and dnaseqlength%60 != 0:
          out_file.write((" " * (10 - dnaseqlength%grouplen) + " " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
        elif dnaseqlength%60 == 0:
          out_file.write((" " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
        elif w%6 == 5 and dnaseqlength%grouplen == 0:
          out_file.write(("           " + " " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
        elif dnaseqlength%grouplen != 0:
          out_file.write(" " * (10 - dnaseqlength%grouplen) + "          " * (6 - len(parts)%6) + " " * (6 - len(parts)%6) + (" " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
        else:
          out_file.write("          " * (6 - len(parts)%6) + " " * (5 - len(parts)%6) + (" " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
      elif w%6 == 0:
        out_file.write((" " * (10 - len(str(w * 10)))) + str(w * 10) + "\n     ")
      w += 1
    out_file.close()
    os.chdir("../../")
    infile = emblfile[6:]
    emblfile = emblfile[6:]
    if taxon == "p":
      glimmeroutputfolder = genomename + "/glimmer/"
    elif taxon == "e":
      glimmeroutputfolder = genomename + "/glimmerhmm/"
    try:
      os.mkdir(glimmeroutputfolder)
    except(IOError,OSError):
      pass
    proteins = embl2proteins(infile,sequence)
    genomic_accnr = proteins[1]
    dnaseqlength = proteins[2]
    proteins = proteins[0]
    writefasta(proteins[0],proteins[1],genomename + "/genome_proteins.fasta")
  else:
    #print "Reading embl/gbk file and creating input FASTA file for gene cluster detection..."
    logfile.write("Reading embl/gbk file and creating input FASTA file for gene cluster detection...\n")
    if infile.split(".")[-1] == "embl" or infile.split(".")[-1] == "EMBL" or infile.split(".")[-1] == "emb" or infile.split(".")[-1] == "EMB":
      sequence = ""
      proteins = embl2proteins(infile,sequence)
      genomic_accnr = proteins[1]
      dnaseqlength = proteins[2]
      proteins = proteins[0]
      writefasta(proteins[0],proteins[1],genomename + "/genome_proteins.fasta")
    elif infile.split(".")[-1] == "gbk" or infile.split(".")[-1] == "GBK" or infile.split(".")[-1] == "gb" or infile.split(".")[-1] == "GB" or infile.split(".")[-1] == "genbank" or infile.split(".")[-1] == "GENBANK":
      proteins = gbk2proteins(infile)
      genomic_accnr = proteins[1]
      dnaseqlength = proteins[2]
      proteins = proteins[0]
      writefasta(proteins[0],proteins[1],genomename + "/genome_proteins.fasta")
  accessiondict = proteins[4]
  seqdict = {}
  fullnamedict = {}
  strandsdict = {}
  z = 0
  for i in proteins[0]:
    name = i.split("|")[4]
    seq = proteins[1][z]
    seqdict[name] = seq
    strand = i.split("|")[3]
    strandsdict[name] = strand
    fullnamedict[name] = i
    z += 1

  elapsed = (time.time() - starttime)
  #print "2968Time since start: " + str(elapsed)

  #Run hmmsearch on proteins from input file and parse output
  #print "Performing HMM search on proteins for detection of signature genes..."
  logfile.write("Performing HMM search on proteins for detection of signature genes...\n")
  hmmslist = ["AMP-binding.hmm","BLS.hmm","CAS.hmm","Chal_sti_synt_C.hmm","Chal_sti_synt_N.hmm","Condensation.hmm","ene_KS.hmm","hyb_KS.hmm","itr_KS.hmm","mod_KS.hmm","tra_KS.hmm","LANC_like.hmm","ATd.hmm","PKS_AT.hmm","PKS_KS.hmm","PP-binding.hmm","t2clf.hmm","t2ks.hmm","t2ks2.hmm","Terpene_synth.hmm","Terpene_synth_C.hmm","strH_like.hmm","neoL_like.hmm","DOIS.hmm","valA_like.hmm","spcFG_like.hmm","spcDK_like_cou.hmm","spcDK_like_glyc.hmm","strK_like1.hmm","strK_like2.hmm","bt1fas.hmm","ft1fas.hmm","t2fas.hmm","hglD.hmm","hglE.hmm","fabH.hmm","AfsA.hmm","IucA_IucC.hmm","ectoine_synt.hmm","phytoene_synt.hmm","Lant_dehyd_N.hmm","Lant_dehyd_C.hmm","Antimicrobial18.hmm","Gallidermin.hmm","L_biotic_typeA.hmm","LE-DUF.hmm","LE-LAC481.hmm","LE-LanBC.hmm","LE-MER+2PEP.hmm","MA-2PEPA.hmm","MA-DUF.hmm","MA-EPI.hmm","MA-LAC481.hmm","MA-NIS+EPI.hmm","MA-NIS.hmm","indsynth.hmm","A-OX.hmm","LmbU.hmm","MoeO5.hmm","LipM.hmm","LipU.hmm","LipV.hmm","ToyB.hmm","TunD.hmm","melC.hmm","strepbact.hmm","goadsporin_like.hmm","Antimicrobial14.hmm","Bacteriocin_IId.hmm","BacteriocIIc_cy.hmm","Bacteriocin_II.hmm","Lactococcin.hmm","Antimicrobial17.hmm","Lactococcin_972.hmm","Bacteriocin_IIc.hmm","LcnG-beta.hmm","Bacteriocin_IIi.hmm","Subtilosin_A.hmm","Cloacin.hmm","Neocarzinostat.hmm","Linocin_M18.hmm","TIGR03603.hmm","TIGR03604.hmm","TIGR03605.hmm","TIGR03731.hmm","TIGR03651.hmm","TIGR03678.hmm","TIGR03693.hmm","TIGR03798.hmm","TIGR03882.hmm","TIGR03601.hmm","TIGR03602.hmm","tabtoxin.hmm","cycdipepsynth.hmm","cyanobactin_synth.hmm","fom1.hmm","bcpB.hmm","frbD.hmm","mitE.hmm",'Lycopene_cycl.hmm','terpene_cyclase.hmm','NapT7.hmm','fung_ggpps.hmm','fung_ggpps2.hmm','dmat.hmm','trichodiene_synth.hmm','novK.hmm','novJ.hmm','novI.hmm','novH.hmm','pur6.hmm','pur10.hmm','nikJ.hmm','nikO.hmm','mvnA.hmm','thiostrepton.hmm','NAD_binding_4.hmm','vlmB.hmm','salQ.hmm','prnB.hmm']
  for i in hmmslist:
    hmmsearch = hmmsearch_path + " " + "--cpu " + str(nrcpus) + " -o " + genomename + "/hmmoutput/" + i.split(".")[0] + "_output.txt" + " --noali --tblout " + genomename + "/hmmoutput/" + i.split(".")[0] + ".txt " + hmms_path + i + " " + genomename + "/genome_proteins.fasta"
    os.system(hmmsearch)
  #print "Parsing HMM outputs..."
  logfile.write("Parsing HMM outputs...\n")
  detecteddomainsdict = {}
  #Extract type I PKS proteins, KS cut-off: 50; AT cut-off: 20; exclude those sequences that score higher on type I FAS HMMs, type IV hglE-like KS domains
  t1pksprots = []
  transatpksprots = []
  if 1 in geneclustertypes or 2 in geneclustertypes or 3 in geneclustertypes or 4 in geneclustertypes:
    ks = parsehmmoutput(50,hmmoutputfolder + "PKS_KS.txt")
    at = parsehmmoutput(50,hmmoutputfolder + "PKS_AT.txt")
    ft1fasks = parsehmmoutput(50,hmmoutputfolder + "ft1fas.txt")
    bt1fasks = parsehmmoutput(50,hmmoutputfolder + "bt1fas.txt")
    hgleks = parsehmmoutput(50,hmmoutputfolder + "hglE.txt")
    hgldks = parsehmmoutput(50,hmmoutputfolder + "hglD.txt")
    fabhks = parsehmmoutput(50,hmmoutputfolder + "fabH.txt")
    pksksprots = ks[0]
    pksatprots = at[0]
    pksatscores = at[1]
    pksksscores = ks[1]
    bt1fasprots = bt1fasks[0]
    bt1fasscores = bt1fasks[1]
    ft1fasprots = ft1fasks[0]
    ft1fasscores = ft1fasks[1]
    hgleprots = hgleks[0]
    hglescores = hgleks[1]
    hgldprots = hgldks[0]
    hgldscores = hgldks[1]
    fabhprots = fabhks[0]
    fabhscores = fabhks[1]
    for i in pksksprots:
      exclude = "n"
      score = pksksscores[pksksprots.index(i)]
      if i in bt1fasprots:
        bt1fasscore = bt1fasscores[bt1fasprots.index(i)]
        if float(score) < float(bt1fasscore):
          exclude = "y"
      if i in ft1fasprots:
        ft1fasscore = ft1fasscores[ft1fasprots.index(i)]
        if float(score) < float(ft1fasscore):
          exclude = "y"
      if i in hgldprots:
        hgldscore = hgldscores[hgldprots.index(i)]
        if float(score) < float(hgldscore):
          exclude = "y"
      if i in hgleprots:
        hglescore = hglescores[hgleprots.index(i)]
        if float(score) < float(hglescore):
          exclude = "y"
      if i in fabhprots:
        fabhscore = fabhscores[fabhprots.index(i)]
        if float(score) < float(fabhscore):
          exclude = "y"
      if i in pksatprots and exclude == "n":
        t1pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["PKS ketosynthase domain",pksksscores[pksksprots.index(i)]])
          detdomlist.append(["PKS acyltransferase domain",pksatscores[pksatprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["PKS ketosynthase domain",pksksscores[pksksprots.index(i)]],["PKS acyltransferase domain",pksatscores[pksatprots.index(i)]]]
    #Extract trans-AT PKSs: proteins with KS hits but without AT hits, and with trans-AT specific ATd-hits
    atd = parsehmmoutput(65,hmmoutputfolder + "ATd.txt")
    traks = parsehmmoutput(50,hmmoutputfolder + "tra_KS.txt")
    traksprots = traks[0]
    atdprots = atd[0]
    atdscores = atd[1]
    for i in pksksprots:
      if i in atdprots and i in traksprots and i not in t1pksprots:
        transatpksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["PKS ketosynthase domain",pksksscores[pksksprots.index(i)]])
          detdomlist.append(["Trans-AT PKS AT-docking domain",atdscores[atdprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["PKS ketosynthase domain",pksksscores[pksksprots.index(i)]],["Trans-AT PKS AT-docking domain",atdscores[atdprots.index(i)]]]
  #Extract type II PKS & CLF proteins, KS-cut-off: 50, t2KS/clf score > modKS,eneKS,itrKS,traKS,t1fas,t2fas,hgle scores
  t2pksprots = []
  if 1 in geneclustertypes or 2 in geneclustertypes or 3 in geneclustertypes or 4 in geneclustertypes:
    t2ks = parsehmmoutput(50,hmmoutputfolder + "t2ks.txt")
    t2ks2 = parsehmmoutput(450,hmmoutputfolder + "t2ks2.txt")
    t2clf = parsehmmoutput(50,hmmoutputfolder + "t2clf.txt")
    eneks = parsehmmoutput(50,hmmoutputfolder + "ene_KS.txt")
    hybks = parsehmmoutput(50,hmmoutputfolder + "hyb_KS.txt")
    modks = parsehmmoutput(50,hmmoutputfolder + "mod_KS.txt")
    itrks = parsehmmoutput(50,hmmoutputfolder + "itr_KS.txt")
    traks = parsehmmoutput(50,hmmoutputfolder + "tra_KS.txt")
    t2fasks = parsehmmoutput(50,hmmoutputfolder + "t2fas.txt")
    ft1fasks = parsehmmoutput(50,hmmoutputfolder + "ft1fas.txt")
    bt1fasks = parsehmmoutput(50,hmmoutputfolder + "bt1fas.txt")
    hgleks = parsehmmoutput(50,hmmoutputfolder + "hglE.txt")
    hgldks = parsehmmoutput(50,hmmoutputfolder + "hglD.txt")
    fabhks = parsehmmoutput(50,hmmoutputfolder + "fabH.txt")
    t2ksprots = t2ks[0]
    t2ks2prots = t2ks2[0]
    t2clfprots = t2clf[0]
    eneksprots = eneks[0]
    hybksprots = hybks[0]
    modksprots = modks[0]
    itrksprots = itrks[0]
    traksprots = traks[0]
    t2fasprots = t2fasks[0]
    t2ksscores = t2ks[1]
    t2ks2scores = t2ks2[1]
    t2clfscores = t2clf[1]
    eneksscores = eneks[1]
    hybksscores = hybks[1]
    modksscores = modks[1]
    itrksscores = itrks[1]
    traksscores = traks[1]
    t2fasscores = t2fasks[1]
    bt1fasprots = bt1fasks[0]
    bt1fasscores = bt1fasks[1]
    ft1fasprots = ft1fasks[0]
    ft1fasscores = ft1fasks[1]
    hgleprots = hgleks[0]
    hglescores = hgleks[1]
    hgldprots = hgldks[0]
    hgldscores = hgldks[1]
    fabhprots = fabhks[0]
    fabhscores = fabhks[1]
    for i in t2ksprots:
      type2 = "y"
      score = t2ksscores[t2ksprots.index(i)]
      if i in eneksprots:
        enescore = eneksscores[eneksprots.index(i)]
        if float(enescore) > float(score):
          type2 = "n"
      if i in hybksprots:
        hybscore = hybksscores[hybksprots.index(i)]
        if float(hybscore) > float(score):
          type2 = "n"
      if i in modksprots:
        modscore = modksscores[modksprots.index(i)]
        if float(modscore) > float(score):
          type2 = "n"
      if i in itrksprots:
        itrscore = itrksscores[itrksprots.index(i)]
        if float(itrscore) > float(score):
          type2 = "n"
      if i in traksprots:
        trascore = traksscores[traksprots.index(i)]
        if float(trascore) > float(score):
          type2 = "n"
      if i in bt1fasprots:
        bt1fasscore = bt1fasscores[bt1fasprots.index(i)]
        if float(bt1fasscore) > float(score):
          type2 = "n"
      if i in ft1fasprots:
        ft1fasscore = ft1fasscores[ft1fasprots.index(i)]
        if float(ft1fasscore) > float(score):
          type2 = "n"
      if i in t2fasprots:
        t2fasscore = t2fasscores[t2fasprots.index(i)]
        if float(t2fasscore) > float(score):
          type2 = "n"
      if i in hgleprots:
        hglescore = hglescores[hgleprots.index(i)]
        if float(hglescore) > float(score):
          type2 = "n"
      if i in fabhprots:
        fabhscore = fabhscores[fabhprots.index(i)]
        if float(fabhscore) > float(score):
          type2 = "n"
      if type2 == "y" and i not in t2pksprots and i not in t1pksprots:
        t2pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Type II ketosynthase",t2ksscores[t2ksprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Type II ketosynthase",t2ksscores[t2ksprots.index(i)]]]
    for i in t2clfprots:
      type2 = "y"
      score = t2clfscores[t2clfprots.index(i)]
      if i in eneksprots:
        enescore = eneksscores[eneksprots.index(i)]
        if float(enescore) > float(score):
          type2 = "n"
      if i in hybksprots:
        hybscore = hybksscores[hybksprots.index(i)]
        if float(hybscore) > float(score):
          type2 = "n"
      if i in modksprots:
        modscore = modksscores[modksprots.index(i)]
        if float(modscore) > float(score):
          type2 = "n"
      if i in itrksprots:
        itrscore = itrksscores[itrksprots.index(i)]
        if float(itrscore) > float(score):
          type2 = "n"
      if i in traksprots:
        trascore = traksscores[traksprots.index(i)]
        if float(trascore) > float(score):
          type2 = "n"
      if i in bt1fasprots:
        bt1fasscore = bt1fasscores[bt1fasprots.index(i)]
        if float(bt1fasscore) > float(score):
          type2 = "n"
      if i in ft1fasprots:
        ft1fasscore = ft1fasscores[ft1fasprots.index(i)]
        if float(ft1fasscore) > float(score):
          type2 = "n"
      if i in t2fasprots:
        t2fasscore = t2fasscores[t2fasprots.index(i)]
        if float(t2fasscore) > float(score):
          type2 = "n"
      if i in hgleprots:
        hglescore = hglescores[hgleprots.index(i)]
        if float(hglescore) > float(score):
          type2 = "n"
      if i in fabhprots:
        fabhscore = fabhscores[fabhprots.index(i)]
        if float(fabhscore) > float(score):
          type2 = "n"
      if type2 == "y" and i not in t2pksprots and i not in t1pksprots:
        t2pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Chain length factor",t2clfscores[t2clfprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Chain length factor",t2clfscores[t2clfprots.index(i)]]]
    for i in t2ks2prots:
      type2 = "y"
      score = t2ks2scores[t2ks2prots.index(i)]
      if i in eneksprots:
        enescore = eneksscores[eneksprots.index(i)]
        if float(enescore) > float(score):
          type2 = "n"
      if i in hybksprots:
        hybscore = hybksscores[hybksprots.index(i)]
        if float(hybscore) > float(score):
          type2 = "n"
      if i in modksprots:
        modscore = modksscores[modksprots.index(i)]
        if float(modscore) > float(score):
          type2 = "n"
      if i in itrksprots:
        itrscore = itrksscores[itrksprots.index(i)]
        if float(itrscore) > float(score):
          type2 = "n"
      if i in traksprots:
        trascore = traksscores[traksprots.index(i)]
        if float(trascore) > float(score):
          type2 = "n"
      if i in bt1fasprots:
        bt1fasscore = bt1fasscores[bt1fasprots.index(i)]
        if float(bt1fasscore) > float(score):
          type2 = "n"
      if i in ft1fasprots:
        ft1fasscore = ft1fasscores[ft1fasprots.index(i)]
        if float(ft1fasscore) > float(score):
          type2 = "n"
      if i in t2fasprots:
        t2fasscore = t2fasscores[t2fasprots.index(i)]
        if float(t2fasscore) > float(score):
          type2 = "n"
      if i in hgleprots:
        hglescore = hglescores[hgleprots.index(i)]
        if float(hglescore) > float(score):
          type2 = "n"
      if i in fabhprots:
        fabhscore = fabhscores[fabhprots.index(i)]
        if float(fabhscore) > float(score):
          type2 = "n"
      if type2 == "y" and i not in t2pksprots and i not in t1pksprots:
        t2pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Type II ketosynthase, model 2",t2ks2scores[t2ks2prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Type II ketosynthase, model 2",t2ks2scores[t2ks2prots.index(i)]]]
  #Extract type III PKS proteins
  t3pksprots = []
  if 1 in geneclustertypes or 2 in geneclustertypes or 3 in geneclustertypes or 4 in geneclustertypes:
    t3n = parsehmmoutput(63,hmmoutputfolder + "Chal_sti_synt_N.txt")
    t3c = parsehmmoutput(35,hmmoutputfolder + "Chal_sti_synt_C.txt")
    t3nprots = t3n[0]
    t3nscores = t3n[1]
    t3cprots = t3c[0]
    t3cscores = t3c[1]
    for i in t3cprots:
      if i not in t3pksprots and i not in t1pksprots and i not in t2pksprots:
        t3pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Chalcone/stilbene synthase,C-terminus",t3cscores[t3cprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Chalcone/stilbene synthase,C-terminus",t3cscores[t3cprots.index(i)]]]
    for i in t3nprots:
      if i not in t3pksprots and i not in t1pksprots and i not in t2pksprots:
        t3pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Chalcone/stilbene synthase,N-terminus",t3nscores[t3nprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Chalcone/stilbene synthase,N-terminus",t3nscores[t3nprots.index(i)]]]
  #Extract 'type IV' hglE-like PKS proteins, cut-off:50; only if not already scored as type 1-3 PKS, and not if FAS HMM has higher score
  t4pksprots = []
  if 1 in geneclustertypes or 2 in geneclustertypes or 3 in geneclustertypes or 4 in geneclustertypes:
    t2fasks = parsehmmoutput(50,hmmoutputfolder + "t2fas.txt")
    t2fasprots = t2fasks[0]
    t2fasscores = t2fasks[1]
    for i in hgleprots:
      type4 = "y"
      score = hglescores[hgleprots.index(i)]
      if i in bt1fasprots:
        bt1fasscore = bt1fasscores[bt1fasprots.index(i)]
        if float(bt1fasscore) > float(score):
          type4 = "n"
      if i in ft1fasprots:
        ft1fasscore = ft1fasscores[ft1fasprots.index(i)]
        if float(ft1fasscore) > float(score):
          type4 = "n"
      if i in t2fasprots:
        t2fasscore = t2fasscores[t2fasprots.index(i)]
        if float(t2fasscore) > float(score):
          type4 = "n"
      if i in fabhprots:
        fabhscore = fabhscores[fabhprots.index(i)]
        if float(fabhscore) > float(score):
          type4 = "n"
      if i not in t1pksprots and i not in t2pksprots and i not in t3pksprots and i not in transatpksprots and type4 == "y":
        t4pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Atypical PKS domain, HglE-like",hglescores[hgleprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Atypical PKS domain, HglE-like",hglescores[hgleprots.index(i)]]]
    for i in hgldprots:
      type4 = "y"
      score = hgldscores[hgldprots.index(i)]
      if i in bt1fasprots:
        bt1fasscore = bt1fasscores[bt1fasprots.index(i)]
        if float(bt1fasscore) > float(score):
          type4 = "n"
      if i in ft1fasprots:
        ft1fasscore = ft1fasscores[ft1fasprots.index(i)]
        if float(ft1fasscore) > float(score):
          type4 = "n"
      if i in t2fasprots:
        t2fasscore = t2fasscores[t2fasprots.index(i)]
        if float(t2fasscore) > float(score):
          type4 = "n"
      if i in fabhprots:
        fabhscore = fabhscores[fabhprots.index(i)]
        if float(fabhscore) > float(score):
          type4 = "n"
      if i not in t1pksprots and i not in t2pksprots and i not in t3pksprots and i not in transatpksprots and type4 == "y" and i not in t4pksprots:
        t4pksprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Atypical PKS domain, HglD-like",hgldscores[hgldprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Atypical PKS domain, HglD-like",hgldscores[hgldprots.index(i)]]]
  #Extract NRPS proteins, C cut-off: 20; A cut-off:20, both should be there, or single domain proteins C,A, or T should be within 20kb of each other or a full NRPS
  nrpsprots = []
  if 1 in geneclustertypes or 5 in geneclustertypes:
    cond = parsehmmoutput(20,hmmoutputfolder + "Condensation.txt")
    amp = parsehmmoutput(20,hmmoutputfolder + "AMP-binding.txt")
    ampox = parsehmmoutput(50,hmmoutputfolder + "A-OX.txt")
    ampoxprots = ampox[0]
    ampoxscores = ampox[1]
    for i in ampox[0]:
      if i not in amp:
        amp.append(i)
    cprots = cond[0]
    cscores = cond[1]
    aprots = amp[0]
    ascores = amp[1]
    nrpsprots = []
    for i in cprots:
      if i in aprots:
        nrpsprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Condensation domain",cscores[cprots.index(i)]])
          if i in aprots:
            detdomlist.append(["Adenylation domain",ascores[aprots.index(i)]])
          elif i in ampoxprots:
            detdomlist.append(["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in aprots:
            detecteddomainsdict[i] = [["Condensation domain",cscores[cprots.index(i)]],["Adenylation domain",ascores[aprots.index(i)]]]
          elif i in ampoxprots:
            detecteddomainsdict[i] = [["Condensation domain",cscores[cprots.index(i)]],["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]]]
    for i in t1pksprots:
      if i in aprots:
        nrpsprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in aprots:
            detdomlist.append(["Adenylation domain",ascores[aprots.index(i)]])
          elif i in ampoxprots:
            detdomlist.append(["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in aprots:
            detecteddomainsdict[i] = [["Adenylation domain",ascores[aprots.index(i)]]]
          elif i in ampoxprots:
            detecteddomainsdict[i] = [["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]]]
    single_aprots = []
    single_cprots = []
    single_pptprots = []
    pptprots = parsehmmoutput(20,hmmoutputfolder + "PP-binding.txt")[0]
    for i in aprots:
      if i not in nrpsprots:
        single_aprots.append(i)
    for i in cprots:
      if i not in nrpsprots:
        single_cprots.append(i)
    for i in pptprots:
      if i not in nrpsprots:
        single_pptprots.append(i)
    genelist = proteins[2]
    genedict = proteins[3]
    single_aprots_positions = {}
    single_cprots_positions = {}
    single_pptprots_positions = {}
    nrpsprots_positions = {}
    for j in single_aprots:
      if j in genelist:
        protstart_abs = min([int(genedict[j][0]),int(genedict[j][1])])
        protend_abs = max([int(genedict[j][0]),int(genedict[j][1])])
        single_aprots_positions[j] = int((protend_abs + protstart_abs) / 2)
    for j in single_cprots:
      if j in genelist:
        protstart_abs = min([int(genedict[j][0]),int(genedict[j][1])])
        protend_abs = max([int(genedict[j][0]),int(genedict[j][1])])
        single_cprots_positions[j] = int((protend_abs + protstart_abs) / 2)
    for j in single_pptprots:
      if j in genelist:
        protstart_abs = min([int(genedict[j][0]),int(genedict[j][1])])
        protend_abs = max([int(genedict[j][0]),int(genedict[j][1])])
        single_pptprots_positions[j] = int((protend_abs + protstart_abs) / 2)
    for j in nrpsprots:
      if j in genelist:
        protstart_abs = min([int(genedict[j][0]),int(genedict[j][1])])
        protend_abs = max([int(genedict[j][0]),int(genedict[j][1])])
        nrpsprots_positions[j] = int((protend_abs + protstart_abs) / 2)
    nrpsprots2 = []
    for i in nrpsprots:
      nrpsprots2.append(i)
    for j in single_aprots:
      include = "n"
      pos = single_aprots_positions[j]
      for i in single_cprots:
        pos2 = single_cprots_positions[i]
        if abs(pos - pos2) < 20000:
          include = "y"
      for i in nrpsprots2:
        pos2 = nrpsprots_positions[i]
        if abs(pos - pos2) < 20000:
          include = "y"
      if include == "y":
        nrpsprots.append(j)
        if detecteddomainsdict.has_key(j):
          detdomlist = detecteddomainsdict[j]
          if j in aprots:
            detdomlist.append(["Adenylation domain",ascores[aprots.index(j)]])
          elif j in ampoxprots:
            detdomlist.append(["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(j)]])
          detecteddomainsdict[j] = detdomlist
        else:
          if j in aprots:
            detecteddomainsdict[j] = [["Adenylation domain",ascores[aprots.index(j)]]]
          elif j in ampoxprots:
            detecteddomainsdict[j] = [["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(j)]]]
    for j in single_cprots:
      include = "n"
      pos = single_cprots_positions[j]
      for i in single_aprots:
        pos2 = single_aprots_positions[i]
        if abs(pos - pos2) < 20000:
          include = "y"
      for i in nrpsprots2:
        pos2 = nrpsprots_positions[i]
        if abs(pos - pos2) < 20000:
          include = "y"
      if include == "y":
        nrpsprots.append(j)
        if detecteddomainsdict.has_key(j):
          detdomlist = detecteddomainsdict[j]
          detdomlist.append(["Condensation domain",cscores[cprots.index(j)]])
          detecteddomainsdict[j] = detdomlist
        else:
          detecteddomainsdict[j] = [["Condensation domain",cscores[cprots.index(j)]]]
  #Extract Terpene synthase proteins, various cut-offs
  terpeneprots = []
  if 1 in geneclustertypes or 6 in geneclustertypes:
    terpene = parsehmmoutput(23,hmmoutputfolder + "Terpene_synth_C.txt")
    terpeneprots = terpene[0]
    terpenescores = terpene[1]
    for i in terpeneprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["Terpene synthase, C-terminus",terpenescores[terpeneprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["Terpene synthase, C-terminus",terpenescores[terpeneprots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    physqualdata = parsehmmoutput(20,hmmoutputfolder + "phytoene_synt.txt")
    physqualprots = physqualdata[0]
    physqualscores = physqualdata[1]
    for i in physqualprots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Phytoene/squalene synthase",physqualscores[physqualprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Phytoene/squalene synthase",physqualscores[physqualprots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    lycopenedata = parsehmmoutput(80,hmmoutputfolder + "Lycopene_cycl.txt")
    lycopeneprots = lycopenedata[0]
    lycopenescores = lycopenedata[1]
    for i in lycopeneprots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Lycopene cyclase",lycopenescores[lycopeneprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Lycopene cyclase",lycopenescores[lycopeneprots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    terpene_cyclasesdata = parsehmmoutput(50,hmmoutputfolder + "terpene_cyclase.txt")
    terpene_cyclases = terpene_cyclasesdata[0]
    terpene_cyclases_scores = terpene_cyclasesdata[1]
    for i in terpene_cyclases:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Terpene cyclase",terpene_cyclases_scores[terpene_cyclases.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Terpene cyclase",terpene_cyclases_scores[terpene_cyclases.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    NapT7 = parsehmmoutput(250,hmmoutputfolder + "NapT7.txt")
    NapT7prots = NapT7[0]
    NapT7scores = NapT7[1]
    for i in NapT7prots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["NapT7",NapT7scores[NapT7prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["NapT7",NapT7scores[NapT7prots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    fung_ggpps = parsehmmoutput(420,hmmoutputfolder + "fung_ggpps.txt")
    fung_ggppsprots = fung_ggpps[0]
    fung_ggppsscores = fung_ggpps[1]
    for i in fung_ggppsprots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Fungal geranylgeranyl pyrophosphate synthase, model 1",fung_ggppsscores[fung_ggppsprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Fungal geranylgeranyl pyrophosphate synthase, model 1",fung_ggppsscores[fung_ggppsprots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    fung_ggpps2 = parsehmmoutput(312,hmmoutputfolder + "fung_ggpps2.txt")
    fung_ggpps2prots = fung_ggpps2[0]
    fung_ggpps2scores = fung_ggpps2[1]
    for i in fung_ggpps2prots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Fungal geranylgeranyl pyrophosphate synthase, model 2",fung_ggpps2scores[fung_ggpps2prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Fungal geranylgeranyl pyrophosphate synthase, model 2",fung_ggpps2scores[fung_ggpps2prots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    dmat = parsehmmoutput(200,hmmoutputfolder + "dmat.txt")
    dmatprots = dmat[0]
    dmatscores = dmat[1]
    for i in dmatprots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Dimethylallyl tryptophan synthase",dmatscores[dmatprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Dimethylallyl tryptophan synthase",dmatscores[dmatprots.index(i)]]]
  if 1 in geneclustertypes or 6 in geneclustertypes:
    trichodiene_synth = parsehmmoutput(150,hmmoutputfolder + "trichodiene_synth.txt")
    trichodiene_synthprots = trichodiene_synth[0]
    trichodiene_synthscores = trichodiene_synth[1]
    for i in trichodiene_synthprots:
      if i not in terpeneprots:
        terpeneprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Trichodiene synthase",trichodiene_synthscores[trichodiene_synthprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Trichodiene synthase",trichodiene_synthscores[trichodiene_synthprots.index(i)]]]
  #Extract lantibiotic proteins, LanC cut-off: 80, Lant_dehN & Lant_dehC combination cut-off: 20 each
  lantprots = []
  if 1 in geneclustertypes or 7 in geneclustertypes:
    lantc = parsehmmoutput(80,hmmoutputfolder + "LANC_like.txt")
    lancprots = lantc[0]
    lancscores = lantc[1]
    landehn = parsehmmoutput(20,hmmoutputfolder + "Lant_dehyd_N.txt")
    landehnprots = landehn[0]
    landehnscores = landehn[1]
    landehc = parsehmmoutput(20,hmmoutputfolder + "Lant_dehyd_C.txt")
    landehcprots = landehc[0]
    landehcscores = landehc[1]
    lanti1 = parsehmmoutput(20,hmmoutputfolder + "Antimicrobial18.txt")
    lanti1prots = lanti1[0]
    lanti1scores = lanti1[1]
    lanti2 = parsehmmoutput(20,hmmoutputfolder + "Gallidermin.txt")
    lanti2prots = lanti2[0]
    lanti2scores = lanti2[1]
    lanti3 = parsehmmoutput(20,hmmoutputfolder + "L_biotic_typeA.txt")
    lanti3prots = lanti3[0]
    lanti3scores = lanti3[1]
    lanti4 = parsehmmoutput(20,hmmoutputfolder + "LE-DUF.txt")
    lanti4prots = lanti4[0]
    lanti4scores = lanti4[1]
    lanti5 = parsehmmoutput(20,hmmoutputfolder + "LE-LAC481.txt")
    lanti5prots = lanti5[0]
    lanti5scores = lanti5[1]
    lanti6 = parsehmmoutput(20,hmmoutputfolder + "LE-LanBC.txt")
    lanti6prots = lanti6[0]
    lanti6scores = lanti6[1]
    lanti7 = parsehmmoutput(20,hmmoutputfolder + "LE-MER+2PEP.txt")
    lanti7prots = lanti7[0]
    lanti7scores = lanti7[1]
    lanti8 = parsehmmoutput(20,hmmoutputfolder + "MA-2PEPA.txt")
    lanti8prots = lanti8[0]
    lanti8scores = lanti8[1]
    lanti9 = parsehmmoutput(20,hmmoutputfolder + "MA-DUF.txt")
    lanti9prots = lanti9[0]
    lanti9scores = lanti9[1]
    lanti10 = parsehmmoutput(20,hmmoutputfolder + "MA-EPI.txt")
    lanti10prots = lanti10[0]
    lanti10scores = lanti10[1]
    lanti11 = parsehmmoutput(20,hmmoutputfolder + "MA-LAC481.txt")
    lanti11prots = lanti11[0]
    lanti11scores = lanti11[1]
    lanti12 = parsehmmoutput(20,hmmoutputfolder + "MA-NIS+EPI.txt")
    lanti12prots = lanti12[0]
    lanti12scores = lanti12[1]
    lanti13 = parsehmmoutput(20,hmmoutputfolder + "MA-NIS.txt")
    lanti13prots = lanti13[0]
    lanti13scores = lanti13[1]
    lanti14 = parsehmmoutput(18,hmmoutputfolder + "TIGR03731.txt")
    lanti14prots = lanti14[0]
    lanti14scores = lanti14[1]
    lantiprots = lanti1prots + lanti2prots + lanti3prots + lanti4prots + lanti5prots + lanti6prots + lanti7prots + lanti8prots + lanti9prots + lanti10prots + lanti11prots + lanti12prots + lanti13prots + lanti14prots
    lantiprots2 = []
    for i in lantiprots:
      if i not in lantiprots2:
        lantiprots2.append(i)
    lantiprots = lantiprots2
    for i in lancprots:
      lantprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["LanC lanthionine synthase domain",lancscores[lancprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["LanC lanthionine synthase domain",lancscores[lancprots.index(i)]]]
    for i in landehnprots:
      if i in landehcprots and i not in lantprots:
        lantprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Lantibiotic dehydratase, N-terminus",landehnscores[landehnprots.index(i)]])
          detdomlist.append(["Lantibiotic dehydratase, C-terminus",landehcscores[landehcprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Lantibiotic dehydratase, N-terminus",landehnscores[landehnprots.index(i)]],["Lantibiotic dehydratase, C-terminus",landehcscores[landehcprots.index(i)]]]
    for i in lantiprots:
      if i not in lantprots:
        lantprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti1prots:
            detdomlist.append(["Antimicrobial18 domain",lanti1scores[lanti1prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti1prots:
            detecteddomainsdict[i] = [["Antimicrobial18 domain",lanti1scores[lanti1prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti2prots:
            detdomlist.append(["Gallidermin domain",lanti2scores[lanti2prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti2prots:
            detecteddomainsdict[i] = [["Gallidermin domain",lanti2scores[lanti2prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti3prots:
            detdomlist.append(["L_biotic_typeA domain",lanti3scores[lanti3prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti3prots:
            detecteddomainsdict[i] = [["L_biotic_typeA domain",lanti3scores[lanti3prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti4prots:
            detdomlist.append(["LE-DUF domain",lanti4scores[lanti4prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti4prots:
            detecteddomainsdict[i] = [["LE-DUF domain",lanti4scores[lanti4prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti5prots:
            detdomlist.append(["LE-LAC481 domain",lanti5scores[lanti5prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti5prots:
            detecteddomainsdict[i] = [["LE-LAC481 domain",lanti5scores[lanti5prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti6prots:
            detdomlist.append(["LE-LanBC domain",lanti6scores[lanti6prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti6prots:
            detecteddomainsdict[i] = [["LE-LanBC domain",lanti6scores[lanti6prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti7prots:
            detdomlist.append(["LE-MER+2PEP domain",lanti7scores[lanti7prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti7prots:
            detecteddomainsdict[i] = [["LE-MER+2PEP domain",lanti7scores[lanti7prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti8prots:
            detdomlist.append(["MA-2PEPA domain",lanti8scores[lanti8prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti8prots:
            detecteddomainsdict[i] = [["MA-2PEPA domain",lanti8scores[lanti8prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti9prots:
            detdomlist.append(["MA-DUF domain",lanti9scores[lanti9prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti9prots:
            detecteddomainsdict[i] = [["MA-DUF domain",lanti9scores[lanti9prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti10prots:
            detdomlist.append(["MA-EPI domain",lanti10scores[lanti10prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti10prots:
            detecteddomainsdict[i] = [["MA-EPI domain",lanti10scores[lanti10prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti11prots:
            detdomlist.append(["MA-LAC481 domain",lanti11scores[lanti11prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti11prots:
            detecteddomainsdict[i] = [["MA-LAC481 domain",lanti11scores[lanti11prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti12prots:
            detdomlist.append(["MA-NIS+EPI domain",lanti12scores[lanti12prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti12prots:
            detecteddomainsdict[i] = [["MA-NIS+EPI domain",lanti12scores[lanti12prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti13prots:
            detdomlist.append(["MA-NIS domain",lanti13scores[lanti13prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti13prots:
            detecteddomainsdict[i] = [["MA-NIS domain",lanti13scores[lanti13prots.index(i)]]]
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          if i in lanti14prots:
            detdomlist.append(["TIGR03731: lantibiotic, gallidermin/nisin family",lanti14scores[lanti14prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in lanti14prots:
            detecteddomainsdict[i] = [["TIGR03731: lantibiotic, gallidermin/nisin family",lanti14scores[lanti14prots.index(i)]]]
  #Bacteriocin proteins, various cut-offs
  bcinprots = []
  if 1 in geneclustertypes or 8 in geneclustertypes:
    bcin1prots = parsehmmoutput(50,hmmoutputfolder + "strepbact.txt")[0]
    bcin2prots = parsehmmoutput(90,hmmoutputfolder + "Antimicrobial14.txt")[0]
    bcin3prots = parsehmmoutput(23,hmmoutputfolder + "Bacteriocin_IId.txt")[0]
    bcin4prots = parsehmmoutput(92,hmmoutputfolder + "BacteriocIIc_cy.txt")[0]
    bcin5prots = parsehmmoutput(40,hmmoutputfolder + "Bacteriocin_II.txt")[0]
    bcin6prots = parsehmmoutput(24,hmmoutputfolder + "Lactococcin.txt")[0]
    bcin7prots = parsehmmoutput(31,hmmoutputfolder + "Antimicrobial17.txt")[0]
    bcin8prots = parsehmmoutput(25,hmmoutputfolder + "Lactococcin_972.txt")[0]
    bcin9prots = parsehmmoutput(27,hmmoutputfolder + "Bacteriocin_IIc.txt")[0]
    bcin10prots = parsehmmoutput(78,hmmoutputfolder + "LcnG-beta.txt")[0]
    bcin11prots = parsehmmoutput(56,hmmoutputfolder + "Bacteriocin_IIi.txt")[0]
    bcin12prots = parsehmmoutput(98,hmmoutputfolder + "Subtilosin_A.txt")[0]
    bcin13prots = parsehmmoutput(27,hmmoutputfolder + "Cloacin.txt")[0]
    bcin14prots = parsehmmoutput(25,hmmoutputfolder + "Linocin_M18.txt")[0]
    bcin15prots = parsehmmoutput(150,hmmoutputfolder + "TIGR03603.txt")[0]
    bcin16prots = parsehmmoutput(440,hmmoutputfolder + "TIGR03604.txt")[0]
    bcin17prots = parsehmmoutput(200,hmmoutputfolder + "TIGR03605.txt")[0]
    bcin18prots = parsehmmoutput(18,hmmoutputfolder + "TIGR03651.txt")[0]
    bcin19prots = parsehmmoutput(35,hmmoutputfolder + "TIGR03678.txt")[0]
    bcin20prots = parsehmmoutput(400,hmmoutputfolder + "TIGR03693.txt")[0]
    bcin21prots = parsehmmoutput(16,hmmoutputfolder + "TIGR03798.txt")[0]
    bcin22prots = parsehmmoutput(150,hmmoutputfolder + "TIGR03882.txt")[0]
    bcin23prots = parsehmmoutput(50,hmmoutputfolder + "TIGR03601.txt")[0]
    bcin24prots = parsehmmoutput(50,hmmoutputfolder + "TIGR03602.txt")[0]
    bcin25prots = parsehmmoutput(20,hmmoutputfolder + "mvnA.txt")[0]
    bcin26prots = parsehmmoutput(20,hmmoutputfolder + "thiostrepton.txt")[0]
    bcin1scores = parsehmmoutput(50,hmmoutputfolder + "strepbact.txt")[1]
    bcin2scores = parsehmmoutput(90,hmmoutputfolder + "Antimicrobial14.txt")[1]
    bcin3scores = parsehmmoutput(23,hmmoutputfolder + "Bacteriocin_IId.txt")[1]
    bcin4scores = parsehmmoutput(92,hmmoutputfolder + "BacteriocIIc_cy.txt")[1]
    bcin5scores = parsehmmoutput(40,hmmoutputfolder + "Bacteriocin_II.txt")[1]
    bcin6scores = parsehmmoutput(24,hmmoutputfolder + "Lactococcin.txt")[1]
    bcin7scores = parsehmmoutput(31,hmmoutputfolder + "Antimicrobial17.txt")[1]
    bcin8scores = parsehmmoutput(25,hmmoutputfolder + "Lactococcin_972.txt")[1]
    bcin9scores = parsehmmoutput(27,hmmoutputfolder + "Bacteriocin_IIc.txt")[1]
    bcin10scores = parsehmmoutput(78,hmmoutputfolder + "LcnG-beta.txt")[1]
    bcin11scores = parsehmmoutput(56,hmmoutputfolder + "Bacteriocin_IIi.txt")[1]
    bcin12scores = parsehmmoutput(98,hmmoutputfolder + "Subtilosin_A.txt")[1]
    bcin13scores = parsehmmoutput(27,hmmoutputfolder + "Cloacin.txt")[1]
    bcin14scores = parsehmmoutput(25,hmmoutputfolder + "Linocin_M18.txt")[1]
    bcin15scores = parsehmmoutput(150,hmmoutputfolder + "TIGR03603.txt")[1]
    bcin16scores = parsehmmoutput(440,hmmoutputfolder + "TIGR03604.txt")[1]
    bcin17scores = parsehmmoutput(200,hmmoutputfolder + "TIGR03605.txt")[1]
    bcin18scores = parsehmmoutput(18,hmmoutputfolder + "TIGR03651.txt")[1]
    bcin19scores = parsehmmoutput(35,hmmoutputfolder + "TIGR03678.txt")[1]
    bcin20scores = parsehmmoutput(400,hmmoutputfolder + "TIGR03693.txt")[1]
    bcin21scores = parsehmmoutput(16,hmmoutputfolder + "TIGR03798.txt")[1]
    bcin22scores = parsehmmoutput(150,hmmoutputfolder + "TIGR03882.txt")[1]
    bcin23scores = parsehmmoutput(50,hmmoutputfolder + "TIGR03601.txt")[1]
    bcin24scores = parsehmmoutput(50,hmmoutputfolder + "TIGR03602.txt")[1]
    bcin25scores = parsehmmoutput(20,hmmoutputfolder + "mvnA.txt")[1]
    bcin26scores = parsehmmoutput(20,hmmoutputfolder + "thiostrepton.txt")[1]
    bcinprots = bcin1prots + bcin2prots + bcin3prots + bcin4prots + bcin5prots + bcin6prots + bcin7prots + bcin8prots + bcin9prots + bcin10prots + bcin11prots + bcin12prots + bcin13prots + bcin14prots + bcin15prots + bcin16prots + bcin17prots + bcin18prots + bcin19prots + bcin20prots + bcin21prots + bcin22prots + bcin23prots + bcin24prots + bcin25prots + bcin26prots
    bcinprots2 = []
    for i in bcinprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin1prots:
          detdomlist.append(["Putative Streptomyces bacteriocin",bcin1scores[bcin1prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin1prots:
          detecteddomainsdict[i] = [["Putative Streptomyces bacteriocin",bcin1scores[bcin1prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin2prots:
          detdomlist.append(["Antimicrobial14 domain",bcin2scores[bcin2prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin2prots:
          detecteddomainsdict[i] = [["Antimicrobial14 domain",bcin2scores[bcin2prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin3prots:
          detdomlist.append(["Bacteriocin_IId domain",bcin3scores[bcin3prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin3prots:
          detecteddomainsdict[i] = [["Bacteriocin_IId domain",bcin3scores[bcin3prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin4prots:
          detdomlist.append(["BacteriocIIc_cy domain",bcin4scores[bcin4prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin4prots:
          detecteddomainsdict[i] = [["BacteriocIIc_cy domain",bcin4scores[bcin4prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin5prots:
          detdomlist.append(["Bacteriocin_II domain",bcin5scores[bcin5prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin5prots:
          detecteddomainsdict[i] = [["Bacteriocin_II domain",bcin5scores[bcin5prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin6prots:
          detdomlist.append(["Lactococcin",bcin6scores[bcin6prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin6prots:
          detecteddomainsdict[i] = [["Lactococcin",bcin6scores[bcin6prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin7prots:
          detdomlist.append(["Antimicrobial17 domain",bcin7scores[bcin7prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin7prots:
          detecteddomainsdict[i] = [["Antimicrobial17 domain",bcin7scores[bcin7prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin8prots:
          detdomlist.append(["Lactococcin_972 domain",bcin8scores[bcin8prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin8prots:
          detecteddomainsdict[i] = [["Lactococcin_972 domain",bcin8scores[bcin8prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin9prots:
          detdomlist.append(["Bacteriocin_IIc domain",bcin9scores[bcin9prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin9prots:
          detecteddomainsdict[i] = [["Bacteriocin_IIc domain",bcin9scores[bcin9prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin10prots:
          detdomlist.append(["LcnG-beta domain",bcin10scores[bcin10prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin10prots:
          detecteddomainsdict[i] = [["LcnG-beta domain",bcin10scores[bcin10prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin11prots:
          detdomlist.append(["Bacteriocin_IIi domain",bcin11scores[bcin11prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin11prots:
          detecteddomainsdict[i] = [["Bacteriocin_IIi domain",bcin11scores[bcin11prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin12prots:
          detdomlist.append(["Subtilosin_A domain",bcin12scores[bcin12prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin12prots:
          detecteddomainsdict[i] = [["Subtilosin_A domain",bcin12scores[bcin12prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin13prots:
          detdomlist.append(["Cloacin domain",bcin13scores[bcin13prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin13prots:
          detecteddomainsdict[i] = [["Cloacin domain",bcin13scores[bcin13prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin14prots:
          detdomlist.append(["Linocin_M18 domain",bcin14scores[bcin14prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin14prots:
          detecteddomainsdict[i] = [["Linocin_M18 domain",bcin14scores[bcin14prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin15prots:
          detdomlist.append(["TIGR03603: bacteriocin biosynthesis cyclodehydratase",bcin15scores[bcin15prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin15prots:
          detecteddomainsdict[i] = [["TIGR03603: bacteriocin biosynthesis cyclodehydratase",bcin15scores[bcin15prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin16prots:
          detdomlist.append(["TGIR03604: bacteriocin biosynthesis docking scaffold",bcin16scores[bcin16prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin16prots:
          detecteddomainsdict[i] = [["TGIR03604: bacteriocin biosynthesis docking scaffold",bcin16scores[bcin16prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin17prots:
          detdomlist.append(["TGIR03605: SagB-type dehydrogenase",bcin17scores[bcin17prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin17prots:
          detecteddomainsdict[i] = [["TGIR03605: SagB-type dehydrogenase",bcin17scores[bcin17prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin18prots:
          detdomlist.append(["TIGR03651: bacteriocin, circularin A/uberolysin family",bcin18scores[bcin18prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin18prots:
          detecteddomainsdict[i] = [["TIGR03651: bacteriocin, circularin A/uberolysin family",bcin18scores[bcin18prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin19prots:
          detdomlist.append(["TIGR03678: bacteriocin, microcyclamide/patellamide family",bcin19scores[bcin19prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin19prots:
          detecteddomainsdict[i] = [["TIGR03678: bacteriocin, microcyclamide/patellamide family",bcin19scores[bcin19prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin20prots:
          detdomlist.append(["TIGR03693: thiazole-containing bacteriocin maturation protein",bcin20scores[bcin20prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin20prots:
          detecteddomainsdict[i] = [["TIGR03693: thiazole-containing bacteriocin maturation protein",bcin20scores[bcin20prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin21prots:
          detdomlist.append(["TIGR03798: bacteriocin propeptide",bcin21scores[bcin21prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin21prots:
          detecteddomainsdict[i] = [["TIGR03798: bacteriocin propeptide",bcin21scores[bcin21prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin22prots:
          detdomlist.append(["TIGR03882: bacteriocin biosynthesis cyclodehydratase",bcin22scores[bcin22prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin22prots:
          detecteddomainsdict[i] = [["TIGR03882: bacteriocin biosynthesis cyclodehydratase",bcin22scores[bcin22prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin23prots:
          detdomlist.append(["TIGR03601: bacteriocin, BA_2677 family",bcin23scores[bcin23prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin23prots:
          detecteddomainsdict[i] = [["TIGR03601: bacteriocin, BA_2677 family",bcin23scores[bcin23prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin24prots:
          detdomlist.append(["TIGR03602: bacteriocin protoxin, streptolysin S family",bcin24scores[bcin24prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin24prots:
          detecteddomainsdict[i] = [["TIGR03602: bacteriocin protoxin, streptolysin S family",bcin24scores[bcin24prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin25prots:
          detdomlist.append(["Bacteriocin, microviridin family",bcin25scores[bcin25prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin25prots:
          detecteddomainsdict[i] = [["Bacteriocin, microviridin family",bcin25scores[bcin25prots.index(i)]]]
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        if i in bcin26prots:
          detdomlist.append(["Thiopeptide, thiostrepton-like",bcin26scores[bcin26prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        if i in bcin26prots:
          detecteddomainsdict[i] = [["Thiopeptide, thiostrepton-like",bcin26scores[bcin26prots.index(i)]]]
      if i not in bcinprots2:
        bcinprots2.append(i)
    bcinprots = bcinprots2
  #Extract beta-lactam synthetase proteins, cut-off: 250
  lactamprots = []
  if 1 in geneclustertypes or 9 in geneclustertypes:
    bls = parsehmmoutput(250,hmmoutputfolder + "BLS.txt")
    blsprots = bls[0]
    blsscores = bls[1]
    for i in bls[0]:
      lactamprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["Beta-lactam synthase",blsscores[blsprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["Beta-lactam synthase",blsscores[blsprots.index(i)]]]
    cas = parsehmmoutput(250,hmmoutputfolder + "CAS.txt")
    casprots = cas[0]
    casscores = cas[1]
    for i in cas[0]:
      if i not in lactamprots:
        lactamprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Clavulanic acid synthase-like",casscores[casprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Clavulanic acid synthase-like",casscores[casprots.index(i)]]]
    tabtoxin = parsehmmoutput(500,hmmoutputfolder + "tabtoxin.txt")
    tabtoxinprots = tabtoxin[0]
    tabtoxinscores = tabtoxin[1]
    for i in tabtoxin[0]:
      if i not in lactamprots:
        lactamprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Tabtoxin synthase-like",tabtoxinscores[tabtoxinprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Tabtoxin synthase-like",tabtoxinscores[tabtoxinprots.index(i)]]]
  #Extract aminoglycoside / aminocyclitol biosynthesis clusters, clusters taken from Flatt & Mahmud et al. 2007
  amglyccyclprots = []
  if 1 in geneclustertypes or 10 in geneclustertypes:
    strH = parsehmmoutput(200,hmmoutputfolder + "strH_like.txt")
    strhprots = strH[0]
    strhscores = strH[1]
    for i in strH[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["StrH-like glycosyltransferase",strhscores[strhprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["StrH-like glycosyltransferase",strhscores[strhprots.index(i)]]]
    strK1 = parsehmmoutput(800,hmmoutputfolder + "strK_like1.txt")
    strk1prots = strK1[0]
    strk1scores = strK1[1]
    for i in strK1[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["StrK-like phosphatase",strk1scores[strk1prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["StrK-like phosphatase",strk1scores[strk1prots.index(i)]]]
    strK2 = parsehmmoutput(650,hmmoutputfolder + "strK_like2.txt")
    strk2prots = strK2[0]
    strk2scores = strK2[1]
    for i in strK2[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["StrK-like phosphatase, model 2",strk2scores[strk2prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["StrK-like phosphatase, model 2",strk2scores[strk2prots.index(i)]]]
    neoL = parsehmmoutput(50,hmmoutputfolder + "neoL_like.txt")
    neolprots = neoL[0]
    neolscores = neoL[1]
    for i in neoL[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["NeoL-like deacetylase",neolscores[neolprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["NeoL-like deacetylase",neolscores[neolprots.index(i)]]]
    DOIS = parsehmmoutput(500,hmmoutputfolder + "DOIS.txt")
    doisprots = DOIS[0]
    doisscores = DOIS[1]
    for i in DOIS[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["2-deoxy-scyllo-inosose synthase",doisscores[doisprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["2-deoxy-scyllo-inosose synthase",doisscores[doisprots.index(i)]]]
    valA = parsehmmoutput(600,hmmoutputfolder + "valA_like.txt")
    valaprots = valA[0]
    valascores = valA[1]
    for i in valA[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["2-epi-5-epi-valiolone synthase, ValA-like",valascores[valaprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["2-epi-5-epi-valiolone synthase, ValA-like",valascores[valaprots.index(i)]]]
    spcFG = parsehmmoutput(200,hmmoutputfolder + "spcFG_like.txt")
    spcfgprots = spcFG[0]
    spcfgscores = spcFG[1]
    for i in spcFG[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["SpcF/SpcG-like glycosyltransferase",spcfgscores[spcfgprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["SpcF/SpcG-like glycosyltransferase",spcfgscores[spcfgprots.index(i)]]]
    spcDK_glyc = parsehmmoutput(600,hmmoutputfolder + "spcDK_like_glyc.txt")
    spcdkglycprots = spcDK_glyc[0]
    spcdkglycscores = spcDK_glyc[1]
    for i in spcDK_glyc[0]:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["SpcD/SpcK-like thymidylyltransferase",spcdkglycscores[spcdkglycprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["SpcD/SpcK-like thymidylyltransferase",spcdkglycscores[spcdkglycprots.index(i)]]]
    salQ = parsehmmoutput(480,hmmoutputfolder + "salQ.txt")
    salqprots = salQ[0]
    salqscores = salQ[1]
    for i in salqprots:
      amglyccyclprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["2-epi-5-epi-valiolone synthase, SalQ-like",salqscores[salqprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["2-epi-5-epi-valiolone synthase, SalQ-like",salqscores[salqprots.index(i)]]]
  #Extract aminocoumarin biosynthesis clusters
  aminocoumarinprots = []
  if 1 in geneclustertypes or 11 in geneclustertypes:
    novK = parsehmmoutput(200,hmmoutputfolder + "novK.txt")
    novkprots = novK[0]
    novkscores = novK[1]
    for i in novkprots:
      aminocoumarinprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["NovK-like reductase",novkscores[novkprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["NovK-like reductase",novkscores[novkprots.index(i)]]]
    novJ = parsehmmoutput(350,hmmoutputfolder + "novJ.txt")
    novjprots = novJ[0]
    novjscores = novJ[1]
    for i in novjprots:
      aminocoumarinprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["NovJ-like reductase",novjscores[novjprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["NovJ-like reductase",novjscores[novjprots.index(i)]]]
    novI = parsehmmoutput(600,hmmoutputfolder + "novI.txt")
    noviprots = novI[0]
    noviscores = novI[1]
    for i in noviprots :
      aminocoumarinprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["NovI-like cytochrome P450",noviscores[noviprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["NovI-like cytochrome P450",noviscores[noviprots.index(i)]]]
    novH = parsehmmoutput(750,hmmoutputfolder + "novH.txt")
    novhprots = novH[0]
    novhscores = novH[1]
    for i in novhprots:
      aminocoumarinprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["NovH-like protein",novhscores[novhprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["NovH-like protein",novhscores[novhprots.index(i)]]]
    spcDK_like_cou = parsehmmoutput(600,hmmoutputfolder + "spcDK_like_cou.txt")
    spcDK_like_cou_prots = spcDK_like_cou[0]
    spcDK_like_cou_scores = spcDK_like_cou[1]
    for i in spcDK_like_cou_prots:
      aminocoumarinprots.append(i)
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["SpcD/SpcK-like thymidylyltransferase, aminocoumarins group",spcDK_like_cou_scores[spcDK_like_cou_prots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["SpcD/SpcK-like thymidylyltransferase, aminocoumarins group",spcDK_like_cou_scores[spcDK_like_cou_prots.index(i)]]]
  #Extract siderophores biosynthesis proteins, IucA/C and AlcB
  siderophoreprots = []
  if 1 in geneclustertypes or 12 in geneclustertypes:
    siderophore = parsehmmoutput(30,hmmoutputfolder + "IucA_IucC.txt")
    siderophoreprots = siderophore[0]
    siderophorescores = siderophore[1]
    for i in siderophoreprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["IucA-IucC domain",siderophorescores[siderophoreprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["IucA-IucC domain",siderophorescores[siderophoreprots.index(i)]]]
  #Extract ectoine biosynthesis proteins
  ectprots = []
  if 1 in geneclustertypes or 13 in geneclustertypes:
    ect = parsehmmoutput(35,hmmoutputfolder + "ectoine_synt.txt")
    ectprots = ect[0]
    ectscores = ect[1]
    for i in ectprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["Ectoine synthase",ectscores[ectprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["Ectoine synthase",ectscores[ectprots.index(i)]]]
  #Extract butyrolactone biosynthesis proteins
  butyrprots = []
  if 1 in geneclustertypes or 14 in geneclustertypes:
    butyr= parsehmmoutput(25,hmmoutputfolder + "AfsA.txt")
    butyrprots = butyr[0]
    butyrscores = butyr[1]
    for i in butyrprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["AfsA butyrolactone synthesis domain",butyrscores[butyrprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["AfsA butyrolactone synthesis domain",butyrscores[butyrprots.index(i)]]]
  #Extract indole biosynthesis proteins
  indoleprots = []
  if 1 in geneclustertypes or 15 in geneclustertypes:
    indole = parsehmmoutput(100,hmmoutputfolder + "indsynth.txt")
    indoleprots = indole[0]
    indolescores = indole[1]
    for i in indoleprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["StaD-like chromopyrrolic acid synthase domain",indolescores[indoleprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["StaD-like chromopyrrolic acid synthase domain",indolescores[indoleprots.index(i)]]]
  #Extract nucleoside antibiotic biosynthesis proteins
  nucleoprots = []
  if 1 in geneclustertypes or 16 in geneclustertypes:
    nucleoprots = []
    lipm = parsehmmoutput(50,hmmoutputfolder + "LipM.txt")
    lipmprots = lipm[0]
    lipmscores = lipm[1]
    lipu = parsehmmoutput(30,hmmoutputfolder + "LipU.txt")
    lipuprots = lipu[0]
    lipuscores = lipu[1]
    lipv = parsehmmoutput(375,hmmoutputfolder + "LipV.txt")
    lipvprots = lipv[0]
    lipvscores = lipv[1]
    toyb = parsehmmoutput(175,hmmoutputfolder + "ToyB.txt")
    toybprots = toyb[0]
    toybscores = toyb[1]
    tund = parsehmmoutput(200,hmmoutputfolder + "TunD.txt")
    tundprots = tund[0]
    tundscores = tund[1]
    pur6 = parsehmmoutput(200,hmmoutputfolder + "pur6.txt")
    pur6prots = pur6[0]
    pur6scores = pur6[1]
    pur10 = parsehmmoutput(600,hmmoutputfolder + "pur10.txt")
    pur10prots = pur10[0]
    pur10scores = pur10[1]
    nikj = parsehmmoutput(200,hmmoutputfolder + "nikJ.txt")
    nikjprots = nikj[0]
    nikjscores = nikj[1]
    niko = parsehmmoutput(400,hmmoutputfolder + "nikO.txt")
    nikoprots = niko[0]
    nikoscores = niko[1]
    for i in lipmprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["LipM-like nucleotidyltransferase",lipmscores[lipmprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["LipM-like nucleotidyltransferase",lipmscores[lipmprots.index(i)]]]
    for i in lipuprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["LipU-like protein",lipuscores[lipuprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["LipU-like protein",lipuscores[lipuprots.index(i)]]]
    for i in lipvprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["LipV-like dehydrogenase",lipvscores[lipvprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["LipV-like dehydrogenase",lipvscores[lipvprots.index(i)]]]
    for i in toybprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["ToyB-like synthase",toybscores[toybprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["ToyB-like synthase",toybscores[toybprots.index(i)]]]
    for i in tundprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["TunD-like putative N-acetylglucosamine transferase",tundscores[tundprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["TunD-like putative N-acetylglucosamine transferase",tundscores[tundprots.index(i)]]]
    for i in pur6prots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Pur6-like synthetase",pur6scores[pur6prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Pur6-like synthetase",pur6scores[pur6prots.index(i)]]]
    for i in pur10prots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Pur10-like oxidoreductase",pur10scores[pur10prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Pur10-like oxidoreductase",pur10scores[pur10prots.index(i)]]]
    for i in nikjprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["NikJ-like protein",nikjscores[nikjprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["NikJ-like protein",nikjscores[nikjprots.index(i)]]]
    for i in nikoprots:
      if i not in nucleoprots:
        nucleoprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["NikO-like enolpyruvyl transferase",nikoscores[nikoprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:

          detecteddomainsdict[i] = [["NikO-like enolpyruvyl transferase",nikoscores[nikoprots.index(i)]]]
  #Extract phosphoglycolipid biosynthesis proteins
  phosphoprots = []
  if 1 in geneclustertypes or 17 in geneclustertypes:
    phosphogl = parsehmmoutput(65,hmmoutputfolder + "MoeO5.txt")
    phosphoprots = phosphogl[0]
    phosphoscores = phosphogl[1]
    for i in phosphoprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["MoeO5-like prenyl-3-phosphoglycerate synthase",phosphoscores[phosphoprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["MoeO5-like prenyl-3-phosphoglycerate synthase",phosphoscores[phosphoprots.index(i)]]]
  #Extract melanin biosynthesis proteins
  melaninprots = []
  if 1 in geneclustertypes or 18 in geneclustertypes:
    melanin = parsehmmoutput(40,hmmoutputfolder + "melC.txt")
    melaninprots = melanin[0]
    melaninscores = melanin[1]
    for i in melaninprots:
      if detecteddomainsdict.has_key(i):
        detdomlist = detecteddomainsdict[i]
        detdomlist.append(["MelC-like melanin synthase",melaninscores[melaninprots.index(i)]])
        detecteddomainsdict[i] = detdomlist
      else:
        detecteddomainsdict[i] = [["MelC-like melanin synthase",melaninscores[melaninprots.index(i)]]]
  #Extract other putative secondary metabolite biosynthesis proteins
  otherprots = []
  amp_t_prots = []
  if 1 in geneclustertypes or 19 in geneclustertypes:
    pptb = parsehmmoutput(20,hmmoutputfolder + "PP-binding.txt")
    pptbprots = pptb[0]
    pptbscores = pptb[1]
    cond = parsehmmoutput(20,hmmoutputfolder + "Condensation.txt")
    amp = parsehmmoutput(20,hmmoutputfolder + "AMP-binding.txt")
    ampprots = amp[0]
    ampscores = amp[1]
    ampox = parsehmmoutput(50,hmmoutputfolder + "A-OX.txt")
    ampoxprots = ampox[0]
    ampoxscores = ampox[1]
    nad4 = parsehmmoutput(40,hmmoutputfolder + "NAD_binding_4.txt")
    nad4prots = nad4[0]
    nad4scores = nad4[1]
    cprots = cond[0]
    aprots = amp[0]
    for i in ampox[0]:
      if i not in aprots:
        aprots.append(i)
    nrpsprots2 = []
    for i in cprots:
      if i in aprots:
        nrpsprots2.append(i)
    tprots = pptb[0]
    for i in tprots:
      if i in aprots and i not in nrpsprots2 and i not in aminocoumarinprots:
        otherprots.append(i)
        amp_t_prots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["PP-binding domain",pptbscores[pptbprots.index(i)]])
          if i in ampprots:
            detdomlist.append(["Adenylation domain",ampscores[ampprots.index(i)]])
          elif i in ampoxprots:
            detdomlist.append(["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in ampprots:
            detecteddomainsdict[i] = [["PP-binding domain",pptbscores[pptbprots.index(i)]],["Adenylation domain",ampscores[ampprots.index(i)]]]
          elif i in ampoxprots:
            detecteddomainsdict[i] = [["PP-binding domain",pptbscores[pptbprots.index(i)]],["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]]]
    for i in nad4prots:
      if i in aprots and i not in aminocoumarinprots:
        otherprots.append(i)
        amp_t_prots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["NAD-binding domain 4",nad4scores[nad4prots.index(i)]])
          if i in ampprots:
            detdomlist.append(["Adenylation domain",ampscores[ampprots.index(i)]])
          elif i in ampoxprots:
            detdomlist.append(["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          if i in ampprots:
            detecteddomainsdict[i] = [["NAD-binding domain 4",nad4scores[nad4prots.index(i)]],["Adenylation domain",ampscores[ampprots.index(i)]]]
          elif i in ampoxprots:
            detecteddomainsdict[i] = [["NAD-binding domain 4",nad4scores[nad4prots.index(i)]],["Adenylation domain with integrated oxidase",ampoxscores[ampoxprots.index(i)]]]
    lmbu = parsehmmoutput(50,hmmoutputfolder + "LmbU.txt")
    lmbuprots = lmbu[0]
    lmbuscores = lmbu[1]
    for i in lmbuprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["LmbU-like protein",lmbuscores[lmbuprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["LmbU-like protein",lmbuscores[lmbuprots.index(i)]]]
    goadsporin = parsehmmoutput(500,hmmoutputfolder + "goadsporin_like.txt")
    goadsporinprots = goadsporin[0]
    goadsporinscores = goadsporin[1]
    for i in goadsporinprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Goadsporin-like protein",goadsporinscores[goadsporinprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Goadsporin-like protein",goadsporinscores[goadsporinprots.index(i)]]]
    neocarzinostat = parsehmmoutput(28,hmmoutputfolder + "Neocarzinostat.txt")
    neocarzinostatprots = neocarzinostat[0]
    neocarzinostatscores = neocarzinostat[1]
    for i in neocarzinostatprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Neocarzinostatin-like protein",neocarzinostatscores[neocarzinostatprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Neocarzinostatin-like protein",neocarzinostatscores[neocarzinostatprots.index(i)]]]
    cyanobactin = parsehmmoutput(80,hmmoutputfolder + "cyanobactin_synth.txt")
    cyanobactinprots = cyanobactin[0]
    cyanobactinscores = cyanobactin[1]
    for i in cyanobactinprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Cyanobactin protease",cyanobactinscores[cyanobactinprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Cyanobactin protease",cyanobactinscores[cyanobactinprots.index(i)]]]
    cycdipeptide = parsehmmoutput(110,hmmoutputfolder + "cycdipepsynth.txt")
    cycdipeptideprots = cycdipeptide[0]
    cycdipeptidescores = cycdipeptide[1]
    for i in cycdipeptideprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Cyclodipeptide synthase",cycdipeptidescores[cycdipeptideprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Cyclodipeptide synthase",cycdipeptidescores[cycdipeptideprots.index(i)]]]
    fom1 = parsehmmoutput(750,hmmoutputfolder + "fom1.txt")
    fom1prots = fom1[0]
    fom1scores = fom1[1]
    for i in fom1prots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Fom1-like phosphomutase",fom1scores[fom1prots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Fom1-like phosphomutase",fom1scores[fom1prots.index(i)]]]
    bcpb = parsehmmoutput(400,hmmoutputfolder + "bcpB.txt")
    bcpbprots = bcpb[0]
    bcpbscores = bcpb[1]
    for i in bcpbprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["BcpB-like phosphomutase",bcpbscores[bcpbprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["BcpB-like phosphomutase",bcpbscores[bcpbprots.index(i)]]]
    frbd = parsehmmoutput(350,hmmoutputfolder + "frbD.txt")
    frbdprots = frbd[0]
    frbdscores = frbd[1]
    for i in frbdprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["FrbD-like phosphomutase",frbdscores[frbdprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["FrbD-like phosphomutase",frbdscores[frbdprots.index(i)]]]
    mite = parsehmmoutput(400,hmmoutputfolder + "mitE.txt")
    miteprots = mite[0]
    mitescores = mite[1]
    for i in miteprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["MitE-like CoA-ligase",mitescores[miteprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["MitE-like CoA-ligase",mitescores[miteprots.index(i)]]]
    vlmb = parsehmmoutput(250,hmmoutputfolder + "vlmB.txt")
    vlmbprots = vlmb[0]
    vlmbscores = vlmb[1]
    for i in vlmbprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Valanimycin biosynthesis VlmB domain",vlmbscores[vlmbprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Valanimycin biosynthesis VlmB domain",vlmbscores[vlmbprots.index(i)]]]
    prnb = parsehmmoutput(200,hmmoutputfolder + "prnB.txt")
    prnbprots = prnb[0]
    prnbscores = prnb[1]
    for i in prnbprots:
      if i not in otherprots:
        otherprots.append(i)
        if detecteddomainsdict.has_key(i):
          detdomlist = detecteddomainsdict[i]
          detdomlist.append(["Pyrrolnitrin biosynthesis PrnB domain",prnbscores[prnbprots.index(i)]])
          detecteddomainsdict[i] = detdomlist
        else:
          detecteddomainsdict[i] = [["Pyrrolnitrin biosynthesis PrnB domain",prnbscores[prnbprots.index(i)]]]
  if 5 not in geneclustertypes and 1 not in geneclustertypes:
    nrpsprots = []
  if 4 not in geneclustertypes and 1 not in geneclustertypes:
    t3pksprots = []
  if 3 not in geneclustertypes and 1 not in geneclustertypes:
    t2pksprots = []
  if 2 not in geneclustertypes and 1 not in geneclustertypes:
    t1pksprots = []
    t4pksprots = []
    transatpksprots = []
  #Assemble all core sec met proteins
  allsecmetprots = []
  for i in t1pksprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in transatpksprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in t2pksprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in t3pksprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in t4pksprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in nrpsprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in terpeneprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in lantprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in bcinprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in lactamprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in amglyccyclprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in siderophoreprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in ectprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in butyrprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in indoleprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in nucleoprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in phosphoprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in melaninprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in aminocoumarinprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  for i in otherprots:
    if i not in allsecmetprots:
      allsecmetprots.append(i)
  allsecmetprots.sort()

  if len(allsecmetprots) == 0:
    logfile.write("No secondary metabolite biosynthesis gene clusters detected in this nucleotide file.\n")
    logfile.close()
    print >> sys.stderr, "No secondary metabolite biosynthesis gene clusters detected in this nucleotide file."
    sys.exit(1)

  elapsed = (time.time() - starttime)
  #print "4713Time since start: " + str(elapsed)

  #Extract approximate gene clusters based on hmmsearch results, create list of core PKS / NRPS genes for further analysis (use less strict parameters for this then in gene cluster detection to include all PKS/NRPS domains)
  #Create nucleotide fasta files with sec met gene clusters
  #print "Extracting gene clusters from gbk/embl file using detected signature genes..."
  logfile.write("Extracting gene clusters from gbk/embl file using detected signature genes...\n")
  fastafile = open(genomename + "/clusterblast/geneclusterprots.fasta","w")
  txtfile = open(genomename + "/clusterblast/geneclusters.txt","w")
  wb = Workbook()
  font1 = Font()
  style1 = XFStyle()
  style1.font = font1
  font1.bold = True
  ws0 = wb.add_sheet('0')
  ws0.write(0,0,"Input accession number",style1)
  ws0.write(0,1,"Input name",style1)
  ws0.write(0,2,"Gene cluster type",style1)
  ws0.write(0,3,"Gene cluster genes",style1)
  if clusterblast == "y":
    ws0.write(0,4,"Compound with gene cluster of highest homology",style1)
  protcodes = allsecmetprots
  nuccode = genomename
  gbkfile = open(infile,"r")
  output = gbkfile.read()
  output = output.replace("\r","\n")
  #Extract description of nucleotide from gbk/embl file
  if ".gbk" in infile or ".GBK" in infile or ".gb" in infile or ".GB" in infile or ".genbank" in infile or ".GENBANK" in infile:
    try:
      nucname1 = output.split("ACCESSION   ")[0]
      nucname2 = nucname1.split("DEFINITION  ")[1]
      nucname3 = nucname2.replace("\n","")
      while "  " in nucname3:
        nucname3 = nucname3.replace("  "," ")
      nucname = nucname3
    except(KeyError,IOError,IndexError):
      nucname = "input_nucleotide"
  elif ".embl" in infile or ".EMBL" in infile or ".emb" in infile or ".EMB" in infile:
    try:
      nucname1 = output.split("DE   ")[1]
      nucname2 = nucname1.split("\n")[0]
      nucname3 = nucname2.replace("\n","")
      while "  " in nucname3:
        nucname3 = nucname3.replace("  "," ")
      nucname = nucname3
    except(KeyError,IOError,IndexError):
      nucname = "input_nucleotide"
  protstartlocations = []
  protendlocations = []
  genelist = proteins[2]
  genedict = proteins[3]
  #Save all locations of query proteins on the nucleotide in a list
  for j in protcodes:
    if j in genelist:
      protstart_abs = min([int(genedict[j][0]),int(genedict[j][1])])
      protend_abs = max([int(genedict[j][0]),int(genedict[j][1])])
      protstartlocations.append(protstart_abs)
      protendlocations.append(protend_abs)
  #Identify clusters of genes based on protein locations on the nucleotide
  clusterstarts = []
  clusterends = []
  protstartlocations.sort()
  protendlocations.sort()
  nrlocations = len(protstartlocations)
  a = 0
  for i in protstartlocations:
    if a == 0:
      start = str(i)
      clusterstarts.append(start)
      if len(protendlocations) == 1:
        clusterends.append(protendlocations[a])
    elif a == nrlocations - 1:
      if i < ((protendlocations[a - 1]) + 20000):
        clusterends.append(str(protendlocations[a]))
      else:
        end = str(protendlocations[a - 1])
        clusterends.append(end)
        clusterstarts.append(str(i))
        clusterends.append(str(protendlocations[a]))
    else:
      if i > ((protendlocations[a - 1]) + 20000):
        clusterends.append(str(protendlocations[a - 1]))
        start = str(i)
        clusterstarts.append(start)
      else:
        pass
    a += 1
    lastendlocation = i
  #Extend clusters with 20kb on each side of the identified core genes
  clusterstarts2 = []
  for i in clusterstarts:
    j = int(i) - 20000
    if j < 0:
      j = 0
    clusterstarts2.append(j)
  clusterstarts = clusterstarts2
  clusterends2 = []
  for i in clusterends:
    j = int(i) + 20000
    clusterends2.append(j)
  clusterends = clusterends2
  #For each genbank secondary metabolite gene cluster: extract all proteins and write to fasta, 
  a = 0
  clusterinfo = {}
  geneclusters = []
  geneclustergenes = []
  allcoregenes = []
  for i in clusterstarts:
    cstart = int(i)
    cend = int(clusterends[a])
    a += 1
    clusternr = a
    geneclusters.append(clusternr)
    coregenes = []
    clustergenes = []
    #For each gene in nucleotide, check if it is inside this cluster; if, so append info to list of clustergenes
    if a == 1:
      for i in genelist:
        geneinfo = genedict[i][:-1]
        geneinfo.append(i)
        genedict[i] = geneinfo
    for i in genelist:
      geneinfo = genedict[i]
      genestart = int(geneinfo[0])
      geneend = int(geneinfo[1])
      if (genestart > cstart and genestart < cend) or (geneend > cstart and geneend < cend):
        clustergenes.append(geneinfo)
    #Determine type of cluster
    type = "other"
    z = 0
    for k in clustergenes:
      i = k[4]
      if i in t1pksprots:
        if z == 0:
          type = "t1pks"
        elif "t1pks" not in type:
          type = type + "-t1pks"
        z = 1
      if i in transatpksprots:
        if z == 0:
          type = "transatpks"
        elif "transatpks" not in type:
          type = type + "-transatpks"
        z = 1
      if i in t2pksprots:
        if z == 0:
          type = "t2pks"
        elif "t2pks" not in type:
          type = type + "-t2pks"
        z = 1
      if i in t3pksprots:
        if z == 0:
          type = "t3pks"
        elif "t3pks" not in type:
          type = type + "-t3pks"
        z = 1
      if i in t4pksprots:
        if z == 0:
          type = "t1pks"
        elif "t1pks" not in type:
          type = type + "-t1pks"
        z = 1
      if i in nrpsprots:
        if z == 0:
          type = "nrps"
        elif "nrps" not in type:
          type = type + "-nrps"
        z = 1
      if i in terpeneprots:
        if z == 0:
          type= "terpene"
        elif "terpene" not in type:
          type = type + "-terpene"
        z = 1
      if i in lantprots:
        if z == 0:
          type= "lant"
        elif "lant" not in type:
          type = type + "-lant"
        z = 1
      if i in bcinprots:
        if z == 0:
          type= "bcin"
        elif "bcin" not in type:
          type = type + "-bcin"
        z = 1
      if i in lactamprots:
        if z == 0:
          type = "blactam"
        elif "blactam" not in type:
          type = type + "-blactam"
        z = 1
      if i in amglyccyclprots:
        if z == 0:
          type = "amglyccycl"
        elif "amglyccycl" not in type:
          type = type + "-amglyccycl"
        z = 1
      if i in siderophoreprots:
        if z == 0:
          type = "siderophore"
        elif "siderophore" not in type:
          type = type + "-siderophore"
        z = 1
      if i in ectprots:
        if z == 0:
          type = "ectoine"
        elif "ectoine" not in type:
          type = type + "-ectoine"
        z = 1
      if i in indoleprots:
        if z == 0:
          type = "indole"
        elif "indole" not in type:
          type = type + "-indole"
        z = 1
      if i in nucleoprots:
        if z == 0:
          type = "nucleoside"
        elif "nucleoside" not in type:
          type = type + "-nucleoside"
        z = 1
      if i in phosphoprots:
        if z == 0:
          type = "phosphoglycolipid"
        elif "phosphoglycolipid" not in type:
          type = type + "-phosphoglycolipid"
        z = 1
      if i in butyrprots:
        if z == 0:
          type = "butyrolactone"
        elif "butyrolactone" not in type:
          type = type + "-butyrolactone"
        z = 1
      if i in melaninprots:
        if z == 0:
          type = "melanin"
        elif "melanin" not in type:
          type = type + "-melanin"
        z = 1
      if i in aminocoumarinprots:
        if z == 0:
          type = "aminocoumarin"
        elif "aminocoumarin" not in type:
          type = type + "-aminocoumarin"
        z = 1
    if "other-" in type[:6]:
      type = type[6:]
    #Shorten gene cluster if type is among typically short gene cluster types
    if cend > dnaseqlength:
      cend = dnaseqlength
    if type == "t3pks" or type == "t2pks":
      if cstart != 0:
        cstart = cstart + 5000
      if cend != dnaseqlength:
        cend = cend - 5000
      clustergenes2 = []
      for i in clustergenes:
        start = int(i[0])
        end = int(i[1])
        if (start > cstart and start < cend) or (end > cstart and end < cend):
          clustergenes2.append(i)
      clustergenes = clustergenes2
    if type == "bcin" or type == "siderophore" or type == "lant" or type == "terpene":
      if cstart != 0:
        cstart = cstart + 10000
      if cend != dnaseqlength:
        cend = cend - 10000
      clustergenes2 = []
      for i in clustergenes:
        start = int(i[0])
        end = int(i[1])
        if (start > cstart and start < cend) or (end > cstart and end < cend):
          clustergenes2.append(i)
      clustergenes = clustergenes2
    if type == "butyrolactone" or type == "melanin" or type == "ectoine":
      if cstart != 0:
        cstart = cstart + 17000
      if cend != dnaseqlength:
        cend = cend - 17000
      clustergenes2 = []
      for i in clustergenes:
        start = int(i[0])
        end = int(i[1])
        if (start > cstart and start < cend) or (end > cstart and end < cend):
          clustergenes2.append(i)
      clustergenes = clustergenes2
    #For all clustergenes, write info to fasta
    for i in clustergenes:
      start = str(i[0])
      end = str(i[1])
      strand = i[2]
      seq = seqdict[i[4]]
      ann = i[3].replace(" ","_")
      accession = i[4]
      name = nuccode + "|c" + str(a) + "|" + start + "-" + end + "|" + strand + "|" + accession + "|" + ann
      fastafile.write(">" + name + "\n" + seq + "\n")
      if accession not in geneclustergenes:
        geneclustergenes.append(accession)
    #Write gene cluster info to separate txt file
    txtfile.write(nuccode + "\t" + nucname + "\t" + "c" + str(a) + "\t" + type + "\t")
    ws0.write(a,0,genomic_accnr)
    try:
      ws0.write(a,1,nucname)
    except:
      ws0.write(a,1,"Name to long to be contained in Excel cell; see txt file in downloadable zip archive.")
    ws0.write(a,2,type)
    xlsgenesfield = ""
    for i in clustergenes:
      txtfile.write(i[4] + ";")
      xlsgenesfield = xlsgenesfield + i[4] + ";"
    txtfile.write("\t")
    for i in clustergenes:
      txtfile.write(accessiondict[i[4]] + ";")
    xlsgenesfield = xlsgenesfield[:-1]
    try:
      ws0.write(a,3,xlsgenesfield)
    except:
      ws0.write(a,3,"Too many genes to be contained in Excel cell; see txt file in downloadable zip archive.")
    txtfile.write("\n")
    #Write gene cluster info to clusterinfo dictionary
    for i in clustergenes:
      if i[4] in allsecmetprots:
        coregenes.append(i[4])
        allcoregenes.append(i[4])
    clusterinfo[clusternr] = [type,cstart,cend,coregenes,clustergenes]
  #Close xls, fasta and txt files
  fastafile.close()
  txtfile.close()

  #Analysis of core PKS/NRPS genes (separate py), detect subgroups and predict specificities and final products
  #Make list of PKS / NRPS gene clusters to be analysed
  #print "Analysing core PKS/NRPS genes..."
  logfile.write("Analysing core PKS/NRPS genes...\n")
  pksnrpsgeneclusters = []
  pksnrpscoregenes = []
  for i in geneclusters:
    if "t1pks" in clusterinfo[i][0] or "t4pks" in clusterinfo[i][0] or "transatpks" in clusterinfo[i][0] or "nrps" in clusterinfo[i][0]:
      pksnrpsgeneclusters.append(i)
  for i in t1pksprots:
    pksnrpscoregenes.append(i)
  for i in transatpksprots:
    pksnrpscoregenes.append(i)  
  for i in t4pksprots:
    pksnrpscoregenes.append(i)
  for i in nrpsprots:
    pksnrpscoregenes.append(i)
  for i in amp_t_prots:
    pksnrpscoregenes.append(i)
  pksnrpsgenestartdict = {}
  for i in pksnrpscoregenes:
    start = int(genedict[i][0])
    pksnrpsgenestartdict[i] = start
  pksnrpscoregenes = sortdictkeysbyvalues(pksnrpsgenestartdict)
  nrpsnames = []
  nrpsseqs = []
  pksnrpsnames = []
  pksnrpsseqs = []
  pksnames = []
  pksseqs = []
  calnames = []
  calseqs = []
  krnames = []
  krseqs = []
  nrpspkstypedict = {}
  domaindict = {}
  if len(pksnrpscoregenes) > 0:
    #Write PKS / NRPS core genes to FASTA file
    for i in pksnrpscoregenes:
      name = i
      seq = seqdict[i]
      pksnrpsnames.append(name)
      pksnrpsseqs.append(seq)
    writefasta(pksnrpsnames,pksnrpsseqs,genomename + "/nrpspks_proteins.fasta")
    #Analyse for abMotifs
    hmmsearch = hmmscan_path + " --cpu " + str(nrcpus) + " -E 0.1 -o " + genomename + "/nrpspks/abmotifshmm_output.txt" + " --noali --tblout " + genomename + "/nrpspks/abmotifshmm.txt "+ hmms_path +"abmotifs.hmm " + genomename + "/nrpspks_proteins.fasta"
    os.system(hmmsearch)
    mhmmlengthsdict = hmmlengths(hmms_path+"abmotifs.hmm")
    motifdict = hmmscanparse(genomename + "/nrpspks/abmotifshmm_output.txt",mhmmlengthsdict)
    #Analyse for C/A/PCP/E/KS/AT/ATd/DH/KR/ER/ACP/TE/TD/COM/Docking/MT/CAL domains
    hmmsearch = hmmscan_path + " --cut_tc --cpu " + str(nrcpus) + " -o " + genomename + "/nrpspks/nrpspkshmm_output.txt" + " --noali --tblout " + genomename + "/nrpspks/nrpspkshmm.txt "+ hmms_path +"nrpspksdomains.hmm " + genomename + "/nrpspks_proteins.fasta"
    os.system(hmmsearch)
    hmmlengthsdict = hmmlengths(hmms_path+"nrpspksdomains.hmm")
    domaindict = hmmscanparse(genomename + "/nrpspks/nrpspkshmm_output.txt",hmmlengthsdict)
    nrpspksdomainsfile = open(genomename + "/nrpspks/nrpspksdomains.txt","w")
    #Analyse KS domains & PKS/NRPS protein domain composition to detect NRPS/PKS types
    kshmmsearch = hmmscan_path + " --cut_tc --cpu " + str(nrcpus) + " -o " + genomename + "/nrpspks/kshmm_output.txt" + " --noali --tblout " + genomename + "/nrpspks/kshmm.txt " + hmms_path + "ksdomains.hmm " + genomename + "/nrpspks_proteins.fasta"
    os.system(kshmmsearch)
    kshmmlengthsdict = hmmlengths(hmms_path+"ksdomains.hmm")
    ksdomaindict = hmmscanparse(genomename + "/nrpspks/kshmm_output.txt",kshmmlengthsdict)
    for k in pksnrpscoregenes:
      #structure of domaindict: domaindict[genename] = [[name,start,end,evalue,score],[name,start,end,evalue,score], etc.]
      domainlist = []
      nrKSdomains = 0
      for i in domaindict[k]:
        domainlist.append(i[0])
        if i[0] == "PKS_KS":
          nrKSdomains += 1
      modKSscore = 0
      traKSscore = 0
      eneKSscore = 0
      iterKSscore = 0
      for i in ksdomaindict[k]:
        if i[0] == "Trans-AT-KS":
          traKSscore += 1
        if i[0] == "Modular-KS":
          modKSscore += 1
        if i[0] == "Enediyne-KS":
          eneKSscore += 1
        if i[0] == "Iterative-KS":
          iterKSscore += 1
      for i in domaindict[k]:
        if "Cglyc" in domainlist and "Epimerization" in domainlist and "AMP-binding" in domainlist and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
          type = "Glycopeptide NRPS"
        elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist) and "AMP-binding" in domainlist and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
          type = "NRPS"
        elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist) or "AMP-binding" in domainlist and ("PKS_KS" in domainlist or "PKS_AT" in domainlist):
          type = "Hybrid PKS-NRPS"
        elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" not in domainlist and "Trans-AT_docking" in domainlist and traKSscore > modKSscore and traKSscore > iterKSscore and traKSscore > eneKSscore:
          type = "Type I Trans-AT PKS"
        elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist and iterKSscore > modKSscore and iterKSscore > traKSscore and iterKSscore > eneKSscore and nrKSdomains < 3:
          type = "Type I Iterative PKS"
        elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist and eneKSscore > modKSscore and eneKSscore > traKSscore and eneKSscore > iterKSscore and nrKSdomains < 3:
          type = "Type I Enediyne PKS"
        elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist and ((modKSscore > eneKSscore and modKSscore > traKSscore and modKSscore > iterKSscore) or nrKSdomains > 3):
          type = "Type I Modular PKS"
        elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist:
          type = "PKS-like protein"
        elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist or "AMP-binding" in domainlist) and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
          type = "NRPS-like protein"
        else:
          type = "PKS/NRPS-like protein"
      nrpspkstypedict[k] = type
    #Write data to output file
    for k in pksnrpscoregenes:
      j = domaindict[k]
      l = motifdict[k]
      nrpspksdomainsfile.write(">> " + k + "\n")
      nrpspksdomainsfile.write(">> " + nrpspkstypedict[k] + "\n")
      nrpspksdomainsfile.write("name\tstart\tend\te-value\tscore\n")
      for i in j:
        #nrpspksdomainsfile.write(str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\n")
        nrpspksdomainsfile.write("%s\t%s\t%s\t%s\t%s\n" % (i[0], i[1], i[2], i[3], i[4]) )
      nrpspksdomainsfile.write("** Motifs: **\n")
      for i in l:
          #nrpspksdomainsfile.write(str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\n")
          nrpspksdomainsfile.write("%s\t%s\t%s\t%s\t%s\n" % (i[0], i[1], i[2], i[3], i[4]) )
      nrpspksdomainsfile.write("\n\n")
    nrpspksdomainsfile.close()

    elapsed = (time.time() - starttime)
    #print "5163Time since start: " + str(elapsed)

    #Predict NRPS A domain specificities with NRPSPredictor and Minowa et al. method
    #print "Predicting NRPS A domain substrate specificities by NRPSPredictor"
    logfile.write("Predicting NRPS A domain substrate specificities by NRPSPredictor\n")
    #NRPSPredictor: extract AMP-binding + 120 residues N-terminal of this domain, extract 8 Angstrom residues and insert this into NRPSPredictor
    for k in pksnrpscoregenes:
      j = domaindict[k]
      nr = 0
      for i in j:
        if i[0] == "AMP-binding" or i[0] == "A-OX":
          nr += 1
          start = int(i[1])
          end = int(i[2]) + 120
          seq = seqdict[k][start:end]
          name = k + "_A" + str(nr)
          nrpsnames.append(name)
          nrpsseqs.append(seq)
    if len(nrpsnames) > 0:
      writefasta(nrpsnames,nrpsseqs,"NRPSPredictor2/nrpsseqs.fasta")
      #nrpspredcommand = "perl nrpsSpecPredictor.pl nrpsseqs.fasta ../" + nrpspredictoroutputfolder + " ." #OLD NRPSPREDICTOR1 command
      os.chdir("NRPSPredictor2/")
      #Get NRPSPredictor2 code predictions, output sig file for input for NRPSPredictor2 SVMs
      if sys.platform == ('win32'):
        nrpspred2codecommand = 'nrpscodepred nrpsseqs.fasta input.sig nrpscodes.txt  > nul'
      if sys.platform == ('linux2'):
        nrpspred2codecommand = 'python nrpscodepred.py nrpsseqs.fasta input.sig nrpscodes.txt > /dev/null'
      os.system(nrpspred2codecommand)
      #Run NRPSPredictor2 SVM
      currentdir = os.getcwd()
      if sys.platform == ('win32'):
        nrpspred2command = 'java -Ddatadir="' + currentdir + '\\data" -cp build/NRPSpredictor2.jar;lib/java-getopt-1.0.13.jar;lib/Utilities.jar;lib/libsvm.jar org.roettig.NRPSpredictor2.NRPSpredictor2 -i input.sig -r ..\\' + nrpspredictoroutputfolder + 'nrpspredictor2.out -s 1'
      if sys.platform == ('linux2'):
        nrpspred2command = './NRPSpredictor2.sh -i input.sig -r ../' + nrpspredictoroutputfolder + 'nrpspredictor2.out -s 1'
      os.popen(nrpspred2command)
      #Copy NRPSPredictor results
      if sys.platform == ('win32'):
        copycommand = 'copy/y nrpscodes.txt ..\\' + nrpspredictoroutputfolder.replace("/","\\") + ' > nul'
      if sys.platform == ('linux2'):
        copycommand = 'cp nrpscodes.txt ../' + nrpspredictoroutputfolder + " > /dev/null"
      os.system(copycommand)
      os.chdir("..")
    elapsed = (time.time() - starttime)
    #print "5206Time since start: " + str(elapsed)
    # folgendes bis zum naechsten time braucht 500s, liegt wohl haupsaechlich an schlechtem minowa_A code
    #Minowa method: extract AMP-binding domain, and run Minowa_A
    if len(nrpsnames) > 0:
      #print "Predicting NRPS A domain substrate specificities by Minowa et al. method\n"
      logfile.write("Predicting NRPS A domain substrate specificities by Minowa et al. method")
      nrpsnames2 = []
      nrpsseqs2 = []
      for k in pksnrpscoregenes:
        j = domaindict[k]
        nr = 0
        for i in j:
          if i[0] in ["AMP-binding", "A-OX"]:
            nr += 1
            start = int(i[1])
            end = int(i[2])
            seq = seqdict[k][start:end]
            name = k + "_A" + str(nr)
            nrpsnames2.append(name)
            nrpsseqs2.append(seq)
      writefasta(nrpsnames2,nrpsseqs2,minowanrpsoutputfolder + "nrpsseqs.fasta")
      if sys.platform == ('win32'):
        minowanrpscommand = "minowa_A ../" + minowanrpsoutputfolder + "nrpsseqs.fasta ../" + minowanrpsoutputfolder + "nrpspredoutput.txt"
      if sys.platform == ('linux2'):
        minowanrpscommand = "python minowa_A.py ../" + minowanrpsoutputfolder + "nrpsseqs.fasta ../" + minowanrpsoutputfolder + "nrpspredoutput.txt"
      os.chdir("Minowa/")
      os.system(minowanrpscommand)
      os.chdir("..")

    elapsed = (time.time() - starttime)
    #print "5235Time since start: " + str(elapsed)
    #Predict PKS AT domain specificities with Minowa et al. method and PKS code (NP searcher / ClustScan / own?)
    for k in pksnrpscoregenes:
      j = domaindict[k]
      nr = 0
      for i in j:
        if i[0] == "PKS_AT":
          nr += 1
          start = int(i[1])
          end = int(i[2])
          seq = seqdict[k][start:end]
          name = k + "_AT" + str(nr)
          pksnames.append(name)
          pksseqs.append(seq)  
    if len(pksnames) > 0:
      writefasta(pksnames,pksseqs,pkssignatureoutputfolder + "pksseqs.fasta")
      writefasta(pksnames,pksseqs,minowapksoutputfolder + "pksseqs.fasta")
      #Run PKS signature analysis
      elapsed = (time.time() - starttime)
      #print "5254Time since start: " + str(elapsed)
      print "Predicting PKS AT domain substrate specificities by Yadav et al. PKS signature sequences"
      logfile.write("Predicting PKS AT domain substrate specificities by Yadav et al. PKS signature sequences\n")
      if sys.platform == ('win32'):
        pkspredcommand = "PKS_analysis ../" + pkssignatureoutputfolder + "pksseqs.fasta ../" + pkssignatureoutputfolder + "pkspredoutput.txt"
      if sys.platform == ('linux2'):
        pkspredcommand = "python PKS_analysis.py ../" + pkssignatureoutputfolder + "pksseqs.fasta ../" + pkssignatureoutputfolder + "pkspredoutput.txt"
      os.chdir("pkssignatures/")
      os.system(pkspredcommand)
      os.chdir("..")
      #Minowa method: run Minowa_AT
      elapsed = (time.time() - starttime)
      #print "5266Time since start: " + str(elapsed)
      print "Predicting PKS AT domain substrate specificities by Minowa et al. method"
      logfile.write("Predicting PKS AT domain substrate specificities by Minowa et al. method\n")
      if sys.platform == ('win32'):
        minowapkscommand = "minowa_AT ../" + minowapksoutputfolder + "pksseqs.fasta ../" + minowapksoutputfolder + "pkspredoutput.txt"
      if sys.platform == ('linux2'):
        minowapkscommand = "python minowa_AT.py ../" + minowapksoutputfolder + "pksseqs.fasta ../" + minowapksoutputfolder + "pkspredoutput.txt"
      os.chdir("Minowa/")
      os.system(minowapkscommand)
      os.chdir("..")

    #Predict PKS CAL domain specificities with Minowa et al. method
    elapsed = (time.time() - starttime)
    #print "5279Time since start: " + str(elapsed)
    print "Predicting CAL domain substrate specificities by Minowa et al. method"
    logfile.write("Predicting CAL domain substrate specificities by Minowa et al. method\n")
    for k in pksnrpscoregenes:
      j = domaindict[k]
      nr = 0
      for i in j:
        if i[0] == "CAL_domain":
          nr += 1
          start = int(i[1])
          end = int(i[2])
          seq = seqdict[k][start:end]
          name = k + "_CAL" + str(nr)
          calnames.append(name)
          calseqs.append(seq)  
    if len(calnames) > 0:
      writefasta(calnames,calseqs,minowacaloutputfolder + "calseqs.fasta")
      if sys.platform == ('win32'):
        minowacalcommand = "minowa_CAL ../" + minowacaloutputfolder + "calseqs.fasta ../" + minowacaloutputfolder + "calpredoutput.txt"
      if sys.platform == ('linux2'):
        minowacalcommand = "python minowa_CAL.py ../" + minowacaloutputfolder + "calseqs.fasta ../" + minowacaloutputfolder + "calpredoutput.txt"
      os.chdir("Minowa/")
      os.system(minowacalcommand)
      os.chdir("..")

    elapsed = (time.time() - starttime)
    #print "5305Time since start: " + str(elapsed)
    #Predict PKS KR domain stereochemistry using pattern as published in ClustScan
    print "Predicting PKS KR activity and stereochemistry using KR fingerprints from Starcevic et al."
    logfile.write("Predicting PKS KR activity and stereochemistry using KR fingerprints from Starcevic et al.\n")
    for k in pksnrpscoregenes:
      j = domaindict[k]
      nr = 0
      for i in j:
        if i[0] == "PKS_KR":
          nr += 1
          start = int(i[1])
          end = int(i[2])
          seq = seqdict[k][start:end]
          name = k + "_KR" + str(nr)
          krnames.append(name)
          krseqs.append(seq) 
    if len(krnames) > 0:
      writefasta(krnames,krseqs,kranalysisoutputfolder + "krseqs.fasta")
      if sys.platform == ('win32'):
        kranalysiscommand = "kr_analysis ../" + kranalysisoutputfolder + "krseqs.fasta ../" + kranalysisoutputfolder + "krpredoutput.txt"
      if sys.platform == ('linux2'):
        kranalysiscommand = "python kr_analysis.py ../" + kranalysisoutputfolder + "krseqs.fasta ../" + kranalysisoutputfolder + "krpredoutput.txt"
      os.chdir("kr_analysis/")
      os.system(kranalysiscommand)
      os.chdir("..")

  #Read and parse all substrate specificity prediction output files
  minowa_nrps_preds = {}
  minowa_nrps_preds_details = {}
  nrps_svm_preds = {}
  nrps_svm_preds_details = {}
  nrps_code_preds = {}
  nrps_code_preds_details = {}
  substratetransdict2 = {'pipecolate':'pip','fOHOrn':'orn','beta-Lys':'blys','5NhOrn':'orn','OHOrn':'orn','Aad':'Aaa','bOHTyr':'bht'}
  if len(nrpsnames) > 0:
    minowa_a_file = open(minowanrpsoutputfolder + "nrpspredoutput.txt","r")
    minowa_a_file = minowa_a_file.read()
    minowa_a_file = minowa_a_file.replace("\r","\n")
    parts = minowa_a_file.split("\\\\\n")[1:]
    for i in parts:
      partlines = i.split("\n")
      acc = partlines[0]
      tophit = partlines[2].split("\t")[0]
      if tophit in substratetransdict2.keys():
        tophit = substratetransdict2[tophit]
      minowa_nrps_preds[acc] = tophit.lower()
      minowa_nrps_preds_details[acc] = "<b>Minowa HMM method A-domain<br>Substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"
    nrpspredictorfile1 = open(nrpspredictoroutputfolder + "nrpspredictor2.out","r")
    nrpspredictorfile2 = open(nrpspredictoroutputfolder + "nrpscodes.txt","r")
    nrpspredictorfile1 = nrpspredictorfile1.read()
    nrpspredictorfile1 = nrpspredictorfile1.replace("\r","\n")
    lines = nrpspredictorfile1.split("\n")[1:-1]
    for k in lines:
      tabs = k.split("\t")
      nrps_svm_preds[tabs[0]] = tabs[6]
      nrps_svm_preds_details[tabs[0]] = "<b> NRPSPredictor2 SVM prediction details:</b><br>\n8 Angstrom 34 AA code:<br>\n" + tabs[1] + "<br>\nPredicted physicochemical class:<br>\n" + tabs[3] + "<br>\nLarge clusters prediction:<br>\n" + tabs[4] + "<br>\nSmall clusters prediction:<br>\n" + tabs[5] + "<br>\nSingle AA prediction:<br>\n" + tabs[6] + "<br><br>\n\n"
    nrpspredictorfile2 = nrpspredictorfile2.read()
    nrpspredictorfile2 = nrpspredictorfile2.replace("\r","\n")
    lines = nrpspredictorfile2.split("\n")[:-1]
    for k in lines:
      tabs = k.split("\t")
      nrps_code_preds[tabs[0]] = tabs[1]
      nrps_code_preds_details[tabs[0]] = "<b> NRPSPredictor2 Stachelhaus code prediction:</b><br>\n" + tabs[1] + "<br><br>\n\n"
  minowa_pks_preds_details = {}
  minowa_pks_preds = {}
  pks_code_preds ={}
  pks_code_preds_details ={}
  substratetransdict = {'Malonyl-CoA':'mal','Methylmalonyl-CoA':'mmal','Methoxymalonyl-CoA':'mxmal','Ethylmalonyl-CoA':'emal','Isobutyryl-CoA':'isobut','2-Methylbutyryl-CoA':'2metbut','trans-1,2-CPDA':'trans-1,2-CPDA','Acetyl-CoA':'Acetyl-CoA','Benzoyl-_CoA':'benz','Propionyl-CoA':'prop','3-Methylbutyryl-CoA':'3metbut','Ethylmalonyl-CoA':'Ethyl_mal','CE-Malonyl-CoA':'cemal','2-Rhyd-Malonyl-CoA':'2Rhydmal','CHC-CoA':'CHC-CoA','inactive':'inactive'}
  if len(pksnames) > 0:
    minowa_at_file = open(minowapksoutputfolder + "pkspredoutput.txt","r")
    minowa_at_file = minowa_at_file.read()
    minowa_at_file = minowa_at_file.replace("\r","\n")
    parts = minowa_at_file.split("\\\\\n")[1:]
    for i in parts:
      partlines = i.split("\n")
      acc = partlines[0]
      if substratetransdict.has_key(partlines[2].split("\t")[0]):
        tophit = substratetransdict[partlines[2].split("\t")[0]]
      else:
        tophit = "pk"
      minowa_pks_preds[acc] = tophit
      minowa_pks_preds_details[acc] = "<b>Minowa HMM method AT-domain<br>Substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"  
    pkssignaturefile = open(pkssignatureoutputfolder + "pkspredoutput.txt","r")
    pkssignaturefile = pkssignaturefile.read()
    pkssignaturefile = pkssignaturefile.replace("\r","\n")
    parts = pkssignaturefile.split("//\n")[1:]
    for i in parts:
      partlines = i.split("\n")
      partlines2 = []
      for j in partlines:
        if j != "":
          partlines2.append(j)
      partlines = partlines2
      acc = partlines[0].split("\t")[0]
      if len(partlines) > 2:
        tophit = (partlines[1].split("\t")[0]).split("__")[1]
        pks_code_preds[acc] = tophit
        codes = []
        prots = []
        scores = []
        for i in partlines[1:4]:
          codes.append(i.split("\t")[0])
          prot = i.split("\t")[1]
          prot = prot.replace("_AT"," (AT")
          prot = prot.replace("__","): ")
          prots.append(prot)
          scores.append(i.split("\t")[2])
        if len(prots) >= 3:
          pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[0] + "% identity)<br>\n" + codes[1] + " - " + prots[1] + " : (" + scores[1] + "% identity)<br>\n" + codes[2] + " - " + prots[2] + " : (" + scores[2] + "% identity)<br><br>\n\n"
        elif len(prots) == 2:
          pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[0] + "% identity)<br>\n" + codes[1] + " - " + prots[1] + " : (" + scores[1] + "% identity)<br><br>\n\n"
        elif len(prots) == 1:
          pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[0] + "% identity)<br><br>\n\n"
      else:
        pks_code_preds[acc] = "N/A"
        pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>No AT-domain substrate specificity prediction hits above 40% identity.<br>\n\n"
  minowa_cal_preds = {}
  minowa_cal_preds_details = {}
  if len(calnames) > 0:
    minowa_cal_file = open(minowacaloutputfolder + "calpredoutput.txt","r")
    minowa_cal_file = minowa_cal_file.read()
    minowa_cal_file = minowa_cal_file.replace("\r","\n")
    parts = minowa_cal_file.split("\\\\\n")[1:]
    for i in parts:
      partlines = i.split("\n")
      acc = partlines[0]
      tophit = partlines[2].split("\t")[0]
      minowa_cal_preds[acc] = tophit
      minowa_cal_preds_details[acc] = "<b>Minowa HMM method<br>CAL-domain substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"  
  kr_activity_preds = {}
  kr_stereo_preds = {}
  if len(krnames) > 0:
    krfile = open(kranalysisoutputfolder + "krpredoutput.txt","r")
    krfile = krfile.read()
    krfile = krfile.replace("\r","\n")
    krlines = krfile.split("\n")[:-1]
    for i in krlines:
      tabs = i.split("\t")
      kr_activity_preds[tabs[0]] = tabs[1]
      kr_stereo_preds[tabs[0]] = tabs[2]

  #Combine substrate specificity predictions into consensus prediction
  consensuspreds = {}
  #available_smiles_parts = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','PHE','PRO','SER','THR','TRP','TYR','VAL','MET','ORN','ala','arg','asn','asp','cys','gln','glu','gly','his','ile','leu','lys','phe','pro','ser','thr','trp','tyr','val','met','orn','Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Phe','Pro','Ser','Thr','Trp','Tyr','Val','Met','Orn','MPRO','23DHB','34DHB','2HIVA','PGLY','DAB','BALA','AEO','4MHA','PICO','AAA','DHA','SCY','PIP','BMT','ADDS','mpro','23dhb','34dhb','2hiva','pgly','dab','bala','aeo','4mha','pico','aaa','dha','scy','pip','bmt','adds','Mpro','23Dhb','34Dhb','2Hiva','Pgly','Dab','Bala','Aeo','4Mha','Pico','Aaa','Dha','Scy','Pip','Bmt','Adds','mal','mmal','omal','emal','nrp','pk']
  available_smiles_parts = ['GLY','ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP','SER','THR','ASN','GLN','TYR','CYS','LYS','ARG','HIS','ASP','GLU','MPRO','ORN','PGLY','DAB','BALA','AEO','DHA','PIP','BMT','gly','ala','val','leu','ile','met','pro','phe','trp','ser','thr','asn','gln','tyr','cys','lys','arg','his','asp','glu','aaa','mpro','dhb','2hiva','orn','pgly','dab','bala','aeo','4mha','pico','phg','dha','scy','pip','bmt','adds','aad','abu','hiv','dhpg','bht','3-me-glu','4pPro','ala-b','ala-d','dht','Sal','tcl','lys-b','hpg','hyv-d','iva','vol','mal','mmal','mxmal','emal','nrp','pk','Gly','Ala','Val','Leu','Ile','Met','Pro','Phe','Trp','Ser','Thr','Asn','Gln','Tyr','Cys','Lys','Arg','His','Asp','Glu','Mpro','23Dhb','34Dhb','2Hiva','Orn','Pgly','Dab','Bala','Aeo','4Mha','Pico','Aaa','Dha','Scy','Pip','Bmt','Adds','DHpg','DHB','nrp','pk']
  for i in pksnrpscoregenes:
    nra = 0
    nrat = 0
    nrcal = 0
    j = domaindict[i]
    for k in j:
      if k[0] == "PKS_AT":
        nrat += 1
        preds = []
        preds.append(minowa_pks_preds[i + "_AT" + str(nrat)])
        preds.append(pks_code_preds[i + "_AT" + str(nrat)])
        cpred = "n"
        for l in preds:
          if preds.count(l) > 1:
            if l in available_smiles_parts:
              consensuspreds[i + "_AT" + str(nrat)] = l
            else:
              consensuspreds[i + "_AT" + str(nrat)] = "pk"
            cpred = "y"
        if cpred == "n":
          consensuspreds[i + "_AT" + str(nrat)] = "pk"
      if k[0] == "AMP-binding" or k[0] == "A-OX":
        nra +=1
        preds = []
        preds.append(minowa_nrps_preds[i + "_A" + str(nra)])
        preds.append(nrps_svm_preds[i + "_A" + str(nra)])
        preds.append(nrps_code_preds[i + "_A" + str(nra)])
        cpred = "n"
        for l in preds:
          if preds.count(l) > 1:
            if l in available_smiles_parts:
              consensuspreds[i + "_A" + str(nra)] = l
            else:
              consensuspreds[i + "_A" + str(nra)] = "nrp"
            cpred = "y"
        if cpred == "n":
          consensuspreds[i + "_A" + str(nra)] = "nrp"
      if k[0] == "CAL_domain":
        nrcal += 1
        if minowa_cal_preds[i + "_CAL" + str(nrcal)] in available_smiles_parts:
          consensuspreds[i + "_CAL" + str(nrcal)] = minowa_cal_preds[i + "_CAL" + str(nrcal)]
        else:
          consensuspreds[i + "_CAL" + str(nrcal)] = "pk"

  #Write all prediction details to HTML files for each gene to be used as pop-up window
  domainnamesdict = {}
  for i in pksnrpscoregenes:
    j = domaindict[i]
    domainnames = []
    for k in j:
      domainnames.append(k[0])
    domainnamesdict[i] = domainnames
  for i in pksnrpscoregenes:
    if "PKS_AT" in domainnamesdict[i] or "AMP-binding" in domainnamesdict[i] or "A-OX" in domainnamesdict[i] or "CAL_domain" in domainnamesdict[i]:
      j = domaindict[i]
      nrat = 0
      nra = 0
      nrcal = 0
      nrkr = 0
      for k in j:
        if k[0] == "PKS_AT":
          nrat += 1
          domainname = i + "_AT" + str(nrat)
          htmloutfile = open(substrspecsfolder + domainname + ".html","w")
          htmloutfile.write('<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
          htmloutfile.write(minowa_pks_preds_details[domainname])
          htmloutfile.write(pks_code_preds_details[domainname])
          htmloutfile.write("<b><i>Consensus Predictions: " + consensuspreds[domainname] + "</b></i>")
          htmloutfile.write('\n</body>\n</html>')
          htmloutfile.close()
        if k[0] == "AMP-binding" or k[0] == "A-OX":
          nra += 1
          domainname = i + "_A" + str(nra)
          htmloutfile = open(substrspecsfolder + domainname + ".html","w")
          htmloutfile.write('<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
          htmloutfile.write(nrps_svm_preds_details[domainname])
          htmloutfile.write(nrps_code_preds_details[domainname])
          htmloutfile.write(minowa_nrps_preds_details[domainname])
          htmloutfile.write("<b><i>Consensus Prediction: '" + consensuspreds[domainname] + "'</b></i>")
          htmloutfile.write('\n</body>\n</html>')
          htmloutfile.close()
        if k[0] == "CAL_domain":
          nrcal += 1
          domainname = i + "_CAL" + str(nrcal)
          htmloutfile = open(substrspecsfolder + domainname + ".html","w")
          htmloutfile.write('<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
          htmloutfile.write(minowa_cal_preds_details[domainname])
          htmloutfile.write('\n</body>\n</html>')
          htmloutfile.close()

  elapsed = (time.time() - starttime)
  #print "5541Time since start: " + str(elapsed)
  #Predict biosynthetic gene order in gene cluster using starter domains, thioesterase domains, gene order and docking domains
  compound_pred_dict = {}
  dockingdomainanalysis = []
  nrpspksclusters = []
  a = 1
  for i in geneclusters:
    genecluster = i
    clustercoregenes = clusterinfo[i][3]
    clusterpksnrpsgenes = []
    for j in clustercoregenes:
      if j in pksnrpscoregenes:
        clusterpksnrpsgenes.append(j)
    if len(clusterpksnrpsgenes) > 0:
      nrpspksclusters.append(genecluster)
      pksgenes = 0
      clusterpksgenes = []
      nrpsgenes = 0
      clusternrpsgenes = []
      hybridgenes = 0
      clusterhybridgenes = []
      for j in clusterpksnrpsgenes:
        k = nrpspkstypedict[j]
        if "PKS" in k and "NRPS" not in k:
          pksgenes += 1
          clusterpksgenes.append(j)
        elif "PKS" not in k and "NRPS" in k:
          nrpsgenes += 1
          clusternrpsgenes.append(j)
        elif "PKS/NRPS" in k:
          if ("PKS_KS" in domainnamesdict[j] or "PKS_AT" in domainnamesdict[j]) and ("AMP-binding" not in domainnamesdict[j] and "A-OX" not in domainnamesdict[j] and "Condensation" not in domainnamesdict[j]):
            pksgenes += 1
            clusterpksgenes.append(j)
          elif ("PKS_KS" not in domainnamesdict[j] and  "PKS_AT" not in domainnamesdict[j]) and ("AMP-binding" in domainnamesdict[j] or "A-OX" in domainnamesdict[j] or "Condensation" in domainnamesdict[j]):
            nrpsgenes += 1
            clusternrpsgenes.append(j)
        elif "PKS" in k and "NRPS" in k:
          hybridgenes += 1
          clusterhybridgenes.append(j)
      #If more than three PKS genes, use dock_dom_analysis if possible to identify order
      dock_dom_analysis = "failed"
      if pksgenes > 3 and nrpsgenes == 0 and hybridgenes == 0:
        #print "Predicting PKS gene order by docking domain sequence analysis"
        logfile.write("Predicting PKS gene order by docking domain sequence analysis")
        dockhtmlfile = open(htmlfolder + "docking_analysis" + str(genecluster) + ".html","w")
        #Find first and last genes based on starter module and TE / TD
        startergene = ""
        endinggene = ""
        for k in clusterpksgenes:
          if "Thioesterase" in domainnamesdict[k] or "TD" in domainnamesdict[k]:
            if endinggene == "":
              endinggene = k
            else:
              endinggene = ""
          if len(domainnamesdict[k]) >=2 and  "PKS_AT" == domainnamesdict[k][0] and "ACP" == domainnamesdict[k][1]:
            if startergene == "":
              startergene = k
            else:
              startergene = ""
        if startergene == "":
          for k in clusterpksgenes:
            if len(domainnamesdict[k]) >=3 and "PKS_KS" == domainnamesdict[k][0] and "PKS_AT" == domainnamesdict[k][1] and "ACP" == domainnamesdict[k][2]:
              if startergene == "":
                startergene = k
              else:
                startergene = ""
                break
        #Extract N-terminal 50 residues of each non-starting protein, scan for docking domains using hmmsearch, parse output to locate interacting residues
        ntermintresdict = {}
        ntermnames = []
        ntermseqs = []
        for k in clusterpksgenes:
          if k != startergene:
            ntermnames.append(k)
            seq = seqdict[k]
            ntermseqs.append(seq[:50])
        ntermfasta = "docking_analysis/input.fasta"
        z = 0
        for k in ntermnames:
          writefasta([ntermnames[z]],[ntermseqs[z]],ntermfasta)
          os.chdir("docking_analysis")
          os.system("muscle -profile -quiet -in1 nterm.fasta -in2 input.fasta -out muscle.fasta")
          intresidues = extractpositions("nterm.fasta","muscle.fasta",[2,15],"EryAIII_5_6_ref",ntermnames[z])
          ntermintresdict[ntermnames[z]] = intresidues
          os.chdir("..")
          z += 1
        #Extract C-terminal 100 residues of each non-ending protein, scan for docking domains using hmmsearch, parse output to locate interacting residues
        ctermintresdict = {}
        ctermnames = []
        ctermseqs = []
        for k in clusterpksgenes:
          if k != endinggene:
            ctermnames.append(k)
            seq = seqdict[k]
            ctermseqs.append(seq[-100:])
        ctermfasta = "docking_analysis/input.fasta"
        z = 0
        for k in ctermnames:
          writefasta([ctermnames[z]],[ctermseqs[z]],ctermfasta)
          os.chdir("docking_analysis")
          os.system("muscle -profile -quiet -in1 cterm.fasta -in2 input.fasta -out muscle.fasta")
          intresidues = extractpositions("cterm.fasta","muscle.fasta",[55,64],"EryAII_ref",ctermnames[z])
          ctermintresdict[ctermnames[z]] = intresidues
          os.chdir("..")
          z += 1
        #If docking domains found in all, check for optimal order using interacting residues
        genes_to_order = []
        z = 0
        for k in clusterpksgenes:
          if k == startergene or k == endinggene:
            pass
          else:
            genes_to_order.append(k)
          z += 1
        possible_orders = list(itertools.permutations(genes_to_order,len(genes_to_order)))
        hydrophobic = ["A","V","I","L","F","W","Y","M"]
        positivecharge = ["H","K","R"]
        negativecharge = ["D","E"]
        other = ["C","G","P","S","T","N","Q","X","U"]
        possible_orders_scoredict = {}
        for k in possible_orders:
          score = 0
          interactions = []
          z = 0
          for l in k[:-1]:
            interactions.append([l,k[z + 1]])
            z += 1
          for l in interactions:
            res1a = ctermintresdict[l[0]][0]
            res1b = ntermintresdict[l[1]][0]
            res2a = ctermintresdict[l[0]][1]
            res2b = ntermintresdict[l[1]][1]
            if (res1a in hydrophobic and res1b in hydrophobic) or (res1a in positivecharge and res1b in negativecharge) or (res1a in negativecharge and res1b in positivecharge):
              score += 1
            if (res1a in positivecharge and res1b in positivecharge) or (res1a in negativecharge and res1b in negativecharge):
              score = score - 1
            if (res2a in hydrophobic and res2b in hydrophobic) or (res2a in positivecharge and res2b in negativecharge) or (res2a in negativecharge and res2b in positivecharge):
              score += 1
            if (res2a in positivecharge and res2b in positivecharge) or (res2a in negativecharge and res2b in negativecharge):
              score = score - 1
          possible_orders_scoredict[k] = score
        ranked_orders = sortdictkeysbyvaluesrev(possible_orders_scoredict)
        ranked_orders_part = []
        ranked_orders2 = []
        a = 0
        ranked_orders_len = len(ranked_orders) - 1
        for i in ranked_orders:
          if a == 0:
            score = possible_orders_scoredict[i]
            ranked_orders_part.append(i)
          elif a == ranked_orders_len:
            ranked_orders_part.append(i)
            ranked_orders2 = ranked_orders2 + ranked_orders_part
          else:
            if possible_orders_scoredict[i] == score:
              ranked_orders_part.append(i)
            else:
              ranked_orders_part.reverse()
              ranked_orders2 = ranked_orders2 + ranked_orders_part
              score = possible_orders_scoredict[i]
              ranked_orders_part = []
              ranked_orders_part.append(i)
          a += 1
        ranked_orders = ranked_orders2[:1000]
        geneorders = ranked_orders
        geneorders2 = []
        for l in geneorders:
          geneorder = []
          if startergene != "":
            geneorder.append(startergene)
          [ geneorder.append(m) for m in l ]
          #for m in l:
          #  geneorder.append(m)
          if endinggene != "":
            geneorder.append(endinggene)
          geneorders2.append(geneorder)
        geneorders = geneorders2
        if len(ranked_orders) == 1000:
          dockhtmlfile.write('<html>\n<head>\n<LINK href="style.css" rel="stylesheet" type="text/css">\n</head>\n<body>\nDocking domain analysis.  Score for 1000 highest scoring gene orders:<br><br><table border=1>\n')
        else:
          dockhtmlfile.write('<html>\n<head>\n<LINK href="style.css" rel="stylesheet" type="text/css">\n</head>\n<body>\nDocking domain analysis. Scores for all possible gene orders:<br><br><table border=1>\n')
        dockhtmlfile.write('<tr><td><b>Gene order</b></td><td><b>Score</b></td></tr>\n')
        for l in geneorders:
          string = "<tr><td>"
          for m in l:
            string = string + m + ","
          if startergene != "" and endinggene != "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l[1:-1])])
          elif startergene == "" and endinggene != "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l[:-1])])
          elif startergene != "" and endinggene == "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l[1:])])
          elif startergene == "" and endinggene == "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l)])
          dockhtmlfile.write(string + "</td></tr>\n")
        dockhtmlfile.write('\n</table></body></html>')
        dockhtmlfile.close()
        #print "Predicting PKS gene order by docking domain sequence analysis succeeded."
        #Write html outfile with docking domain analysis output
        #
        logfile.write("Predicting PKS gene order by docking domain sequence analysis succeeded.")
        dockingdomainanalysis.append(genecluster)
      #If NRPS genes, mixed NRPS/PKS genes, PKS genes without detected docking domains, or clusters with a 1-3 PKS genes, assume colinearity
      direction = 0
      for k in clusterpksnrpsgenes:
        if strandsdict[k] == "+":
          direction += 1
        elif strandsdict[k] == "-":
          direction = direction - 1
      if direction < 0:
        clusterpksnrpsgenes.reverse()
      if "Thioesterase" in domainnamesdict[clusterpksnrpsgenes[0]] or "TD" in domainnamesdict[clusterpksnrpsgenes[0]]:
        clusterpksnrpsgenes.reverse()
      geneorder = clusterpksnrpsgenes
      #Generate substrates order from predicted gene order and consensus predictions
      prediction = ""
      for k in geneorder:
        domains = domainnamesdict[k]
        nra = 0
        nrat = 0
        nrcal = 0
        for l in domains:
          if "PKS_AT" in l:
            nrat += 1
            prediction = prediction + consensuspreds[k + "_AT" + str(nrat)] + " "
          if "AMP-binding" in l or "A-OX" in l:
            nra += 1
            prediction = prediction + consensuspreds[k + "_A" + str(nra)] + " "
          if "CAL_domain" in l:
            nrcal += 1
            prediction = prediction + consensuspreds[k + "_CAL" + str(nrcal)] + " "
      prediction = prediction[:-1]
      compound_pred_dict[genecluster] = prediction
    a += 1

  #Combine predictions into a prediction of the final chemical structure and generate images
  os.chdir("NRPeditor")
  failedstructures = []
  for i in geneclusters:
    genecluster = i
    if compound_pred_dict.has_key(genecluster):
      residues = compound_pred_dict[genecluster]
      nrresidues = len(residues.split(" "))
      if nrresidues > 1:
        if sys.platform == ('win32'):
          structcommand = 'main input 100 4000 1000 AA DDV DIM ' + str(nrresidues + 1) + ' "'
        elif sys.platform == ('linux2'):
          structcommand = './main input 100 4000 1000 AA DDV DIM ' + str(nrresidues + 1) + ' "'
        for i in residues.split(" "):
          structcommand = structcommand + i + " "
        structcommand = structcommand + 'TE"'
        smilesinfo = os.popen(structcommand)
        smilesinfo = smilesinfo.read()
        smiles_string = (smilesinfo.split("core peptide: ")[1]).split("\ntermintype")[0]
        if sys.platform == ('linux2'):
          smiles_string.replace("[X]","[*:X]")
          smiles_string2 = ""
          a = 1
          for k in smiles_string:
            if k == "X":
              smiles_string2 = smiles_string2 + str(a)
              a += 1
            else:
              smiles_string2 = smiles_string2 + k
          smiles_string = smiles_string2
        smilesfile = open("genecluster" + str(genecluster) + ".smi","w")
        smilesfile.write(smiles_string)
        smilesfile.close()
        depictstatus = depict_smile(genecluster,structuresfolder)
        if depictstatus == "failed":
          failedstructures.append(genecluster)
    elif clusterinfo[genecluster][0] == "ectoine":
      smiles_string = "CC1=NCCC(N1)C(=O)O"
      smilesfile = open("genecluster" + str(genecluster) + ".smi","w")
      smilesfile.write(smiles_string)
      smilesfile.close()
      depictstatus = depict_smile(genecluster,structuresfolder)
      if depictstatus == "failed":
        failedstructures.append(genecluster)
      elif genecluster in failedstructures:
        del failedstructures[failedstructures.index(genecluster)]
      compound_pred_dict[genecluster] = "ectoine "
  os.chdir("..")

  elapsed = (time.time() - starttime)
  #print "5826 Time since start: " + str(elapsed)
  #ClusterBlast
  if clusterblast == "y":
    #Load gene cluster database into memory
    #print "ClusterBlast: Loading gene clusters database into memory..."
    logfile.write("ClusterBlast: Loading gene clusters database into memory...\n")

    os.chdir(genomename + "/clusterblast")
    #file = open( os.path.join(antismash_path, "clusterblast/geneclusters.txt") ,"r")
    #filetext = file.read()
    #lines = filetext.split("\n")
    clusters = {}
    #for i in open(os.path.join(antismash_path, "clusterblast/geneclusters.txt")):
    bin_path = os.path.join(antismash_path, "clusterblast/geneclusters.bin")
    if os.path.exists( bin_path ):
      clusters = cPickle.load( open(bin_path) )
      #print clusters
    else:
      for line in open( os.path.join(antismash_path, "clusterblast/geneclusters.txt") ,"r"):
        line = line.strip()
        tabs = line.split("\t")
        accession = tabs[0]
        clusterdescription = tabs[1]
        clusternr = tabs[2]
        clustertype = tabs[3]
        clustername = accession + "_" + clusternr
        clustertags = tabs[4].split(";")
        clusterprots = tabs[5].split(";")
        clusters[clustername] = [clusterprots,clusterdescription,clustertype,clustertags]
      cPickle.dump(clusters, open(bin_path, 'w'), -1)
    #Load gene cluster database proteins info into memory
    #print "ClusterBlast: Loading gene cluster database proteins into memory..."
    logfile.write("ClusterBlast: Loading gene cluster database proteins into memory...\n")
    #file = open( os.path.join(antismash_path, "clusterblast/geneclusterprots.fasta") ,"r")
    #filetext = file.read()
    #filetext = filetext.replace("\r","\n")
    #lines = filetext.split("\n")
    proteingeneclusters = {}
    proteinlocations = {}
    proteinstrands = {}
    proteinannotations = {}
    proteintags = {}
    bin_path = os.path.join(antismash_path, "clusterblast/geneclusterprots.fasta.bin")
    if os.path.exists( bin_path ):
      (proteingeneclusters, proteinlocations, proteinstrands, proteinannotations, proteintags) = cPickle.load( open(bin_path, 'r') )
    else:
      for line in open( os.path.join(antismash_path, "clusterblast/geneclusterprots.fasta") ,"r"):
        line = line.replace('\n', '')
        if line.startswith(">"):
          tabs = line.split("|")
          #print 'Protein:', tabs
          protein = tabs[6]
          locustag = tabs[4]
          if accessiondict.has_key(locustag):
            locustag = "h_" + locustag
          proteintags[protein] = locustag
          clustername = tabs[0] + "_" + tabs[1]
          proteingeneclusters[protein] = clustername
          location = tabs[2]
          proteinlocations[protein] = location
          strand = tabs[3]
          proteinstrands[protein] = strand
          annotation = tabs[5]
          proteinannotations[protein] = annotation
      cPickle.dump([proteingeneclusters, proteinlocations, proteinstrands, proteinannotations, proteintags], open(bin_path, 'w'), -1)
    #Run BLAST on gene cluster proteins of each cluster on itself to find internal homologs, store groups of homologs - including singles - in a dictionary as a list of lists accordingly
    #print "Finding internal homologs in each gene cluster.."
    logfile.write("Finding internal homologs in each gene cluster..\n")
    internalhomologygroupsdict = {}
    for i in geneclusters:
      clusternumber = i
      #Create input fasta files for BLAST search
      queryclusterprotslist = clusterinfo[i][4]
      queryclusterprots = []
      for i in queryclusterprotslist:
        queryclusterprots.append(i[4])
      queryclusternames = []
      queryclusterseqs = []
      for i in queryclusterprots:
        seq = seqdict[i]
        name = fullnamedict[i]
        queryclusterseqs.append(seq)
        queryclusternames.append(name)
      writefasta(queryclusternames,queryclusterseqs,"internal_input.fasta")
      #Run and parse BLAST search
      makeblastdbcommand = "makeblastdb -in internal_input.fasta -out internal_input.fasta -dbtype prot"
      blastsearch = "blastp  -db internal_input.fasta -query internal_input.fasta -outfmt 6 -max_target_seqs 1000 -evalue 1e-05 -out internal_input.out"
      if "--gui" in sys.argv and sys.argv[sys.argv.index("--gui") + 1] == "y":
        os.popen(makeblastdbcommand)
        os.popen(blastsearch)
      else:
        os.system(makeblastdbcommand)
        os.system(blastsearch)
      #print "5920 makeblastdb finised"
      blastoutput = open("internal_input.out","r").read()
      minseqcoverage = 25
      minpercidentity = 30
      seqlengths = fastaseqlengths(proteins)
      iblastinfo = blastparse(blastoutput,minseqcoverage,minpercidentity,seqlengths,geneclustergenes)
      iblastdict = iblastinfo[0]
      iquerylist = iblastinfo[1]
      #find and store internal homologs
      groups = []
      for j in queryclusternames:
        jsplit = j.split("|")[4]
        if iblastdict.has_key(j):
          hits = iblastdict[j][0]
          group = []
          for k in hits:
            if k[:2] == "h_":
              group.append(k[2:])
            elif k.count("|") > 4:
              group.append(k.split("|")[4])
            else:
              group.append(k)
          if jsplit not in group:
            group.append( jsplit )
          x = 0
          for l in groups:
            for m in group:
              if m in l:
                del groups[x]
                [group.append(n) for n in l if n not in group]
                #for n in l:
                #  if n not in group:
                #    group.append(n)
                break
            x += 1
          group.sort()
          groups.append(group)
        else:
          groups.append([ jsplit ])
      internalhomologygroupsdict[clusternumber] = groups

    #Run BLAST on gene cluster proteins of each cluster and parse output
    #print "5961 Running NCBI BLAST+ gene cluster searches.."
    logfile.write("Running NCBI BLAST+ gene cluster searches..\n")
    for i in geneclusters:
      clusternumber = i
      #print "   Gene cluster " + str(clusternumber)
      #Create input fasta files for BLAST search
      queryclusterprotslist = clusterinfo[i][4]
      queryclusterprots = []
      for i in queryclusterprotslist:
        queryclusterprots.append(i[4])
      queryclusternames = []
      queryclusterseqs = []
      for i in queryclusterprots:
        seq = seqdict[i]
        name = fullnamedict[i]
        queryclusterseqs.append(seq)
        queryclusternames.append(name)
      equalpartsizes = int(len(queryclusternames)/nrcpus)
      for i in range(nrcpus):
        if i == 0:
          setnames = queryclusternames[:equalpartsizes]
          setseqs = queryclusterseqs[:equalpartsizes]
        elif i == (nrcpus - 1):
          setnames = queryclusternames[(i*equalpartsizes):]
          setseqs = queryclusterseqs[(i*equalpartsizes):]
        else:
          setnames = queryclusternames[(i*equalpartsizes):((i+1)*equalpartsizes)]
          setseqs = queryclusterseqs[(i*equalpartsizes):((i+1)*equalpartsizes)]
        writefasta(setnames,setseqs,"input" + str(i) + ".fasta")
      processes = []
      processnames = []
      for i in range(nrcpus):
        processes.append(Process(target=runblast, args=["input" + str(i) + ".fasta"]))
      [i.start() for i in processes]
      time.sleep(10)
      while True:
        processrunning = "n"
        for i in processes:
          if i.is_alive():
            processrunning = "y"
        if processrunning == "y":
          time.sleep(5)
        else:
          break
      [i.join() for i in processes]
      blastoutput = ""
      for i in range(nrcpus):
        output = open("input" + str(i) + ".out","r")
        output = output.read()
        blastoutput = blastoutput + output
      os.chdir("..")
      blastoutputfile = open("./clusterblastoutput.txt","w")
      blastoutputfile.write(blastoutput)
      blastoutputfile.close()
      os.chdir("clusterblast")
      #print "   Blast search finished. Parsing results..."
      logfile.write("   Blast search finished. Parsing results...\n")
      minseqcoverage = 25
      minpercidentity = 30
      seqlengths = fastaseqlengths(proteins)
      blastinfo = blastparse(blastoutput,minseqcoverage,minpercidentity,seqlengths,geneclustergenes)
      blastdict = blastinfo[0]
      querylist = blastinfo[1]
      #Remove queries without hits
      querylist2 = []
      for i in querylist:
        if blastdict.has_key(i):
          querylist2.append(i)
        else:
          pass
      querylist = querylist2
      hitclusters = blastinfo[2]
      #Score BLAST output on all gene clusters
      #Rank gene cluster hits based on 1) number of protein hits covering >25% sequence length or at least 100aa alignment, with >30% identity and 2) cumulative blast score
      #Find number of protein hits and cumulative blast score for each gene cluster
      #print "   Scoring Blast outputs on database of gene clusters..."
      logfile.write("   Scoring Blast outputs on database of gene clusters...\n")
      hitclusterdict = {}
      hitclusterdata = {}
      for i in hitclusters:
        hitclusterdatalist = []
        nrhits = float(0)
        nrcoregenehits = float(0)
        cumblastscore = float(0)
        hitpositions = []
        hitposcorelist = []
        for j in querylist:
          querynrhits = 0
          querycumblastscore = float(0)
          nrhitsplus = "n"
          for k in blastdict[j][0]:
            if i == blastdict[j][1][k][0]:
              if [querylist.index(j),clusters[i][0].index(blastdict[j][1][k][9])] not in hitpositions:
                nrhitsplus = "y"
                querynrhits += 1
                blastscore = float(blastdict[j][1][k][6]) / 1000000
                querycumblastscore = querycumblastscore + blastscore
                hitclusterdatalist.append([j,k,blastdict[j][1][k][5],blastdict[j][1][k][6],blastdict[j][1][k][7],blastdict[j][1][k][8]])
                hitclusterdata[i] = hitclusterdatalist
                hitpositions.append([querylist.index(j),clusters[i][0].index(blastdict[j][1][k][9])])
          if nrhitsplus == "y":
            nrhits += 1
            if j.split("|")[4] in allcoregenes:
              nrcoregenehits += 0.1
              for hit in range(querynrhits):
                hitposcorelist.append(1)
            else:
              for hit in range(querynrhits):
                hitposcorelist.append(0)
            cumblastscore = cumblastscore + float(querycumblastscore)
        query_givenscores_querydict = {}
        query_givenscores_hitdict = {}
        #Find groups of hits
        hitgroupsdict = {}
        for p in hitpositions:
          if not hitgroupsdict.has_key(p[0]):
            hitgroupsdict[p[0]] = [p[1]]
          else:
            hitgroupsdict[p[0]].append(p[1])
        #Calculate synteny score; give score only if more than one hits (otherwise no synteny possible), and only once for every query gene and every hit gene
        synteny_score = 0
        z = 1
        if nrhits > 1:
          for p in hitpositions[:-1]:
            tandem = "n"
            #Check if a gene homologous to this gene has already been scored for synteny in the previous entry
            if p[1] in hitgroupsdict[hitpositions[z][0]]:
              tandem = "y"
            #Score entry
            if ((not query_givenscores_querydict.has_key(p[0])) or query_givenscores_querydict[p[0]] == 0) and ((not query_givenscores_hitdict.has_key(p[1])) or query_givenscores_hitdict[p[1]] == 0) and tandem == "n":
              q = hitpositions[z]
              if (abs(p[0] - q[0]) < 2) and abs(p[0]-q[0]) == abs(p[1]-q[1]):
                synteny_score += 1
                if hitposcorelist[z - 1] == 1 or hitposcorelist[z] == 1:
                  synteny_score += 1
                query_givenscores_querydict[p[0]] = 1
                query_givenscores_hitdict[p[1]] = 1
              else:
                query_givenscores_querydict[p[0]] = 0
                query_givenscores_hitdict[p[1]] = 0
            z += 1
        #Give bonus to gene clusters with >0 core gene hits
        if nrcoregenehits > 0:
          corebonus = 3
        else:
          corebonus = 0
        #sorting score is based on number of hits (discrete values) & cumulative blast score (behind comma values)
        sortingscore = nrhits + synteny_score + corebonus + nrcoregenehits + cumblastscore
        hitclusterdict[i] = sortingscore
      #Sort gene clusters
      rankedclusters = sortdictkeysbyvaluesrev(hitclusterdict)
      rankedclustervalues = sortdictkeysbyvaluesrevv(hitclusterdict)
      #Output for each hit: table of genes and locations of input cluster, table of genes and locations of hit cluster, table of hits between the clusters
      #print "   Writing output file..."
      logfile.write("   Writing output file...\n")
      #os.chdir("..")
      #os.chdir(genomename)
      #os.chdir("clusterblast")
      out_file = open("cluster" + str(clusternumber) + ".txt","w")
      out_file.write("ClusterBlast scores for " + infile)
      out_file.write("\n\nTable of genes, locations, strands and annotations of query cluster:\n")
      #out_file.write("\n")
      #out_file.write("Table of genes, locations, strands and annotations of query cluster:")
      #out_file.write("\n")
      for i in queryclusterprots:
        out_file.write("%s\t%s\t%s\t%s\t%s\t\n" % (i, proteins[3][i][0], proteins[3][i][1], proteins[3][i][2], proteins[3][i][3]))
        """out_file.write(i)
        out_file.write("\t")
        out_file.write(proteins[3][i][0])
        out_file.write("\t")
        out_file.write(proteins[3][i][1])
        out_file.write("\t")
        out_file.write(proteins[3][i][2])
        out_file.write("\t")
        out_file.write(proteins[3][i][3])
        out_file.write("\t")
        out_file.write("\n")"""
      out_file.write("\n\nSignificant hits: \n")
      #out_file.write("\n")
      #out_file.write("Significant hits: ")
      #out_file.write("\n")
      z = 0
      for i in rankedclusters[:100]:
        #out_file.write(str(z+1) + ". " + i + "\t" + clusters[i][1])
        #out_file.write("\n")
        out_file.write("%s. %s\t%s\n" % ((z+1), i, clusters[i][1]) )
        z += 1
      out_file.write("\n\n")
      #out_file.write("\n")
      z = 0
      out_file.write("Details:")
      for i in rankedclusters[:100]:
        value = str(rankedclustervalues[z])
        nrhits = value.split(".",1)[0]
        if nrhits > 0:
          cumblastscore = str(int(float(value.split(".")[1])))
          out_file.write("\n\n>>\n\n%s. %s\nSource: %s\nType: %s\nNumber of proteins with BLAST hits to this cluster: %s\nCumulative BLAST score: %s\n\nTable of genes, locations, strands and annotations of subject cluster:\n" % (z+1, i, clusters[i][1], clusters[i][2], nrhits, cumblastscore))
          clusterproteins = clusters[i][0]
          #print 'clusterproteins\n\n', clusterproteins
          """out_file.write("\n\n")
          out_file.write(">>")
          out_file.write("\n")
          cumblastscore = str(int(float(value.split(".")[1])))
          out_file.write("\n")
          out_file.write(str(z+1) + ". " + i)
          out_file.write("\n")
          out_file.write("Source: " + clusters[i][1])
          out_file.write("\n")
          out_file.write("Type: " + clusters[i][2])
          out_file.write("\n")
          out_file.write("Number of proteins with BLAST hits to this cluster: " + nrhits)
          out_file.write("\n")
          out_file.write("Cumulative BLAST score: " + cumblastscore)
          out_file.write("\n")
          out_file.write("\n")
          out_file.write("Table of genes, locations, strands and annotations of subject cluster:")
          out_file.write("\n")
          clusterproteins = clusters[i][0]"""

          for j in clusterproteins:
            #print '##########asdfasdf######', j, '---'+proteinlocations.keys()[0]+ '---', proteinannotations.has_key(j), proteinstrands.has_key(j), proteinlocations.has_key(j)
            if proteinlocations.has_key(j) and proteinannotations.has_key(j) and proteinstrands.has_key(j):
              if proteintags[j] == "no_locus_tag":
                out_file.write(j)
              else:
                out_file.write(proteintags[j])
              out_file.write( "\t%s\t%s\t%s\t%s\t%s\n" % (j, proteinlocations[j].split("-")[0], proteinlocations[j].split("-")[1], proteinstrands[j], proteinannotations[j]) )
              """out_file.write("\t")
              out_file.write(j)
              out_file.write("\t")
              out_file.write(proteinlocations[j].split("-")[0])
              out_file.write("\t")
              out_file.write(proteinlocations[j].split("-")[1])
              out_file.write("\t")
              out_file.write(proteinstrands[j])
              out_file.write("\t")
              out_file.write(proteinannotations[j])
              out_file.write("\n")
              """

          out_file.write("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
          if i in hitclusterdata.keys():
            tabledata = hitclusterdata[i]
            for x in tabledata:
              w = 0
              for y in x:
                if w == 0:
                  out_file.write( "%s\t" % y.split("|")[4] )
                  #out_file.write("\t")
                  w += 1
                else:
                  out_file.write("%s\t" % y)
                  #out_file.write("\t")
              out_file.write("\n") 
          else:
            "data not found"
            out_file.write("\n") 
          out_file.write("\n")
          z += 1
      #os.chdir("..")
      #os.chdir("..")
      #os.chdir("clusterblast")
    os.chdir("..")
    out_file.close()

  elapsed = (time.time() - starttime)
  #print "Time since start: " + str(elapsed)
  #smCOG analysis
  smcogtreedict = {}
  if smcogs == "y":
    #print "Performing smCOG analysis"
    logfile.write("Performing smCOG analysis\n")
    hmmsearch = hmmscan_path + " --cpu " + str(nrcpus) + " -E 1E-6 -o " + "./smcogs/smcogshmm_output.txt" + " --noali --tblout " + "./smcogs/smcogshmm.txt "+ hmms_path +"smcogs.hmm " + "./clusterblast/geneclusterprots.fasta"
    #print hmmsearch
    os.system(hmmsearch)
    #print 'finised'
    smcoghmmlengthsdict = hmmlengths(hmms_path+"smcogs.hmm")
    smcogdict = hmmscanparse("./smcogs/smcogshmm_output.txt", smcoghmmlengthsdict)
    smcogdict2 = {}
    for i in smcogdict.keys():
      newkey = i.split("|")[4]
      smcogdict2[newkey] = smcogdict[i]
    smcogdict = smcogdict2
    #Write output
    #os.chdir(genomename)
    os.chdir("smcogs")
    smcogfile = open("smcogs.txt","w")
    for k in geneclustergenes:
      if k not in pksnrpscoregenes:
        l = smcogdict[k]
        smcogfile.write(">> " + k + "\n")
        smcogfile.write("name\tstart\tend\te-value\tscore\n")
        smcogfile.write("** smCOG hits **\n")
        for i in l:
            smcogfile.write(str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\n")
        smcogfile.write("\n\n")
    smcogfile.close()
    os.chdir("..")
    os.chdir("..")
    #smCOG phylogenetic tree construction
    #print "Calculating and drawing phylogenetic trees of cluster genes with smCOG members"
    logfile.write("Calculating and drawing phylogenetic trees of cluster genes with smCOG members")
    os.chdir("smcogtree")
    smcoganalysisgenes = []
    #for k in geneclustergenes:
    #  if k not in pksnrpscoregenes:
    #    smcoganalysisgenes.append(k)
    [smcoganalysisgenes.append(k) for k in geneclustergenes if k not in pksnrpscoregenes]
    smcogsets = []
    equalpartsizes = int(len(smcoganalysisgenes)/nrcpus)
    for i in range(nrcpus):
      if i == 0:
        geneslist = smcoganalysisgenes[:equalpartsizes]
      elif i == (nrcpus - 1):
        geneslist = smcoganalysisgenes[(i*equalpartsizes):]        
      else:
        geneslist = smcoganalysisgenes[(i*equalpartsizes):((i+1)*equalpartsizes)]        
      smcogsets.append(geneslist)
    processes = []
    processnames = []
    z = 0
    for k in smcogsets:
      processes.append(Process(target=smcog_analysis, args=[k,z,accessiondict,seqdict,smcogdict,smcogsoutputfolder]))
      z += 1
    for k in processes:
      k.start()
    time.sleep(1)
    while True:
      processrunning = "n"
      for k in processes:
        if k.is_alive():
          processrunning = "y"
      if processrunning == "y":
        time.sleep(5)
      else:
        break
    for k in processes:
      k.join()
    os.chdir("..")
    currentpath = os.getcwd()
    os.chdir(smcogsoutputfolder)
    dircontents = getdircontents()
    for k in dircontents:
      #POTENTIAL pERFORMANCE gainfor k in glob.glob('*.png'):
      if ".png" in k:
        tag = k.split(".png")[0]
        smcogtreedict[tag] = tag + ".png"
    os.chdir(currentpath)


  ##Visualization
  #Read in ClusterBlast data
  #Read in PubMed / PubChem links of database gene clusters
  if clusterblast == "y":
    if genomename in os.getcwd():
        os.chdir('..')
    pubmed_dict = {}
    pubchem_dict = {}
    known_compound_dict = {}
    #pubfile = open(antismash_path + "pubmed_pubchem_links.txt","r")
    #pubfile = pubfile.read()
    #publines = pubfile.split("\n")
    #for i in publines:
    bin_path = os.path.join(antismash_path, "pubmed_pubchem_links.bin")
    if os.path.exists( bin_path ):
      (pubmed_dict, pubchem_dict, known_compound_dict) = cPickle.load( open(bin_path) )
    else:
      for line in open(antismash_path + "pubmed_pubchem_links.txt","r"):
        line = line.replace('\n', '')
        tabs = line.split("\t")
        acc = tabs[0]
        if tabs[1] != "":
          pubmed_dict[acc] = tabs[1]
        if tabs[2] != "":
          pubchem_dict[acc] = tabs[2]
        if tabs[3] != "":
          known_compound_dict[acc] = tabs[3]
        cPickle.dump([pubmed_dict, pubchem_dict, known_compound_dict], open(bin_path, 'w'), -1)
  #print "Writing visualization SVGs and XHTML"
  logfile.write("Writing visualization SVGs and XHTML\n")
  queryclusterdata = {}
  nrhitgeneclusters = {}
  cblastclusternr = 1
  #print os.getcwd()
  if clusterblast == "y":
    for x in geneclusters:
      clusterblastfile = open(clusterblastoutputfolder + "cluster" + str(x) + ".txt","r")
      #print clusterblastfile
      clusterblastfile = clusterblastfile.read()
      clusterblastfile = clusterblastfile.replace("\r","\n")
      toptenhitclusters = []
      #Identify top ten hits for visualization
      hitlines = ((clusterblastfile.split("Significant hits: \n")[1]).split("\nDetails:")[0]).split("\n")
      #print '\n\n#######hitlines\n', hitlines
      a = 0
      cb_accessiondict = {}
      b = 1
      for i in hitlines:
        if " " in i:
          cb_accessiondict[b] = (i.split("\t")[0]).split(" ")[1]
        if genomic_accnr == "" or genomic_accnr not in i:
          b += 1
          if a < 10:
            if len(i) < 80:
              toptenhitclusters.append(i)  
            elif len(i) >= 80:
              j = i[0:77] + "..."
              toptenhitclusters.append(j)
          a += 1
      #print clusterblastfile
      details = (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:]
      #print details
      nrhitclusters = len(toptenhitclusters)
      #Save query gene cluster data
      querylines = ((clusterblastfile.split("Table of genes, locations, strands and annotations of query cluster:\n")[1]).split("\n\n\nSignificant hits:")[0]).split("\n")
      queryclustergenes = []
      queryclustergenesdetails = {}
      for i in querylines:
        tabs = i.split("\t")
        queryclustergenes.append(tabs[0])
        queryclustergenesdetails[tabs[0]] = [tabs[1],tabs[2],tabs[3],tabs[4]]
      #For every gene cluster, store hit genes and details
      colorgroupsdict = {}
      hitclusterdata = {}
      hitclusternr = 1
      compound_found = "n"
      nrhitgeneclusters[x] = 0
      for i in details:
        hitclustergenes = []
        hitclustergenesdetails = {}
        #Only calculate for first ten hit gene clusters
        if genomic_accnr == "" or genomic_accnr not in i:
          if hitclusternr <= 10:
            nrhitgeneclusters[x] = hitclusternr
            accession = cb_accessiondict[hitclusternr]
            hitclustergeneslines = ((i.split("Table of genes, locations, strands and annotations of subject cluster:\n")[1]).split("\n\nTable of Blast hits ")[0]).split("\n")
            #print '***********\n', i, '\n'
            #print hitclustergeneslines
            for j in hitclustergeneslines:
              tabs = j.split("\t")
              hitclustergenes.append(tabs[0])
              hitclustergenesdetails[tabs[0]] = [tabs[2],tabs[3],tabs[4],tabs[5],tabs[1]]

            blasthitslines = ((i.split("%coverage, e-value):\n")[1]).split("\n\n")[0]).split("\n")
            querygeneswithhits = []
            coregeneswithhits = []
            
            
            blasthitdict = {}
            blastdetailsdict = {}
            querygenes = []
            revblasthitdict = {}
            hitgenes = []
            
            
            for k in blasthitslines:
              tabs = k.split("\t")
              if tabs[0] not in querygeneswithhits:
                querygeneswithhits.append(tabs[0])
              if tabs[0] in allcoregenes and tabs[0] not in coregeneswithhits:
                coregeneswithhits.append(tabs[0])


              if blasthitdict.has_key(tabs[0]):
                hits = blasthitdict[tabs[0]]
                hits.append(tabs[1])
                blasthitdict[tabs[0]] = hits
                if revblasthitdict.has_key(tabs[1]):
                  revhits = revblasthitdict[tabs[1]]
                  revhits.append(tabs[0])
                  revblasthitdict[tabs[1]] = revhits
                else:
                  revblasthitdict[tabs[1]] = [tabs[0]]
                blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[5],tabs[3]]
                if tabs[0] not in querygenes:
                  querygenes.append(tabs[0])
                hitgenes.append(tabs[1])
              else:
                blasthitdict[tabs[0]] = [tabs[1]]
                if revblasthitdict.has_key(tabs[1]):
                  revhits = revblasthitdict[tabs[1]]
                  revhits.append(tabs[0])
                  revblasthitdict[tabs[1]] = revhits
                else:
                  revblasthitdict[tabs[1]] = [tabs[0]]
                blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[5],tabs[3]]
                if tabs[0] not in querygenes:
                  querygenes.append(tabs[0])
                hitgenes.append(tabs[1])



            for k in known_compound_dict.keys():
              if k in i and compound_found == "n" and len(querygeneswithhits) > 2 and len(coregeneswithhits) > 0:
                ws0.write(x,4,known_compound_dict[k])
                compound_found = "y"
            """blasthitdict = {}
            blastdetailsdict = {}
            querygenes = []
            revblasthitdict = {}
            hitgenes = []
            for i in blasthitslines:
              tabs = i.split("\t")
              if blasthitdict.has_key(tabs[0]):
                hits = blasthitdict[tabs[0]]
                hits.append(tabs[1])
                blasthitdict[tabs[0]] = hits
                if revblasthitdict.has_key(tabs[1]):
                  revhits = revblasthitdict[tabs[1]]
                  revhits.append(tabs[0])
                  revblasthitdict[tabs[1]] = revhits
                else:
                  revblasthitdict[tabs[1]] = [tabs[0]]
                blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[5],tabs[3]]
                if tabs[0] not in querygenes:
                  querygenes.append(tabs[0])
                hitgenes.append(tabs[1])
              else:
                blasthitdict[tabs[0]] = [tabs[1]]
                if revblasthitdict.has_key(tabs[1]):
                  revhits = revblasthitdict[tabs[1]]
                  revhits.append(tabs[0])
                  revblasthitdict[tabs[1]] = revhits
                else:
                  revblasthitdict[tabs[1]] = [tabs[0]]
                blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[5],tabs[3]]
                if tabs[0] not in querygenes:
                  querygenes.append(tabs[0])
                hitgenes.append(tabs[1])
            """
            #Make groups of genes for coloring
            colorgroups = []
            internalgroups = internalhomologygroupsdict[x]
            for i in internalgroups:
              querygenes_and_hits = []
              for j in i:
                #Make list of query gene and its hits
                additionalhits = []
                #For each hit, check if it was also hit by another gene; if so, only add it to the group if this hit had the lowest blast score
                otherscores = []
                queryscore = 0
                if blasthitdict.has_key(j):
                  for k in blasthitdict[j]:
                    for l in blastdetailsdict.keys():
                      if k in l and j in l:
                        queryscore = blastdetailsdict[l][1]
                      elif k in l and j not in l:
                        otherscores.append(blastdetailsdict[l][1])
                    allscores = otherscores + [queryscore]
                    if queryscore == max(allscores):
                      additionalhits.append(k)
                  #Add additional hits to the querygenes_and_hits list that will form a colorgroup
                  querygenes_and_hits = querygenes_and_hits + additionalhits
                  if j not in querygenes_and_hits:
                    querygenes_and_hits.append(j)
              if len(querygenes_and_hits) > 0:
                colorgroups.append(querygenes_and_hits)
            colorgroupsdict[hitclusternr] = colorgroups
            hitclusterdata[hitclusternr] = [colorgroupsdict,hitclustergenes,hitclustergenesdetails,queryclustergenes,queryclustergenesdetails,toptenhitclusters,accession]
            hitclusternr += 1
        elif hitclusternr > 10 and hitclusternr <= 50:
          blasthitslines = ((i.split("%coverage, e-value):\n")[1]).split("\n\n")[0]).split("\n")
          querygeneswithhits = []
          coregeneswithhits = []
          for k in blasthitslines:
            tabs = k.split("\t")
            if tabs[0] not in querygeneswithhits:
              querygeneswithhits.append( tabs[0] )
            if tabs[0] in allcoregenes and tabs[0] not in coregeneswithhits:
              coregeneswithhits.append(tabs[0])
          for k in known_compound_dict.keys():
            if k in i and compound_found == "n" and len(querygeneswithhits) > 2 and len(coregeneswithhits) > 0:
              ws0.write(x,4,known_compound_dict[k])
              compound_found = "y"
          hitclusternr += 1
      queryclusterdata[cblastclusternr] = [nrhitclusters,hitclusterdata]
      cblastclusternr += 1
  wb.save(genomename + "/" + genomename + ".geneclusters.xls")
  #Gather and store data on each gene cluster
  gtrcoglist = ['SMCOG1045','SMCOG1062','SMCOG1102']
  transportercoglist = ['SMCOG1000','SMCOG1005','SMCOG1011','SMCOG1020','SMCOG1029','SMCOG1033','SMCOG1035','SMCOG1044','SMCOG1065','SMCOG1067','SMCOG1069','SMCOG1074','SMCOG1085','SMCOG1096','SMCOG1106','SMCOG1118','SMCOG1131','SMCOG1166','SMCOG1169','SMCOG1184','SMCOG1202','SMCOG1205','SMCOG1214','SMCOG1234','SMCOG1243','SMCOG1245','SMCOG1252','SMCOG1254','SMCOG1288']
  qgeneclusterdata = {}
  if smcogs == "y":
    smcogdict2 = {}
    smcogdescriptions = {}
    for i in smcogdict.keys():
      if len(smcogdict[i]) > 0 and len(smcogdict[i][0]) > 0 and ":" in smcogdict[i][0][0]:
        smcogdict2[i] = (smcogdict[i][0][0]).split(":")[0]
        smcogdescriptions[(smcogdict[i][0][0]).split(":")[0]] = (smcogdict[i][0][0]).split(":")[1]
      elif len(smcogdict[i]) > 0:
        smcogdict2[i] = smcogdict[i][0][0]
    smcogdict = smcogdict2
  for genecluster in geneclusters:
    clustergenes = clusterinfo[genecluster][4]
    clustergenes2 = []
    #for i in clustergenes:
    #  clustergenes2.append(i[4])
    [clustergenes2.append(i[4]) for i in clustergenes]
    clustergenes = clustergenes2
    clusternr = 1
    clustertype = clusterinfo[genecluster][0]
    annotations = {}
    colors = []
    starts = []
    ends = []
    strands = []
    pksnrpsprots = []
    gtrs = []
    transporters = []
    for j in clustergenes:
      annotations[j] = proteins[3][j][3]
      starts.append(int(proteins[3][j][0]))
      ends.append(int(proteins[3][j][1]))
      strands.append(proteins[3][j][2])
      if j in allcoregenes:
        colors.append("#810E15")
      else:
        colors.append("grey")
      if j in pksnrpscoregenes:
        pksnrpsprots.append(j)
      if smcogs == "y":
        if smcogdict.has_key(j) and len(smcogdict[j]) > 0 :
          if smcogdict[j][0] in gtrcoglist:
            gtrs.append(j)
          if smcogdict[j][0] in transportercoglist:
            transporters.append(j)
    clustersize = max(ends) - min(starts)
    if clusterblast == "n":
      nrhitgeneclusters = {}
      for i in geneclusters:
        nrhitgeneclusters[i] = 0
    hitgeneclusters = range(1,(nrhitgeneclusters[genecluster] + 1))
    hitgeneclusterdata = {}
    hitgeneclusterdata[genecluster] = [hitgeneclusters]
    pksnrpsprotsnames = nrpspkstypedict
    pksnrpsdomains = {}
    domlist = []
    domsdetails = {}
    substrspecnrpspredictordict = {}
    substrspecminowadict = {}
    substrspecpkssigdict = {}
    substrspecconsensusdict = {}
    krpredictionsdict = {}
    for i in pksnrpsprots:
      domlist = []
      domsdetails = {}
      doms = domaindict[i]
      for j in doms:
        nr = 1
        while j[0] + str(nr) in domlist:
          nr += 1
        domname = j[0] + str(nr)
        domlist.append(domname)
        domsdetails[domname] = [j[1],j[2]]
        if "AMP-binding" in domname or "A-OX" in domname:
          domname2 = i + "_" + "A" + str(nr)
          substrspecminowadict[domname2] = minowa_nrps_preds[i + "_A" + str(nr)]
          substrspecnrpspredictordict[domname2] = [nrps_code_preds[i + "_A" + str(nr)],nrps_svm_preds[i + "_A" + str(nr)]]
          substrspecconsensusdict[domname2] = consensuspreds[i + "_A" + str(nr)]
        if "PKS_AT" in domname:
          domname2 = i + "_" + "AT" + str(nr)
          substrspecminowadict[domname2] = minowa_pks_preds[i + "_AT" + str(nr)]
          substrspecpkssigdict[domname2] = pks_code_preds[i + "_AT" + str(nr)]
          substrspecconsensusdict[domname2] = consensuspreds[i + "_AT" + str(nr)]
        if "CAL_domain" in domname:
          domname2 = i + "_" + "CAL" + str(nr)
          substrspecminowadict[domname2] = minowa_cal_preds[i + "_CAL" + str(nr)]
          substrspecconsensusdict[domname2] = consensuspreds[i + "_CAL" + str(nr)]
        if "CAL_domain" in domname:
          domname2 = i + "_" + "CAL" + str(nr)
          substrspecminowadict[domname2] = minowa_cal_preds[i + "_CAL" + str(nr)]
          substrspecconsensusdict[domname2] = consensuspreds[i + "_CAL" + str(nr)]
        if "PKS_KR" in domname:
          domname2 = i + "_" + "KR" + str(nr)
          krpredictionsdict[domname2] = [kr_activity_preds[i + "_KR" + str(nr)],kr_stereo_preds[i + "_KR" + str(nr)]]
      pksnrpsdomains[i] = [domlist,domsdetails]
    if compound_pred_dict.has_key(genecluster):
      structpred = compound_pred_dict[genecluster]
    else:
      structpred = "N/A"
    qgeneclusterdata[genecluster] = [clustertype,clustersize,clustergenes,annotations,starts,ends,strands,pksnrpsprots,pksnrpsprotsnames,pksnrpsdomains,substrspecnrpspredictordict,substrspecminowadict,substrspecpkssigdict,substrspecconsensusdict,gtrs,transporters,colors,hitgeneclusterdata,structpred,krpredictionsdict]

  #Create genecluster svg for each gene cluster
  geneposdict = {}
  for qclusternr in geneclusters:
    data = qgeneclusterdata[qclusternr]
    #Some of the below 23 lines may already be internal to script, scan to remove unnecessary data fetching
    clustertype = data[0]
    clustersize = data[1]
    genes = data[2]
    annotations = data[3]
    starts = data[4]
    ends = data[5]
    strands = data[6]
    pksnrpsprots = data[7]
    pksnrpsprotsnames = data[8]
    pksnrpsdomains = data[9]
    substrspecnrpspredictordict = data[10]
    substrspecminowadict = data[11]
    substrspecpkssigdict = data[12]
    substrspecconsensusdict = data[13]
    gtrs = data[14]
    transporters = data[15]
    colors = data[16]
    hitgeneclusterdata = data[17]
    structpred = data[18]
    krpredictionsdict = data[19]
    relpositions = relativepositions(starts,ends,clustersize)
    rel_starts = relpositions[0]
    rel_ends = relpositions[1]
    y = 0
    for i in genes:
      geneposdict[i] = [starts[y],ends[y]]
      y += 1
    s = geneclustersvg(genes,rel_starts,rel_ends,strands,geneposdict,pksnrpsprots,pksnrpsdomains,qclusternr)
    outfile = open(svgfolder + "genecluster" + str(qclusternr) + ".svg","w")
    outfile.write(s.getXML())
    outfile.close()
  #Create ClusterBlast svg
  if clusterblast == "y":
    clusterblastpositiondata = {}
    #Create alignment svg for each pair of hit&query
    for i in geneclusters:
      hitclusters = range(queryclusterdata[i][0] + 1)[1:]
      #Create svgs for pairwise gene cluster alignment
      colorschemedict,rgbcolorscheme = calculate_colorgroups(i,hitclusters,queryclusterdata,internalhomologygroupsdict)
      for k in hitclusters:
        cresults = clusterblastresults(i,[k],queryclusterdata,colorschemedict,rgbcolorscheme)
        s = cresults[0]
        clusterblastpositiondata[str(i) + "_"+str(k)] = cresults[1]
        outfile = open(svgfolder + "clusterblast" + str(i) + "_" + str(k) + ".svg","w")
        outfile.write(s.getXML())
        outfile.close()
      #Create svgs for multiple gene cluster alignment
      cresults = clusterblastresults(i,hitclusters,queryclusterdata,colorschemedict,rgbcolorscheme)
      s = cresults[0]
      clusterblastpositiondata[str(i) + "_all"] = cresults[1]
      outfile = open(svgfolder + "clusterblast" + str(i) + "_all.svg","w")
      outfile.write(s.getXML())
      outfile.close()

  #Create folder for SEARCHGTR HTML files, load search form template
  formtemplate = open("search_form.html","r")
  formtemplate = formtemplate.read()
  formtemplate = formtemplate.replace("\r","\n")
  formtemplateparts = formtemplate.split("FASTASEQUENCE")
  #Create HTML file with gene cluster info in hidden div tags
  htmlfile = open("empty.xhtml","r")
  html = htmlfile.read()
  html = html.replace("\r","\n")
  htmlparts = html.split("<SPLIT HERE>")
  htmloutfile = open(genomename + "/display.xhtml","w")
  htmloutfile.write(htmlparts[0])
  #Add lines toreload all svgs up front
  for qclusternr in geneclusters:
    htmloutfile.write('  loadsvg(' + str(qclusternr) + ');\n')
  if clusterblast == "y":
    cblastclusters = [1,2,3,4,5,6,7,8,9,10]
    for qclusternr in geneclusters:
      nrhitclusters = queryclusterdata[qclusternr][0]
      for j in range(nrhitclusters):
        htmloutfile.write('  loadcblastsvg(' + str(qclusternr) + ',' + str(j+1) + ');\n')
  #For each gene cluster, add hidden div tags for gene names, add hidden div tags for NRPS/PKS domains, add hidden div tags for ClusterBLAST depictions
  htmloutfile.write(htmlparts[1])
  for qclusternr in geneclusters:
    data = qgeneclusterdata[qclusternr]
    pksnrpsprots = data[7]
    pksnrpsprotsnames = data[8]
    pksnrpsdomains = data[9]
    a = 0
    for i in pksnrpsprots:
      for j in pksnrpsdomains[i][0]:
        htmloutfile.write('  $("#b' + str(qclusternr) + '_00' + str(a) + '_div").hide();\n')
        a += 1
  htmloutfile.write(htmlparts[2])
  #Add top menu
  gifdict = {"t1pks":"16","t2pks":"17","t3pks":"18","t4pks":"20","nrps":"10","amglyccycl":"1","bcin":"2","blactam":"3","butyrolactone":"4","ectoine":"5","terpene":"19","indole":"7","lant":"8","melanin":"9","nucleoside":"12","other":"13","phosphoglycolipid":"14","siderophore":"15"}
  htmloutfile.write('<img border="0" align="top" src="images/empty.png" name="img0_" />\n')
  menubutton_nr = 1
  nrclustercolumns = 1
  for i in geneclusters:
    if qgeneclusterdata[i][0] in gifdict.keys():
      typenr = gifdict[qgeneclusterdata[i][0]]
    elif "-" in qgeneclusterdata[i][0]:
      typenr = "6"
    else:
      typenr = "13"
    htmloutfile.write('<a href="javascript:displaycluster(' + str(i) + ')"><img align="top" border="0" src="images/img' + str(i) + '_1.png" name="img' + str(i) + '_" onmouseover="over(' + str(i) + '),over2(0,' + typenr + ')" onmouseout="out(' + str(i) + '),out2(0,' + typenr + ')"/></a>\n')
    if menubutton_nr == 22 or menubutton_nr == 49:
      htmloutfile.write('<br/>')
      nrclustercolumns += 1
    menubutton_nr += 1

  #Add gene cluster description
  htmloutfile.write(htmlparts[3])
  extrapixelsdict = {}
  for qclusternr in geneclusters:
    data = qgeneclusterdata[qclusternr]
    clustertype = data[0]
    clustersize = data[1]
    genes = data[2]
    annotations = data[3]
    starts = data[4]
    ends = data[5]
    strands = data[6]
    pksnrpsprots = data[7]
    pksnrpsprotsnames = data[8]
    pksnrpsdomains = data[9]
    substrspecnrpspredictordict = data[10]
    substrspecminowadict = data[11]
    substrspecpkssigdict = data[12]
    substrspecconsensusdict = data[13]
    gtrs = data[14]
    transporters = data[15]
    colors = data[16]
    hitgeneclusterdata = data[17]
    structpred = data[18]
    krpredictionsdict = data[19]
    relpositions = relativepositions(starts,ends,clustersize)
    rel_starts = relpositions[0]
    rel_ends = relpositions[1]
    #Create genes overview pop-up HTMLs
    genepopupoutfile = open(htmlfolder + "geneclustergenes" + str(qclusternr) + '.html',"w")
    genepopupoutfile.write('<html>\n<head>\n<LINK href="style.css" rel="stylesheet" type="text/css">\n</head>\n<body>\nOverview of gene cluster genes:<br><br><table border=1>\n')
    genepopupoutfile.write('<tr><td><b>Gene</b></td><td><b>Annotation</b></td><td><b>Start position</b></td><td><b>End position</b></td><td><b>Strand</b></td></tr>\n')
    for i in genes:
      genepopupoutfile.write('<tr><td>' + i + '</td><td>' + annotations[i].replace("_"," ") + '</td><td>' + str(starts[genes.index(i)]) + '</td><td>' + str(ends[genes.index(i)]) + '</td><td>' + strands[genes.index(i)] +  '</td></tr>\n')
    genepopupoutfile.write('\n</table><br><br><br>Biosynthetic gene cluster signature gene domains detected: <br><br>\n')
    genepopupoutfile.write('<table border=1><tr><td><b>Gene</b></td><td><b>Detected domains</b></td><td><b>Bit scores</b></td>\n')
    for i in genes:
      if i in allcoregenes:
        detected_doms = detecteddomainsdict[i]
        for j in detected_doms:
          genepopupoutfile.write('<tr><td>' + i + '</td><td>' + str(j[0]) + '</td><td>' + str(j[1]) + '</td>\n')
    genepopupoutfile.write('\n</table><br><br><br>')
    genepopupoutfile.write('\n</body>\n</html>\n')
    genepopupoutfile.close()
    #Add gene cluster description on top
    if qclusternr == 1:
      htmloutfile.write('<div id="genecluster'+ str(qclusternr) + '">')
    else:
      htmloutfile.write('\n\n<div id="genecluster'+ str(qclusternr) + '" style="display:none">')
    #Add menu bars 1 & 2
    htmloutfile.write('<div id="bartext1" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:' + str(113 + nrclustercolumns * 28) + 'px; left:30px;"><b>Gene cluster description</b></div>')
    htmloutfile.write('<div id="bartext2" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:' + str(263 + nrclustercolumns * 28) + 'px; left:30px;"><b>PKS/NRPS domain annotation</b></div>')
    htmloutfile.write('<div id="descrbar1" style="position:absolute; z-index:1; top:' + str(110 + nrclustercolumns * 28) + 'px;"><img src="images/bar.png" height="25" width="' + str(int(0.75 * screenwidth)) + '"/></div>\n')
    htmloutfile.write('<div class="help" id="help1" style="position:absolute; z-index:1; top:' + str(112 + nrclustercolumns * 28) + 'px; left:' + str(int(screenwidth * 0.75) - 20) + 'px;"><a href="http://antismash.secondarymetabolites.org/help.html#panel1" target="_blank"><img border="0" src="images/help.png"/></a></div>\n')
    htmloutfile.write('<div id="descrbar2" style="position:absolute; z-index:1; top:' + str(260 + nrclustercolumns * 28) + 'px;"><img src="images/bar.png" height="25" width="' + str(int(0.75 * screenwidth)) + '"/></div>\n')
    htmloutfile.write('<div class="help" id="help2" style="position:absolute; z-index:1; top:' + str(262 + nrclustercolumns * 28) + 'px; left:' + str(int(screenwidth * 0.75) - 20) + 'px;"><a href="http://antismash.secondarymetabolites.org/help.html#panel2" target="_blank"><img border="0" src="images/help.png"/></a></div>\n')
    if screenwidth < 1280:
      htmloutfile.write('<div class="clusterdescr" style="font-size:0.7em; position:absolute; top:' + str(125 + nrclustercolumns * 28) + 'px; left:' + str(12) + 'px;">\n')
    else:
      htmloutfile.write('<div class="clusterdescr" style="font-size:0.8em; position:absolute; top:' + str(120 + nrclustercolumns * 28) + 'px; left:' + str(12) + 'px;">\n')
    htmloutfile.write("<br/>Gene Cluster " + str(qclusternr) + ". Type = " + clustertype + ". Location: "+ str(starts[0]) + " - " + str(ends[-1]) + " nt. Click on genes for more information.")
    if len(genomic_accnr) > 4:
      htmloutfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + genomic_accnr + '" target="_blank">GBK</a>')
    #Genes overview pop-up.
    if len(clustertype) > 20:
      htmloutfile.write('<br/>')
    htmloutfile.write('&nbsp;&nbsp;&nbsp;&nbsp;<a href="html/geneclustergenes' + str(qclusternr) + '.html" onclick=\'window.open("html/geneclustergenes' + str(qclusternr) + '.html","popup","width=800,height=800,scrollbars=yes,resizable=yes,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'>Genes and detection info overview</a>')
    htmloutfile.write("</div>\n\n")
    htmloutfile.write('<div id="display' + str(qclusternr) + '">\n')
    if nrclustercolumns > 1:
      spacers = nrclustercolumns - 1
      for i in range(spacers):
        htmloutfile.write('<img src="images/spacer.png"/>\n')
    htmloutfile.write('</div>\n')
    #Add gene pop-ups
    a = 0
    for i in genes:
      htmloutfile.write('<div id="a' + str(qclusternr) + '_00' + str(a) + '_div" class="hidden popup" style="position:absolute; z-index:2; top:' + str(185 + nrclustercolumns * 28) + 'px; left:' + str(int(((rel_starts[a] + rel_ends[a])/2)*0.875)) + 'px;">\n')
      htmloutfile.write(annotations[i].replace("_"," ").replace("&","&amp;") + "\n")
      if smcogs == "y":
        if smcogdict.has_key(i):
          smcog = smcogdict[i]
          htmloutfile.write("<br/>smCOG: " + smcog + " (" + smcogdescriptions[smcog].replace("_"," ").replace("&","&amp;") + ")\n")
          if smcog in gtrcoglist:
            formfileloc = searchgtrfolder + i + ".html"
            formfile = open(formfileloc,"w")
            specificformtemplate = formtemplateparts[0].replace("GlycTr",i)
            formfile.write(specificformtemplate)
            formfile.write(i + "\n" + seqdict[i])
            formfile.write(formtemplateparts[1])
            formfile.close()
            htmloutfile.write("<br/><a href=\"searchgtr/" + i + ".html\" target=\"_blank\"> Run SEARCHGTr on this gene </a>\n")
          if smcog in transportercoglist:
            link = "http://blast.jcvi.org/er-blast/index.cgi?project=transporter;program=blastp;sequence=sequence%0A" + seqdict[i]
            htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> TransportDB BLAST on this gene </a>\n")
        else:
          htmloutfile.write("<br/>smCOG: -\n")
      link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + seqdict[i] + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
      htmloutfile.write("<br/>Location: " + str(starts[a]) + "-" + str(ends[a]) + "\n")
      htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a><br/>\n")
      browse_start = starts[a] - 10000
      browse_end = ends[a] + 10000
      if browse_start < 0:
        browse_start = 0
      if browse_end > dnaseqlength:
        browse_end = dnaseqlength
      if genomic_accnr != "none" and genomic_accnr != "":
        htmloutfile.write('<a href="http://www.ncbi.nlm.nih.gov/projects/sviewer/?Db=gene&amp;DbFrom=protein&amp;Cmd=Link&amp;noslider=1&amp;id=' + genomic_accnr + '&amp;from=' + str(browse_start) + '&amp;to=' + str(browse_end) + '" target=\"_blank\">View genomic context</a><br/>\n')
      if smcogs == "y":
        if smcogtreedict.has_key(i.rpartition(".")[0]):
          htmloutfile.write('<a href="smcogs/' + smcogtreedict[i.rpartition(".")[0]] + '" onclick=\'window.open("smcogs/' + smcogtreedict[i.rpartition(".")[0]] + '","popup","width=1280,height=1500,resizable=yes,scrollbars=yes,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'>View smCOG seed phylogenetic tree with this gene</a>\n')
        elif smcogtreedict.has_key(i):
          htmloutfile.write('<a href="smcogs/' + smcogtreedict[i] + '" onclick=\'window.open("smcogs/' + smcogtreedict[i] + '","popup","width=1280,height=1500,resizable=yes,scrollbars=yes,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'>View smCOG seed phylogenetic tree with this gene</a>\n')
      htmloutfile.write("</div>\n\n")
      htmloutfile.write('<div id="a' + str(qclusternr) + '_00' + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(162 + nrclustercolumns * 28) + 'px; left:' + str(float((rel_starts[a]+rel_ends[a])/2)*0.9375) + 'px;">\n')
      htmloutfile.write(i)
      htmloutfile.write("</div>\n\n")
      a += 1
    #Early calculation of nr of domains to be able to fit structure prediction information of large NRPSs/PKSs
    pksnrpsdomainnr = 0
    krdomainnr = 0
    adomainnr = 0
    for i in pksnrpsprots:
      doms = pksnrpsdomains[i][0]
      first = "no"
      nra = 0
      nrat = 0
      nrkr = 0
      nrcal = 0
      for j in doms:
        if "AMP-binding" in j or "A-OX" in j:
          j = "A"
          nra += 1
          adomainnr += 1
          z = nra
        if "KR" in j:
          j = "KR"
          nrkr += 1
          krdomainnr += 1
          z = nrkr
        if "AT" in j and "docking" not in j:
          j = "AT"
          nrat += 1
          pksnrpsdomainnr += 1
          z = nrat
        if "CAL" in j:
          j = "CAL"
          nrcal += 1
          pksnrpsdomainnr += 1
          z = nrcal
    pixels = adomainnr * 50  + pksnrpsdomainnr * 40 + krdomainnr * 30 + (len(pksnrpsprots) * 16) + 375
    extrapixels = pixels - (676 + len(pksnrpsprots) * 99)
    if extrapixels < 0:
      extrapixels = 0
    extrapixelsdict[qclusternr] = extrapixels
    #Add picture of predicted chemical structure
    htmloutfile.write('<div id="verticalbar1" style="position:absolute; left:' + str(int(screenwidth * 0.75) + 12) + 'px; top:' + str(106 + nrclustercolumns * 28) + 'px;"><img src="images/linefill.png" height="' + str(1126 + len(pksnrpsprots) * 99 + extrapixels) + '" width="2"/></div>\n')
    htmloutfile.write('<div id="verticalbar2" style="position:absolute; left:' + str(int(screenwidth * 0.98)) + 'px; top:0px;"><img src="images/linefill.png" height="' + str(1288 + len(pksnrpsprots) * 99 + nrclustercolumns * 28 + extrapixels) + '" width="2"/></div>\n')
    htmloutfile.write('<div id="horizbar1" style="position:absolute; left:0px; top:' + str(92 + nrclustercolumns * 28) + 'px;"><img src="images/linefill.png" height="2" width="' + str(screenwidth * 0.98) + '"/></div>\n')
    htmloutfile.write('<div id="horizbar2" style="position:absolute; left:0px; top:82px;"><img src="images/linefill.png" height="2" width="' + str(screenwidth * 0.98) + '"/></div>\n')
    htmloutfile.write('<div id="horizbar3" style="position:absolute; left:0px; top:' + str(1223 + len(pksnrpsprots) * 99 + nrclustercolumns * 28 + extrapixels) + 'px;"><img src="images/linefill.png" height="2" width="' + str(screenwidth * 0.98) + '"/></div>\n')
    if screenwidth < 1280:
      htmloutfile.write('<div id="bartext4" style="color:#FFFFFF; font-size:0.8em; position:absolute; z-index:2; top:' + str(114 + nrclustercolumns * 28) + 'px; left:' + str(int(screenwidth * 0.75) + 30) + 'px;"><b>Predicted core structure</b></div>\n')
    else:
      htmloutfile.write('<div id="bartext4" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:' + str(113 + nrclustercolumns * 28) + 'px; left:' + str(int(screenwidth * 0.75) + 30) + 'px;"><b>Predicted core structure</b></div>\n')
    htmloutfile.write('<div class="title" style="position:absolute; top:' + str(110 + nrclustercolumns * 28) + 'px; left:' + str(screenwidth * 0.75 + 20) + 'px;">\n')
    htmloutfile.write('<div id="descrbar4" style="right:25px; position:absolute; z-index:1; top:0px; left:0px;"><img src="images/bar.png" height="25" width="' + str(int(0.21 * screenwidth)) + '"/></div>\n')
    htmloutfile.write('<div class="help" id="help4" style="position:absolute; z-index:1; top:2px; left:' + str(int(screenwidth * 0.2) - 20) + 'px;"><a href="http://antismash.secondarymetabolites.org/help.html#sidepanel1" target="_blank"><img border="0" src="images/help.png"/></a></div>\n')
    if qclusternr in failedstructures:
      htmloutfile.write('<br/><br/><img src="images/nostructure_icon.png" border="1" width="' + str(int(screenwidth * 0.19)) + '" height="200" />\n')  
    elif " " in structpred:
      htmloutfile.write('<br/><br/><a href="structures/genecluster' + str(qclusternr) + '.png" onclick=\'window.open("structures/genecluster' + str(qclusternr) + '.png","popup","width=600,height=300,scrollbars=yes,resizable=yes,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'><img src="structures/genecluster' + str(qclusternr) + '_icon.png" border="1" width="' + str(int(screenwidth * 0.19)) + '" height="200" /></a>\n')
    else:
      htmloutfile.write('<br/><br/><img src="images/nostructure_icon.png" border="1" width="' + str(int(screenwidth * 0.19)) + '" height="200" />\n')
    htmloutfile.write('<div class="clusterdescr" style="font-size:0.8em;">\n')
    htmloutfile.write("Monomers prediction: " + structpred + "<br/>\n")
    if qclusternr in dockingdomainanalysis:
      htmloutfile.write('<a href="html/docking_analysis' + str(qclusternr) + '.html" onclick=\'window.open("html/docking_analysis' + str(qclusternr) + '.html","popup","width=600,height=1200,scrollbars=yes,resizable=yes,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'>Docking domain analysis results.</a><br/>\n')
    nrpsfound = "no"
    pksnrpsdomainnr = 0
    adomainnr = 0
    krdomainnr = 0
    for i in pksnrpsprots:
      doms = pksnrpsdomains[i][0]
      first = "no"
      nra = 0
      nrat = 0
      nrkr = 0
      nrcal = 0
      for j in doms:
        if "AMP-binding" in j or "A-OX" in j:
          j = "A"
          nra += 1
          adomainnr += 1
          z = nra
        if "KR" in j:
          j = "KR"
          nrkr += 1
          krdomainnr += 1
          z = nrkr
        if "AT" in j and "docking" not in j:
          j = "AT"
          nrat += 1
          pksnrpsdomainnr += 1
          z = nrat
        if "CAL" in j:
          j = "CAL"
          nrcal += 1
          pksnrpsdomainnr += 1
          z = nrcal
        prediction = "no"
        domname = str(i) + "_" + str(j) + str(z)
        if domname in substrspecnrpspredictordict.keys():
          nrpsfound = "yes"
          prediction = "yes"
          if substrspecnrpspredictordict[domname][0] == "nrp":
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;NRPSPredictor code prediction, '+ str(j) + str(z) + ': ?</font><br/>\n')
          else:
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;NRPSPredictor code prediction, '+ str(j) + str(z) + ': ' + substrspecnrpspredictordict[domname][0] + '</font><br/>\n')
          if substrspecnrpspredictordict[domname][1] == "nrp":
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;NRPSPredictor SVM prediction, '+ str(j) + str(z) + ': ?</font><br/>\n')
          else:
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;NRPSPredictor SVM prediction, '+ str(j) + str(z) + ': ' + substrspecnrpspredictordict[domname][1] + '</font><br/>\n')
        if domname in substrspecminowadict.keys():
          prediction = "yes"
          if substrspecminowadict[domname] == "nrp" or substrspecminowadict[domname] == "pk":
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;Minowa prediction, '+ str(j) + str(z) + ': ?</font><br/>\n')
          else:
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;Minowa prediction, '+ str(j) + str(z) + ': ' + substrspecminowadict[domname] + '</font><br/>\n')
        if domname in substrspecpkssigdict.keys():
          prediction = "yes"
          if substrspecpkssigdict[domname] == "pk":
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;PKS code prediction, '+ str(j) + str(z) + ': ?</font><br/>\n')
          else:
            if first == "no":
              first = "yes"
              htmloutfile.write(i + ':<br/>')  
            htmloutfile.write('<font size="1">&nbsp;&nbsp;PKS code prediction, '+ str(j) + str(z) + ': ' + substrspecpkssigdict[domname] + '</font><br/>\n')
        if domname in krpredictionsdict.keys():
          if first == "no":
            first = "yes"
            htmloutfile.write(i + ':<br/>')  
          htmloutfile.write('<font size="1">&nbsp;&nbsp;KR activity, '+ str(j) + str(z) + ': ' + krpredictionsdict[domname][0] + "</font><br/>\n")
          htmloutfile.write('<font size="1">&nbsp;&nbsp;KR stereochemistry, '+ str(j) + str(z) + ': ' + krpredictionsdict[domname][1] + "</font><br/>\n")
        #Add link to prediction details pop-up
        if prediction == "yes":
          htmloutfile.write('<font size="1">&nbsp;&nbsp;&nbsp;&nbsp;<a href="substrspecs/' + domname + '.html" onclick=\'window.open("substrspecs/' + domname + '.html","popup","width=500,height=400,scrollbars=yes,resizable=no,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'>Prediction details</a></font><br/>\n')
    if nrpsfound == "yes":
      htmloutfile.write('<br/><a href="http://bioinfo.lifl.fr/norine/form2.jsp" target="_blank">Perform Norine peptide search</a>')
    htmloutfile.write('</div>')
    if screenwidth < 1280:
      htmloutfile.write('<div id="bartext5" style="color:#FFFFFF; font-size:0.8em; position:absolute; z-index:2; top:' + str(624 + adomainnr * 50 + pksnrpsdomainnr * 40 + krdomainnr * 30 + (len(pksnrpsprots) * 16)) + 'px; left:10px;"><b>File outputs</b></div>\n')
    else:
      htmloutfile.write('<div id="bartext5" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:' + str(623 + adomainnr * 50  + pksnrpsdomainnr * 40 + krdomainnr * 30 + (len(pksnrpsprots) * 16)) + 'px; left:10px;"><b>Downloadable output files</b></div>\n')
    htmloutfile.write('<div id="descrbar5" style="right:25px; position:absolute; z-index:1; top:' + str(620 + adomainnr * 50  + pksnrpsdomainnr * 40 + krdomainnr * 30 + (len(pksnrpsprots) * 16)) + 'px; left:0px;"><img src="images/bar.png" height="25" width="' + str(int(0.21 * screenwidth)) + '"/></div>\n')
    htmloutfile.write('<div class="help" id="help5" style="position:absolute; z-index:1; top:' + str(622 + adomainnr * 50  + pksnrpsdomainnr * 40 + krdomainnr * 30 + (len(pksnrpsprots) * 16)) + 'px; left:' + str(int(screenwidth * 0.2) - 20) + 'px;"><a href="http://antismash.secondarymetabolites.org/help.html#sidepanel2" target="_blank"><img border="0" src="images/help.png"/></a></div>\n')
    htmloutfile.write('<div class="text" id="outputinfo" style="font-size:0.8em; right:25px; position:absolute; z-index:1; top:' + str(655 + adomainnr * 50  + pksnrpsdomainnr * 40 + krdomainnr * 30 + (len(pksnrpsprots) * 16)) + 'px; left:0px;">')
    if fullhmm == "y" or fullblast == "y":
      htmloutfile.write('<a href="' + oldgenomename + '.final.embl" target="_blank">Open EMBL summary file</a><br/><br/>')
    #htmloutfile.write('<a href="' + genomename + '.final.csv" target="_blank">Download CSV summary file</a><br/><br/>')
    if fullhmm == "y":
      htmloutfile.write('<a href="' + oldgenomename + '.cluster_prediction.png" onclick=\'window.open("' + oldgenomename + '.cluster_prediction.png","popup","width=1024,height=1400,scrollbars=0,resizable=0,toolbar=0,directories=0,location=0,menubar=0,status=0,left=0,top=0"); return false\'>Sec. met. enriched genome regions</a><br/><br/>')
    htmloutfile.write('<a href="' + genomename + '.geneclusters.xls" target="_blank">Open XLS overview table</a><br/><br/>')
    htmloutfile.write('</div>')
    htmloutfile.write("</div>\n\n")
    #Add descriptions of NRPS/PKS genes
    htmloutfile.write('<div class="title" style="position:absolute; top:' + str(180) + 'px; left:' + str(12) + 'px;">\n')
    htmloutfile.write("</div>\n\n")
    z = 1
    for i in pksnrpsprots:
      htmloutfile.write('<div class="text" style="position:absolute; top:' + str(228 + 84 * z + nrclustercolumns * 28) + 'px; left:' + str(12) + 'px;">\n')
      htmloutfile.write(i + " (" + pksnrpsprotsnames[i].lower() + ")")
      htmloutfile.write("</div>\n\n")
      z += 1
    #Add NRPS/PKS domain pop-ups
    longestprot = 0
    protlengthdict = {}
    for i in pksnrpsprots:
        protlength = (geneposdict[i][1] - geneposdict[i][0]) / 3
        protlengthdict[i] = protlength
        if protlength > longestprot:
          longestprot = protlength
    try:
        aa2pixelratio = longestprot * 0.75 / screenwidth
    except:
        aa2pixelratio = 0.1
    a = 0
    z = 1
    for i in pksnrpsprots:
      domainsdict = pksnrpsdomains[i][1]
      nra = 0
      nrat = 0
      nrkr = 0
      nrcal = 0
      for j in pksnrpsdomains[i][0]:
        startpos = domainsdict[j][0]
        endpos = domainsdict[j][1]
        htmloutfile.write('<div id="b' + str(qclusternr) + '_00' + str(a) + '_div" class="hidden popup" style="position:absolute; z-index:2; top:' + str(277 + 84 * z + nrclustercolumns * 28) + 'px; left:' + str( ( ( (endpos+startpos) / 2) / aa2pixelratio) * 0.9375 ) + 'px;">\n')
        htmloutfile.write("Domain " + j + " (" + i + ")")
        link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + seqdict[i][startpos:endpos] + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
        htmloutfile.write("<br/>Location: " + str(startpos) + "-" + str(endpos) + " AA\n")
        domid = i + "_" + j
        if "AMP-binding" in j or "A-OX" in j:
          j = "A"
          nra += 1
          y = nra
        if "PKS_KR" in j:
          j = "KR"
          nrkr += 1
          y = nrkr
        if "PKS_AT" in j:
          j = "AT"
          nrat += 1
          y = nrat
        if "CAL_domain" in j:
          j = "CAL"
          nrcal += 1
          y = nrcal
        prediction = "no"
        domid = str(i) + "_" + str(j) + str(y)
        if substrspecnrpspredictordict.has_key(domid) or substrspecminowadict.has_key(domid) or substrspecpkssigdict.has_key(domid):
          htmloutfile.write("<br/>Predicted substrate: " + substrspecconsensusdict[domid] + "\n")
        if substrspecnrpspredictordict.has_key(domid):
          htmloutfile.write("<br/>-NRPSPredictor code: " + substrspecnrpspredictordict[domid][0] + "\n")
          htmloutfile.write("<br/>-NRPSPredictor SVM: " + substrspecnrpspredictordict[domid][1] + "\n")
        if substrspecminowadict.has_key(domid):
          htmloutfile.write("<br/>-Minowa HMM: " + substrspecminowadict[domid] + "\n")
        if substrspecpkssigdict.has_key(domid):
          htmloutfile.write("<br/>-PKS code: " + substrspecpkssigdict[domid] + "\n")
        if krpredictionsdict.has_key(domid):
          htmloutfile.write("<br/>KR activity: " + krpredictionsdict[domid][0] + "\n")
          htmloutfile.write("<br/>KR stereochemistry: " + krpredictionsdict[domid][1] + "\n")
        htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this domain </a>\n")
        htmloutfile.write("</div>\n\n")
        a += 1
      z += 1 
    htmloutfile.write('</div>\n')

  if clusterblast == "y":
    #Write ClusterBlast divs with pictures and description pop-up tags
    htmloutfile.write('<div id="clusterblastview" class="clusterdescr">\n\n')
    #Add menu bar 3
    htmloutfile.write('<div id="bartext3" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:3px; left:20px;"><b>Homologous gene clusters</b></div>')
    htmloutfile.write('<div id="descrbar3" style="position:absolute; z-index:1; top:0px;"><img src="images/bar.png" height="25" width="' + str(int(0.75*screenwidth)) + '"/></div>')
    htmloutfile.write('<div class="help" id="help3" style="position:absolute; z-index:1; top:2px; left:' + str(int(screenwidth * 0.75) - 30) + 'px;"><a href="http://antismash.secondarymetabolites.org/help.html#panel3" target="_blank"><img border="0" src="images/help.png"/></a></div>')
    for qclusternr in geneclusters:
      nrhitclusters = queryclusterdata[qclusternr][0]
      hitclusterdata = queryclusterdata[qclusternr][1]
      if qclusternr == 1:
        htmloutfile.write('<div id="qcluster' + str(qclusternr) + '">\n<br/><br/>\n<div align="left">\n<form name="clusterform' + str(qclusternr) + '">\n<select name="selection' + str(qclusternr) + '" onchange="javascript:navigate(this);">\n')
      else:
        htmloutfile.write('<div id="qcluster' + str(qclusternr) + '" style="display:none">\n<br/><br/>\n<div align="left">\n<form name="clusterform' + str(qclusternr) + '">\n<select name="selection' + str(qclusternr) + '" onchange="javascript:navigate(this);">\n')
      htmloutfile.write('<option value="">Select gene cluster alignment</option>\n')
      for i in range(nrhitclusters):
        htmloutfile.write('<option value="javascript:displaycblastresults(' + str(qclusternr) + ',' + str(i+1) + ')">' + hitclusterdata[i+1][5][i].replace("&","&amp;") + '</option>\n')
      htmloutfile.write('</select>\n</form>\n\n</div>')
      htmloutfile.write('<div style="position:absolute; top:33px; left:' + str(screenwidth*0.625) + 'px;"><img src="images/button.gif" name="button' + str(qclusternr) + '" onclick="javascript:displaybutton(' + str(qclusternr) + ');"/></div>')
      clustersizes = []
      for i in range(nrhitclusters):
        hitclusterdata = queryclusterdata[qclusternr][1]
        queryclustergenes = hitclusterdata[1][3]
        queryclustergenesdetails = hitclusterdata[1][4]
        hitclusternumber =  i + 1
        cluster_acc = hitclusterdata[hitclusternumber][6]
        hitclustergenes = hitclusterdata[hitclusternumber][1]
        hitclustergenesdetails = hitclusterdata[hitclusternumber][2]
        relpositiondata = clusterblastpositiondata[str(qclusternr) + "_" + str(i+1)]
        qrel_starts = relpositiondata[0][0]
        qrel_ends = relpositiondata[0][1]
        hrel_starts = relpositiondata[1][hitclusternumber ][0]
        hrel_ends = relpositiondata[1][hitclusternumber ][1]
        strandsbalance = relpositiondata[2][hitclusternumber]
        if strandsbalance < 0:
          hitclustergenes.reverse()
        if qclusternr == 1 and (i+1) == 1:
          htmloutfile.write('<div id="hitcluster' + str(qclusternr) + '_' + str(i+1) + '">\n')
        else:
          htmloutfile.write('<div id="hitcluster' + str(qclusternr) + '_' + str(i+1) + '" style="display:none">\n')
        #Insert gene cluster descriptions
        cdescription = hitclusterdata[i+1][5][i].replace("&","&amp;").replace("\t"," ").partition(" ")[2].partition(" ")[2].split(", whole")[0].split(", complete")[0]
        if len(nucname) < 80:
          qdescription = nucname
        else:
          qdescription = nucname[0:77] + "..."
        htmloutfile.write('<div id="descriptionquery" style="text-align:right; position:absolute; top:70px; right:50px; font-size:10px; font-style:italic">' + qdescription + '</div>\n')
        htmloutfile.write('<div id="description' + str(qclusternr) + '" style="text-align:right; position:absolute; top:137px; right:50px; font-size:10px; font-style:italic">' + cdescription + '</div>\n')
        #Insert pubmed/pubchem links
        htmloutfile.write('<div id="pub_pics" style="position:absolute; top:60px; left:' + str(int(screenwidth * 0.0)) + 'px; font-size:10px"> Hit cluster cross-links: \n')
        htmloutfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.split(".")[0] + '" target="_blank"><img align="bottom" border="0" src="images/genbank.gif"/></a>\n')
        present = "n"
        for j in pubmed_dict.keys():
          if j in cluster_acc:
            present = "y"
        for j in pubchem_dict.keys():
          if j in cluster_acc:
            present = "y"
        if present == "y":
          for j in pubmed_dict.keys():
            if j in cluster_acc:
              pubmedstring = pubmed_dict[j]
              htmloutfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/pubmed/' + pubmedstring + '" target="_blank"><img align="bottom" border="0" src="images/pubmed.gif"/></a>\n')
          for j in pubchem_dict.keys():
            if j in cluster_acc:
              pubchemstring = pubchem_dict[j]
              if "," in pubchemstring:
                htmloutfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&amp;term=' + pubchemstring + '" target="_blank"><img align="bottom" border="0" src="images/struct.gif"/></a>\n')
              else:
                htmloutfile.write('&nbsp;&nbsp;<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=' + pubchemstring + '" target="_blank"><img align="bottom" border="0" src="images/struct.gif"/></a>\n')
        htmloutfile.write('</div>\n\n')
        #Create gene pop-ups
        a = 0
        for j in queryclustergenes:
          j_accession = accessiondict[j]
          htmloutfile.write('<div id="q' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(113) + 'px; left:' + str(int(float(qrel_starts[a])*0.875)) + 'px;">\n')
          htmloutfile.write(queryclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
          link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j_accession + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
          htmloutfile.write("<br/>Location: " + str(queryclustergenesdetails[j][0]) + "-" + str(queryclustergenesdetails[j][1]) + "\n")
          htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
          htmloutfile.write("</div>\n\n")
          htmloutfile.write('<div id="q' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(83) + 'px; left:' + str(int(float((float(qrel_starts[a])+float(qrel_ends[a]))/2)*0.9375)) + 'px;">\n')
          htmloutfile.write(j)
          htmloutfile.write("</div>\n\n")        
          a+= 1
        a = 0
        for j in hitclustergenes:
          j_accession = hitclustergenesdetails[j][4]
          htmloutfile.write('<div id="h' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(183) + 'px; left:' + str(int(float(hrel_starts[a])*0.875)) + 'px;">\n')
          htmloutfile.write(hitclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
          link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j_accession + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
          htmloutfile.write("<br/>Location: " + str(hitclustergenesdetails[j][0]) + "-" + str(hitclustergenesdetails[j][1]) + "\n")
          htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
          htmloutfile.write("</div>\n\n")
          htmloutfile.write('<div id="h' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(153) + 'px; left:' + str(int(float((float(hrel_starts[a])+float(hrel_ends[a]))/2)*0.9375)) + 'px;">\n')
          htmloutfile.write(j)
          htmloutfile.write("</div>\n\n")    
          a += 1
        htmloutfile.write('</div>\n')
      #Find new relative positions for display of all gene clusters in one picture
      relpositiondata = clusterblastpositiondata[str(qclusternr) + "_all"]
      qrel_starts = relpositiondata[0][0]
      qrel_ends = relpositiondata[0][1]
      htmloutfile.write('<div id="hitcluster' + str(qclusternr) + '_all" style="display:none">\n')
      if len(nucname) < 80:
        qdescription = nucname
      else:
        qdescription = nucname[0:77] + "..."
      htmloutfile.write('<div id="descriptionquery" style="text-align:right; position:absolute; top:60px; right:50px; font-size:10px; font-style:italic">' + qdescription + '</div>\n')
      for i in range(nrhitclusters):
        hitclusterdata = queryclusterdata[qclusternr][1]
        queryclustergenes = hitclusterdata[1][3]
        queryclustergenesdetails = hitclusterdata[1][4]
        hitclusternumber =  i + 1
        hrel_starts = relpositiondata[1][hitclusternumber][0]
        hrel_ends = relpositiondata[1][hitclusternumber][1]
        cluster_acc = hitclusterdata[hitclusternumber][6]
        hitclustergenes = hitclusterdata[hitclusternumber][1]
        hitclustergenesdetails = hitclusterdata[hitclusternumber][2]
        strandsbalance = relpositiondata[2][hitclusternumber]
        cdescription = hitclusterdata[i+1][5][i].replace("&","&amp;").replace("\t"," ").partition(" ")[2].partition(" ")[2].split(", whole")[0].split(", complete")[0]
        htmloutfile.write('<div id="description' + str(qclusternr) + '" style="text-align:right; position:absolute; top:' + str(60 + (57 * hitclusternumber)) + 'px; right:50px; font-size:10px; font-style:italic">' + cdescription + '</div>\n')
        if hitclusternumber == 1:
          a = 0
          for j in queryclustergenes:
            htmloutfile.write('<div id="all_' + str(qclusternr) + "_0_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(100) + 'px; left:' + str(int(float(qrel_starts[a])*0.875)) + 'px; z-index:2;">\n')
            htmloutfile.write(queryclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
            link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
            htmloutfile.write("<br/>Location: " + str(queryclustergenesdetails[j][0]) + "-" + str(queryclustergenesdetails[j][1]) + "\n")
            htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
            htmloutfile.write("</div>\n\n")
            htmloutfile.write('<div id="all_' + str(qclusternr) + "_0_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(75) + 'px; left:' + str(int(float((float(qrel_starts[a])+float(qrel_ends[a]))/2)*0.9375)) + 'px;">\n')
            htmloutfile.write(j)
            htmloutfile.write("</div>\n\n")        
            a+= 1
        a = 0
        for j in hitclustergenes:
          htmloutfile.write('<div id="all_' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(100 + 57 * hitclusternumber) + 'px; left:' + str(int(float(hrel_starts[a])*0.875)) + 'px; z-index:2;">\n')
          htmloutfile.write(hitclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
          link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
          htmloutfile.write("<br/>Location: " + str(hitclustergenesdetails[j][0]) + "-" + str(hitclustergenesdetails[j][1]) + "\n")
          htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
          htmloutfile.write("</div>\n\n")
          htmloutfile.write('<div id="all_' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(75 + 56.75 * hitclusternumber) + 'px; left:' + str(int(float((float(hrel_starts[a])+float(hrel_ends[a]))/2)*0.9375)) + 'px;">\n')
          htmloutfile.write(j)
          htmloutfile.write("</div>\n\n")    
          a += 1
      htmloutfile.write('</div>\n')
      htmloutfile.write('</div>\n\n')
  if clusterblast == "y":
    htmloutfile.write('</div>\n')
  for i in geneclusters:
    data = qgeneclusterdata[i]
    extrapixels = extrapixelsdict[i]
    pksnrpsprots = data[7]
    if i == 1:
      htmloutfile.write('<div id="creditsbar' + str(i) + '" class="banner" style="position:absolute; width:' + str(int(0.98 * screenwidth)) +'px; align:\'left\'; height:75; top:' + str(1242 + int(len(pksnrpsprots) * 99) + nrclustercolumns * 28 + extrapixels) + 'px; left:0px; color:#810E15; z-index:-1;">')
    else:
      htmloutfile.write('<div id="creditsbar' + str(i) + '" class="banner" style="display:none; position:absolute; width:' + str(int(0.98 * screenwidth)) +'px; align:\'left\'; height:75; top:' + str(1242 + int(len(pksnrpsprots) * 99) + nrclustercolumns * 28 + extrapixels) + 'px; left:0px; color:#810E15; z-index:-1;">')
    htmloutfile.write('<div style="float:center; font-size:0.9em;">\n<div style="position:absolute; top:0px; left:30px;">\n<img src="images/ruglogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n<img src="images/gbblogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n<img src="images/tueblogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n<img src="images/ucsflogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n</div>\n<div style="position:absolute; top:0px; left:600px;">\nantiSMASH: Rapid identification, annotation and analysis of secondary metabolite biosynthesis gene clusters.\n<br/>Marnix H. Medema, Kai Blin, Peter Cimermancic, Victor de Jager, Piotr Zakrzewski, Michael A. Fischbach, Tilmann Weber, Rainer Breitling &amp; Eriko Takano\n<br/><i>Nucleic Acids Research</i> (2011), proposal submitted.\n</div>\n</div>\n</div>')
  #Add final part of HTML file
  htmloutfile.write(htmlparts[-1])
  #Copy accessory files for HTML viewing
  if sys.platform == ('win32'):
    copycommand1 = "copy/y vis\\* " + genomename + " > nul"
    copycommand2 = "copy/y vis\\html\\* " + genomename + "\\html > nul"
    copycommand3 = "copy/y vis\\images\\* " + genomename + "\\images > nul"
  elif sys.platform == ('linux2'):
    copycommand1 = "cp -r vis/* " + genomename + " > /dev/null"
    copycommand2 = "true"
    copycommand3 = "true"
  os.system(copycommand1)
  os.system(copycommand2)
  os.system(copycommand3)

  #Generate EMBL output
  emblfile = open(genomename + "/embl_lines.txt","w")
  for i in geneclustergenes:
    emblfile.write(i + "\t")
    if smcogs == "y":
      if smcogdict.has_key(i):
        emblfile.write("smCOG: " + smcogdict[i] + ":" + smcogdescriptions[smcogdict[i]] + "\t")
    if nrpspkstypedict.has_key(i):
      emblfile.write("NRPS/PKS type: " + nrpspkstypedict[i] + "\t")
    if domaindict.has_key(i):
      domains = domaindict[i]
      for j in domains:
        emblfile.write(j[0] + " (" + str(j[1]) + "-" + str(j[2]) + "); E-value:" + str(j[3]) + "; Bit score: " + str(j[4]) + "\t")
      nrat = 0
      for k in minowa_pks_preds.keys():
        if i in k:
          nrat += 1
          emblfile.write("AT-domain " + str(nrat) + " Minowa substrate specificity prediction: " + minowa_pks_preds[k] + "\t")
      nrat = 0
      for k in pks_code_preds.keys():
        if i in k:
          nrat += 1
          emblfile.write("AT-domain " + str(nrat) + " PKS code substrate specificity prediction: " + pks_code_preds[k] + "\t")
      nrcal = 0
      for k in minowa_cal_preds.keys():
        if i in k:
          nrcal += 1
          emblfile.write("CAL-domain " + str(nrcal) + " Minowa substrate specificity prediction: " + minowa_cal_preds[k] + "\t")
      nra = 0
      for k in minowa_nrps_preds.keys():
        if i in k:
          nra += 1
          emblfile.write("A-domain " + str(nra) + " Minowa substrate specificity prediction: " + minowa_nrps_preds[k] + "\t")
      nra = 0
      for k in nrps_code_preds.keys():
        if i in k:
          nra += 1
          emblfile.write("A-domain " + str(nra) + " Stachelhaus code substrate specificity prediction: " + nrps_code_preds[k] + "\t")
      nra = 0
      for k in nrps_svm_preds.keys():
        if i in k:
          nra += 1
          emblfile.write("A-domain " + str(nra) + " NRPSPredictor2 SVM substrate specificity prediction: " + nrps_svm_preds[k] + "\t")
      nrkr = 0
      for k in kr_activity_preds.keys():
        if i in k:
          nrkr += 1
          emblfile.write("KR-domain " + str(nrat) + " activity prediction: " + kr_activity_preds[k] + "\t")
          emblfile.write("KR-domain " + str(nrat) + " predicted stereochemistry group: " + kr_stereo_preds[k] + "\t")
      if motifdict.has_key(i):
        l = motifdict[i]
        for m in l:
          emblfile.write("Motif " + str(m[0]) + " (" + str(m[1]) + "-" + str(m[2]) + "). E-value: " + str(m[3]) + "; Bit score: " + str(m[4]) + "\t")
    emblfile.write("\n")
  emblfile.write("\n\n>>\n\n")
  #enter separate domain entries
  for i in geneclustergenes:
    strand = strandsdict[i]
    startpos = geneposdict[i][0]
    endpos = geneposdict[i][1]
    if domaindict.has_key(i):
      domains = domaindict[i]
      for j in domains:
        if strand == "+":
          emblfile.write("misc_feature\t" + str(startpos + j[1] * 3) + ".." + str(startpos + j[2] * 3) + "\t" + str(j[0]) + " domain;\tE-value: " + str(j[3]) + "\tBit score: " + str(j[4]) + "\t/colour=2\n")
        elif strand == "-":
          emblfile.write("misc_feature\tcomplement(" + str(endpos - j[2] * 3) + ".." + str(endpos - j[1] * 3) + ")\t" + str(j[0]) + "domain;\tE-value: " + str(j[3]) + "Bit score: " + str(j[4]) + "\t/colour=2\n")
      if motifdict.has_key(i):
        l = motifdict[i]
        for m in l:
          if strand == "+":
            emblfile.write("misc_feature\t" + str(startpos + m[1] * 3) + ".." + str(startpos + m[2] * 3) + "\t" + str(m[0]) + " motif;\tE-value: " + str(m[3]) + "\tBit score: " + str(m[4]) + "\t/colour=6\n")
          elif strand == "-":
            emblfile.write("misc_feature\tcomplement(" + str(endpos - m[2] * 3) + ".." + str(endpos - m[1] * 3) + ")\t" + str(m[0]) + " motif;\tE-value: " + str(m[3]) + "\tBit score: " + str(m[4]) + "\t/colour=6\n")
  emblfile.write("\n\n>>\n\n")
  for i in geneclusters:
    cstart = clusterinfo[i][1]
    if cstart == 0:
      cstart = 1
    cend = clusterinfo[i][2]
    emblfile.write("misc_feature\t" + str(cstart) + ".." + str(cend) + "\t" + clusterinfo[i][0] + " gene cluster\t/colour=13\n")
  emblfile.close()

  #Close open html file
  htmloutfile.close()

  #Run whole-genome BLAST / HMM CLUSEAN modules & ClusterFinder
  if sys.platform == ('win32'):
    copycommand = "copy " + infile + " " + genomename + ' > nul'
  if sys.platform == ('linux2'):
    copycommand = "cp " + infile + " " + genomename
  os.system(copycommand)
  os.chdir(genomename)
  args = "--cpus %s " % nrcpus
  if fullblast == "n":
    args += "--without-blast "
  if fullhmm == "n":
    args += "--without-hmmer "
  if fullhmm == "y":
    args += '--pfamdbpath %s ' % pfamdbpath
  if fullblast == "y":
    args += '--blastdbpath %s ' % blastdbpath
  logfile.write("Running CLUSEAN pipeline modules.\n")
  if sys.platform == ('win32'):
    os.system("python ..\\clusean\\scripts\\runPipeline.py %s" % args)
  if sys.platform == ('linux2'):
    os.system( antismash_path + "clusean/scripts/runPipeline.py %s" % args)
    #print antismash_path + "clusean/scripts/runPipeline.py %s" % args

  os.chdir('..')

  #Close log file
  logfile.write("antiSMASH successfully finished in " + str(elapsed) + " seconds.\n")
  #print "antiSMASH successfully finished in " + str(elapsed) + " seconds.\n"
  logfile.close()
