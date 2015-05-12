import numpy
import sys
import random
import subprocess
import re
import decimal
import math
import os
import shutil
import time
import types
import argparse
#from argparse import RawTextHelpFormatter

#############################
# FUNCTIONS

def print2file(f, i, m):
	"""
		print content i to file f in mode m
	"""
	line = str(i)
	if m == "a":
		call = "echo \"" +  line + "\" >> " + f
	elif m == "w":
		call = "echo \"" +  line + "\" > " + f
	os.system(call)

# checking and correcting the alphabet of the constraint sequence
def checkSequenceConstraint(SC):
	"""
		Checks the Sequence constraint for illegal nucleotide characters
	"""
	out = ""
	for c in SC:
		c = c.upper()
		if c not in "ACGURYSWKMBDHVN": 
# and c!= "R" and c != "Y" and c != "S" and c != "W" and c != "K" and c != "M" and c != "B" and c != "D" and c != "H" and c != "V":
			if c == "T":
				c = "U"
			else:
				print "\tIllegal Character in the constraint sequence!"
				print "\tPlease use the IUPAC nomenclature for defining nucleotides in the constraint sequence!"
				print "\tA   	Adenine"
				print "\tC   	Cytosine"
				print "\tG   	Guanine"
				print "\tT/U 	Thymine/Uracil"
				print "\tR 	A or G"
				print "\tY 	C or T/U"
				print "\tS 	G or C"
				print "\tW 	A or T/U"
				print "\tK 	G or T/U"
				print "\tM 	A or C"
				print "\tB 	C or G or T/U"
				print "\tD 	A or G or T/U"
				print "\tH 	A or C or T/U"
				print "\tV 	A or C or G"
				print "\tN	any base"
				exit(0)
		out += c
	return (1, out) 
  
  
def transform(seq):
	"""
		Transforms "U" to "T" for the processing is done on DNA alphabet
	"""
	S = ""
	for s in seq:
		if s == "T":
			S += "U"
		else:
			S += s
	return S
  
  
def checkSimilarLength(s, SC):
	"""
		Compares sequence and structure constraint length
	"""
	if len(s) == len(SC):
		return 1
	else:
		return 0
  
  
def isStructure(s):
	"""
		Checks if the structure constraint only contains "(", ")", and "." and legal fuzzy structure constraint characters.
	"""
	returnvalue = 1
	for a in range(0,len(s)):
		if s[a] not in  ".()[]{}<>":
			if s[a] not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
				returnvalue = 0
	return returnvalue   


def isBalanced(s):
	"""
		Check if the structure s is of a balanced nature
	"""
	
	balance = 1
	for bracket in ["()", "[]", "{}", "<>"]:
		counter = 0
		for a in xrange(len(s)):
			if s[a] in bracket[0]:
				counter += 1
			elif s[a] in bracket[1]:
				counter -= 1
		if counter != 0:
			balance = 0
	return balance



def fulfillsHairpinRule(s):
	"""
		CHECKING FOR THE 3 nt LOOP INTERSPACE
			for all kind of basepairs, even wihtin the pdeudoknots 
	"""
	
	fulfillsRules = 1
	for bracket in ["()", "[]", "{}", "<>"]:
		last_opening_char = 0
		check = 0
		for a in xrange(len(s)):
			if s[a] == bracket[0]:
				last_opening_char = a
				check = 1
			elif s[a] == bracket[1] and check == 1:
				check = 0
				if a - last_opening_char < 4:
					return 0
	return 1
    
    
def isValidStructure(s):
	"""
		Checks, if the structure s is a valid structure
	"""
	
	Structure = isStructure(s)
	Balanced = isBalanced(s)
	HairpinRule = fulfillsHairpinRule(s)
	
	if Structure == 1 and Balanced == 1 and HairpinRule == 1:
		return 1
	else:
		print Structure, Balanced, HairpinRule
		return 0

def loadIUPACcompatibilities(IUPAC, useGU):
	"""
		Generating a hash containing all compatibilities of all IUPAC RNA NUCLEOTIDES
	"""
	compatible = {}
	for nuc1 in IUPAC: # ITERATING OVER THE DIFFERENT GROUPS OF IUPAC CODE
		sn1 = list(IUPAC[nuc1])
		for nuc2 in IUPAC: # ITERATING OVER THE DIFFERENT GROUPS OF IUPAC CODE
			sn2 = list(IUPAC[nuc2])
			compatib = 0
			for c1 in sn1: # ITERATING OVER THE SINGLE NUCLEOTIDES WITHIN THE RESPECTIVE IUPAC CODE:
				for c2 in sn2: # ITERATING OVER THE SINGLE NUCLEOTIDES WITHIN THE RESPECTIVE IUPAC CODE:
					# CHECKING THEIR COMPATIBILITY
					if useGU == True:
						if (c1 == "A" and c2 == "U") or (c1 == "U" and c2 == "A") or (c1 == "C" and c2 == "G") or (c1 == "G" and c2 == "C") or (c1 == "G" and c2 == "U") or (c1 == "U" and c2 == "G"): 
							compatib = 1
					else:
						if (c1 == "A" and c2 == "U") or (c1 == "U" and c2 == "A") or (c1 == "C" and c2 == "G") or (c1 == "G" and c2 == "C"): 
							compatib = 1
			compatible[nuc1 + "_" + nuc2] = compatib # SAVING THE RESPECTIVE GROUP COMPATIBILITY, REVERSE SAVING IS NOT REQUIRED, SINCE ITERATING OVER ALL AGAINST ALL
	return compatible
	
def isCompatibleToSet(c1, c2, IUPAC_compatibles):
	"""
		Checks compatibility of c1 wihtin c2
	"""
	compatible = True
	for setmember in c2:
		#print setmember
		if isCompatible(c1, setmember, IUPAC_compatibles) == False: 
			return False
	return compatible
	
	
def isCompatible(c1, c2, IUPAC_compatibles):
	"""
		Checks compatibility between character c1 and c2
	"""
	if IUPAC_compatibles[c1 + "_" + c2] == 1:
		return True
	else:
		return False
  
  
def isStructureCompatible(lp1, lp2 ,bp): 
	"""
		Checks, if the region within lp1 and lp2 is structurally balanced
	"""
	x = lp1 + 1
	while (x < lp2):
		if (bp[x] <= lp1 or bp[x] > lp2):
			return False
		if x == bp[x]:
			x += 1
		else:
			x = bp[x] + 1
	return x == lp2
 
 
def checkConstaintCompatibility(basepairstack, sequenceconstraint, IUPAC_compatibles):
	"""
		Checks if the constraints are compatible to each other
	"""
	returnstring = ""
	compatible = 1
	for id1 in basepairstack:  # key = (constraint , (pos, constraint)))
		constr1 = basepairstack[id1][0]
		id2 = basepairstack[id1][1][0]
		constr2 = basepairstack[id1][1][1]
    
		if id1 != id2 and not isCompatible(constr1, constr2, IUPAC_compatibles):
			
			compatible = 0
			returnstring += "nucleotide constraint " + str(constr1) + " at position " + str(id1) + " is not compatible with nucleotide constraint " + str(constr2) + " at position " + str(id2) + "\n"
    #if not isCompatible(basepairstack[basepair][0], basepairstack[basepair][1][1]):
      
      #compatible = 0
    #else:
      #returnstring += "nucleotide constraint " + str(basepairstack[basepair][0]) + " at position " + str(basepair) + " is compatible with nucleotide constraint " + str(basepairstack[basepair][1][1]) + " at position " + str(basepairstack[basepair][1][0]) + "\n"
	return (compatible, returnstring)

    
def getLP(BPSTACK):  
  """
    Retreives valid lonley base pairs from a base pair stack
  """
  #20 ('N', (>BLOCK<, 'N'))

  # geting single base pairs
  stack = {}
  LP = {}
  if type(BPSTACK[random.choice(BPSTACK.keys())]) == types.TupleType:
    for i in BPSTACK.keys():
      #if str(BPSTACK[i][1][0]) not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      stack[i] = int(BPSTACK[i][1][0])
	#print i , BPSTACK[i][1][0]
  else: 
    for i in BPSTACK.keys():
      #if str(BPSTACK[i]) not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      stack[i] = BPSTACK[i]

  # removing redundant base pair indices 
  for i in stack.keys():
    if i >= stack[i]:
	del stack[i]

  # actual checking for single lonley base pairs
  for i in stack.keys():
    if not (i-1 in stack and stack[i-1] == stack[i] + 1) and not (i+1 in stack and stack[i+1] == stack[i] - 1): 
      LP[i] = stack[i]

  ##actual removal of 2er lonley base pairs
  for i in stack.keys():
    if not (i-1 in stack and stack[i-1] == stack[i] + 1) and  (i+1 in stack and stack[i+1] == stack[i] - 1) and not (i+2 in stack and stack[i+2] == stack[i] - 2): 
      LP[i] = stack[i]
      LP[i+1] = stack[i+1]
  
  
  #if type(BPSTACK[random.choice(BPSTACK.keys())]) == types.TupleType:
    #for i in BPSTACK.keys():
      
      ##if str(BPSTACK[i][1][0]) not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      #stack[i] = int(BPSTACK[i][1][0])
	##print i , BPSTACK[i][1][0]
  #else:
    #for i in BPSTACK.keys():
      ##if str(BPSTACK[i]) not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      #stack[i] = BPSTACK[i]

  #for i in stack.keys():
    #if i >= stack[i]:
	#del stack[i]

  
  
  return LP
  
  
def getBPStack(s, seq):
	"""
		Returns a dictionary of the corresponding basepairs of the structure s and the sequence constraint seq.
	"""
	tmp_stack = {"()":[], "{}":[], "[]":[], "<>":[]}
	bpstack = {}
	for i in xrange(len(s)):
		
    # REGULAR SECONDARY STRUCTURE DETECTION
		if s[i] in "(){}[]<>":

			no = 0
			### opening
			if s[i] in "([{<":
				if s[i] == "(":
					tmp_stack["()"].append((i, seq[i]))
				elif s[i] == "[":
					tmp_stack["[]"].append((i, seq[i]))
				elif s[i] == "{":
					tmp_stack["{}"].append((i, seq[i]))
				elif s[i] == "<":
					tmp_stack["<>"].append((i, seq[i]))

			#closing
			elif s[i] in ")]}>":
				if s[i]  == ")":
					no, constr = tmp_stack["()"].pop() 
				elif s[i]  == "]":
					no, constr = tmp_stack["[]"].pop() 
				elif s[i]  == "}":
					no, constr = tmp_stack["{}"].pop() 
				elif s[i]  == ">":
					no, constr = tmp_stack["<>"].pop() 
				bpstack[no] = (constr, (i, seq[i])) 
				bpstack[i] = (seq[i] ,(no, constr)) 

		elif s[i] == ".": 
			bpstack[i] = (seq[i], (i, seq[i])) 
		elif s[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
			bpstack[i] = (seq[i], (i, seq[i])) 

	return (bpstack, getLP(bpstack))


def getbpStack(s): 
	"""
		Returns a dictionary of the corresponding basepairs of the structure s and the sequence constraint seq.
	"""
	tmp_stack = {"()":[], "{}":[], "[]":[], "<>":[]}
	bpstack = {}

	for i in xrange(len(s)):
		if s[i] in "(){}[]<>":

			no = 0
			### opening
			if s[i] in "([{<":
				if s[i] == "(":
					tmp_stack["()"].append(i)
				elif s[i] == "[":
					tmp_stack["[]"].append(i)
				elif s[i] == "{":
					tmp_stack["{}"].append(i)
				elif s[i] == "<":
					tmp_stack["<>"].append(i)

			#closing
			elif s[i] in ")]}>":
				if s[i] == ")":
					no = tmp_stack["()"].pop() 
				elif s[i] == "]":
					no = tmp_stack["[]"].pop() 
				elif s[i] == "}": 
					no = tmp_stack["{}"].pop() 
				elif s[i] == ">": 
					no = tmp_stack["<>"].pop() 
				bpstack[no] = i # save basepair in the format {opening base id (opening seq constr,(closing base id, closing seq constr))}
				bpstack[i] = no # save basepair in the format {closing base id (closing seq constr,(opening base id, opening seq constr))}

		elif s[i] == ".": # no structural constaint given: produce entry, which references itself as a base pair partner....
			bpstack[i] = i
	
		elif s[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
			bpstack[i] = i

    #elif s[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      ## per position, assigned to a certain block, the target nucleotide, with whcih it should interact is marked with the specific
      ## block character
      #bpstack[i] = s[i]
    
	return (bpstack, getLP(bpstack))

def maprange( a, b, s):
  """
   Mapping function
  """
  (a1, a2), (b1, b2) = a, b
  return  b1 + ((s - a1) * (b2 - b1) / (a2 - a1))


def applyGCcontributionPathAdjustment(pathlength, tmpGC, nt):
	"""
		GC path length contribution calculation.
	"""
	GCadjustment = 1.5
	minimum = 0.5
	upper = GCadjustment
	lower = minimum

	if nt == "A" or nt == "U":
		pathlength = pathlength * maprange( (0, 1) , (lower, upper), tmpGC)
    
	if nt == "G" or nt == "C":
		#pathlength = pathlength * (float(1-tmpGC))
		pathlength = pathlength * maprange( (1, 0) , (lower, upper), tmpGC)
	return pathlength

def getConstraint(TE, BPstack, IUPAC, IUPAC_compatibles, IUPAC_reverseComplements):
	"""
	Dependend on the situation in the constraint an the respective path section, setting wether a specific constraint can be given or not (for that path section)
	"""
	# TE :: transition element / path section under dispute
	# id1 :: id of the position of the caharacter to which the transition is leading to
	# id2 :: id of the position of the character, which is listed in the BPinformation, it can be id1 as well, when no bp is present
	# val :: BPstack information of the specific position
	# constr1 :: constraining character of pos id1
	# constr2 :: constraining character of pos id2

	id1 = int(TE.split(".")[0])
	val = BPstack[id1] # check out the value of the destination character in the basepair/constraint stack
	constr1 = val[0] # getting the constraint character of position id1
	id2 = int(val[1][0]) # getting position id2
	constr2 = val[1][1] # getting the sequence constraint for position id2
	targetNucleotide = TE.split(".")[1][-1:] # where the edge is leading to
	
	c1 = set(IUPAC[constr1]) # getting all explicit symbols of c1
	c2 = set(IUPAC_reverseComplements[constr2]) # getting the reverse complement explicit symbols of c2

	if targetNucleotide in c1:
		if id1 == id2:
			return 1
		else:
			if targetNucleotide in c2:
				return 1
			else:
				return 0
	else:
		return 0

"""
def getConstraint(TE, BPstack):
  # TE :: transition element / path section
  # id1 :: id of the position of the caharacter to which the transition is leading to
  # id2 :: id of the position of the character, which is listed in the BPinformation, it can be id1 as well, when no bp is present
  # val :: BPstack information of the specific position
  # constr1 :: constraining character of pos id1
  # constr2 :: constraining character of pos id2
  
  ### BPstack [id1] = (constr1, (id2, constr2))
  
  id1 = TE.split(".")[0]
  #print id1
  #id1 = TE.find(TE.strip("_")) #  strip the path section and getting the position of the section 
  #if len(TE.strip("_")) == 2: # check if the path section is from an internal and not an initial transition
    #id1 += 1 # increase position id1 by 1, since the last character of the section is the destination character
  val = BPstack[int(id1)] # check out the value of the destination character in the basepair/constraint stack
  constr1 = val[0] # getting the constraint character of position id1
  id2 = val[1][0] # getting position id2
  constr2 = val[1][1] # getting the sequence constraint for position id2
  #print TE, id1, constr1, id2, constr2,
  
  #TE.split(".")[1][-1:]
  if id1 == id2: # both ids were the same with either character, sequential or no sequential constraint -> no basepair constraint
    if constr1 == TE.split(".")[1][-1:] and constr2 == TE.split(".")[1][-1:]: # case if the single base constraints on position id1 == id2 are the same as the destination character on id1
      #print 1
      return 1
    elif constr1 == constr2 == "N": # case if the single base constraints on position id1 == id2 has no constraint 
      #print 1
      return 1
    else: # single base sequence constraints differ
      #print 0
      return 0
    
  elif id1 != id2: # showing differentq ids, indicating a bp, (basepair structural constraint)
    if constr1 == "N" and constr2 == "N": # no sequence constraint
      #print 1
      return 1
    if constr1 == "N" and constr2 != "N": # c1 has no constraint, c2 has character constraint (sequence constraint of closing bases)
      if TE.split(".")[1][-1:] == complementBase(constr2): # the current path section destination base is equal to the complement base of the mentioned sequence constraint in constr2
	#print 1
	return 1
      else: # case if the current path section destination base is not equeal to the mentioned complement sequence constraint in constr2
	#print 0
	return 0
    if constr1 != "N" and constr2 == "N": # c1 has character constraint, c2 has no character constraint (sequence constraint in the opening base)
      if TE.split(".")[1][-1:] == constr1:  # the current path section destination base is as constrained with constr1
	#print 1
	return 1
      else: # the current path section destination base is not as constrained in constr1
	#print 0
	return 0
    if constr1 != "N" and constr2 != "N": # both positions have sequential constraint
      if TE.split(".")[1][-1:] == constr1:
	#print 1
	return 1
      else:
	#print 0
	return 0
"""

def applyTerrainModification(terrain, s, tmpGC, SC, BPstack, IUPAC, IUPAC_compatibles, IUPAC_reverseComplements):
	#nucleotides = {'A': 0, 'C': 1,'G': 2,'T': 3}
	
	dels = []
	for terrainelement in sorted(terrain):
		pheromone, pathlength = terrain[terrainelement]
		pheromone = getConstraint(terrainelement, BPstack, IUPAC, IUPAC_compatibles, IUPAC_reverseComplements)
		pathlength = getConstraint(terrainelement, BPstack, IUPAC, IUPAC_compatibles, IUPAC_reverseComplements)
		pathlength = applyGCcontributionPathAdjustment(pathlength, tmpGC,terrainelement.split(".")[1][-1:])
		if pheromone * pathlength == 0: dels.append(terrainelement)
		terrain[terrainelement] = (pheromone, pathlength,[])
	further_dels = {}
	for terrainelement in sorted(dels):
		pos, nucs = terrainelement.split(".")
		if int(pos) < len(s)-1:
			to_nt = nucs[-1:]
			successor_pos = int(pos) + 1
			for i in ["A", "C", "G", "U"]:
				del_element = str(successor_pos) + "." + to_nt + i
				further_dels[del_element] = 1
		further_dels[terrainelement] = 1
	# deleting the inbound and outbound edges, which are forbidden
	for terrainelement in further_dels:
		del terrain[terrainelement]
	# allocate the appropriate children of edges 
	for terrainelement in terrain:
		pheromone, pathlength, children = terrain[terrainelement]
		pos, nucs = terrainelement.split(".")
		if int(pos) < len(s):
			to_nt = nucs[-1:]
			successor_pos = int(pos) + 1
			for i in ["A", "C", "G", "U"]:
				if str(successor_pos) + "." + to_nt + i in terrain:
					children.append(str(successor_pos) + "." + to_nt + i)
		terrain[terrainelement] = (pheromone, pathlength,children)
	starts = []
	for i in ["A", "C", "G", "U"]:
		if str(0) + "." + i in terrain:
			starts.append(str(0) + "." + i)
	terrain["00.XY"] = (1, 1, starts)
	return (terrain, BPstack)


def initTerrain(s): # THE CLASSIC
	"""
		Initialization of the terrain with graph like terrain... vertices are modeled implicitly
	"""
	nt = ["A","C","G","U"] 
	nt2 = ["AA","AC","AG","AU","CA","CC","CG","CU","GA","GC","GG","GU","UA","UC","UG","UU"] # Allowed dinucleotides
	e = {}
	pathlength = 1
	pheromone = 1
	for p in xrange(len(s)):
		if p == 0:
			for i in nt:
				e["%s.%s"%(p,i)] = (pheromone, pathlength)
		elif p > 0:
			for n in nt2:
				e["%s.%s"%(p,n)] = (pheromone, pathlength)
	return e
  
  

def complementBase(c):
	"""
		Returns the complement RNA character of c (without GU base pairs)
	"""
	retChar = ""
	if c == "A" :
		retChar = "U"
	elif c == "U":
		retChar = "A"
	elif c == "C":
		retChar = "G"
	elif c == "G":
		retChar = "C"
	return retChar
  
def printTerrain(terrain):
	#print sorted(terrain.keys())
	tmp_i = "0"
	tmp_c = 0
	terrain = terrain[0]
	
	for a, i in enumerate(sorted(terrain.keys())):
		#print a
		if i.split(".")[0] != tmp_i:
			print "\nElements:", tmp_c,"\n#########################\n", i, terrain[i]
			
			tmp_c = 1
			tmp_i = i.split(".")[0]
		else:
			print i, terrain[i]
			tmp_c += 1
			
	print "\nElements:", tmp_c
	print "#########################"
	print len(terrain)
	
def pickStep(tmp_steps, summe):
	"""
		Selects a step within the terrain
	"""
	if len(tmp_steps) == 1:
		return tmp_steps[0][1] # returning the nucleotide of the only present step
	else:
		rand = random.random() # draw random number
		mainval = 0
		for choice in xrange(len(tmp_steps)):
			val, label = tmp_steps[choice]
			mainval += val/float(summe)
			if mainval > rand: # as soon, as the mainval gets larger than the random value the assignment is done
				return label

def getPath(s, tmp_terrain, tmp_BPstack, alpha, beta, IUPAC, IUPAC_reverseComplements):
	"""
		Performs a walk through the terrain and assembles a sequence, while respecting the structure constraint and IUPAC base complementarity
		of the base pairs GU, GC and AT
	"""
	nt = ["A","C","G","U"]
	prev_edge = "00.XY"
	sequence = ""
	while len(sequence) <  len(s):
		coming_from = sequence[-1:]
		summe = 0
		steps = []
		i = len(sequence)
		allowed_nt = "ACGU"
		# base pair closing case check, with subsequent delivery of a reduced allowed nt set
		
		if i > tmp_BPstack[i][1][0]:
			jump =  tmp_BPstack[i][1][0]
			nuc_at_jump = sequence[jump]
			allowed_nt = IUPAC_reverseComplements[nuc_at_jump]
			
			#allowed_nt = complementBase(nuc_at_jump)
			
		# Checking for every possible nt if it is suitable for the selection procedure
		for edge in tmp_terrain[prev_edge][-1]:
			
			if edge[-1:] in allowed_nt:
				pheromone, PL , children = tmp_terrain[edge]
				#if PL > 0:
				value = ((float(pheromone * alpha)) + ((1/float(PL)) * beta))
				summe += value
				steps.append((value, edge))
		prev_edge = pickStep(steps, summe)
		sequence += prev_edge[-1:]
		
	return sequence


###
# STRUCTURE PREDICTORS
###
def getPKStructure(sequence, temperature, mode = "A"):
	"""
		Initialization pKiss mfe pseudoknot prediction
	"""
	p2p = "pKiss"
	#p2p = "/usr/local/pkiss/2014-03-17/bin/pKiss_mfe"
	strategy = "--strategy "
	t = "--temperature " + str(temperature)
	
	if mode == "A": strategy += "A"
	elif mode == "B": strategy += "B"
	elif mode == "C": strategy += "C"
	elif mode == "D": strategy += "D"
	elif mode == "P": strategy += "P"
	
	p = subprocess.Popen( ([p2p, "--mode mfe", strategy, t]),
				#shell = True,
				stdin = subprocess.PIPE,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				close_fds = True)
	#print p.stderr.readline()

	p.stdin.write(sequence+'\n')
	pks = p.communicate()
	structure = "".join(pks[0].split("\n")[2].split(" ")[-1:])
	return structure
	
def init_RNAfold(version, temperature, paramFile = ""):
	"""
		Initialization RNAfold listener
	"""
	p2p = ""
	t = "-T " + str(temperature)
	P = ""
	if paramFile != "":
		P = "-P " + paramFile
	if version == 185:
		p2p = "/home/rk/Software/ViennaRNA/ViennaRNA-1.8.5/Progs/RNAfold"
		p = subprocess.Popen( ([p2p, '--noPS', '-d 2', t, P]),
					shell = True,
					stdin = subprocess.PIPE,
					stdout = subprocess.PIPE,
					stderr = subprocess.PIPE,
					close_fds = True)
		return p
	elif version == 213:
		p2p = "RNAfold"
		p = subprocess.Popen( ([p2p, '--noPS', '-d 2', t, P]),
					#shell = True,
					stdin = subprocess.PIPE,
					stdout = subprocess.PIPE,
					stderr = subprocess.PIPE,
					close_fds = True)
		return p
	else:
		exit(0)
		
def consult_RNAfold(seq, p):
	"""
		Consults RNAfold listener
	"""
	p.stdin.write(seq+'\n')
	out = ""
	for i in xrange(2):
		out += p.stdout.readline()
	return out


def getRNAfoldStructure(struct2, process1):
	"""
		Retrieves folded structure of a RNAfold call
	"""
  
	RNAfold_pattern = re.compile('.+\n([.()]+)\s.+')
	#RNAdist_pattern = re.compile('.*\s([\d]+)')
	RNAfold_match = RNAfold_pattern.match(consult_RNAfold(struct2, process1))
	current_structure = ""
	#if RNAfold_match:
	return RNAfold_match.group(1)
  
  
def init_RNAdistance():
	"""
		Initialization of RNAdistance listener
	"""
	#p2p = "/home/rk/Software/ViennaRNA/ViennaRNA-1.8.5/Progs/RNAdistance"
	p2p = "RNAdistance"
	p = subprocess.Popen( ([p2p]),
							#shell = True,
							stdin = subprocess.PIPE,
							stdout = subprocess.PIPE,
							stderr = subprocess.PIPE,
							close_fds = True)
	return p


def consult_RNAdistance(s1, s2, p):
	"""
		Consulting the RNAdistance listener
	"""
	p.stdin.write(s1+'\n')
	p.stdin.write(s2+'\n')
	out = ""
	out_tmp = p.stdout.readline().strip()
	if out_tmp != "":
		out += out_tmp
	return out

def getInducingSequencePositions(Cseq, degreeOfSequenceInducement):
	"""
		Delimiting the degree of structure inducement by the supplied sequence constraint.
		0 : no sequence induced structure constraint
		1 : "ACGT" induce structure (explicit nucleotide structure inducement level)
		2 : "MWKSYR" and "ACGT" (explicit and double instances)
		3 : "BDHV" , "MWKSYR" and "ACGT" (explicit, double, and triple instances)
	"""
	setOfNucleotides = "" # resembling the "0"-case
	if degreeOfSequenceInducement == 1:
		setOfNucleotides = "ACGU"
	elif degreeOfSequenceInducement == 2:
		setOfNucleotides = "ACGUMWKSYR"
	elif degreeOfSequenceInducement == 3:
		setOfNucleotides = "ACGUMWKSYRBDHV"
	#elif degreeOfSequenceInducement == 4:
		#setOfNucleotides = "ACGTMWKSYRBDHVN"
		
	tmpSeq = ""
	listset = setOfNucleotides
	for pos in Cseq:
		if pos not in listset:
			tmpSeq += "N"
		else:
			tmpSeq += pos
	
	return setOfNucleotides, tmpSeq

	
def getBPDifferenceDistance(stack1, stack2):
	"""
		Based on the not identical amount of base pairs within both structure stacks
	"""
	d = 0
	for i in stack1.keys():
		# check base pairs in stack 1
		if i < stack1[i] and stack1[i] != stack2[i]:
			d += 1
		# check base pairs in stack 2
	for i in stack2.keys():
		if i < stack2[i] and stack1[i] != stack2[i]:
			d += 1
	return d


def getStructuralDistance(target_structure, Cseq,  path, RNAfold, verbose, LP, BP, RNAfold_pattern, IUPAC_compatibles, degreeOfSequenceInducement, pseudoknots, strategy):
	"""
		Calculator for Structural Distance
	"""
	# fold the current solution's sequence to obtain the structure
	
	current_structure = ""
	
	if pseudoknots:
		current_structure = getPKStructure(path,strategy)
	else:
		RNAfold_match = RNAfold_pattern.match(consult_RNAfold(path, RNAfold))
		current_structure = RNAfold_match.group(1)

	# generate the current structure's base-pair stack
	bp = getbpStack(current_structure)[0]
	# add case-dependend structural constraints in case of lonley basepairs formation
	tmp_target_structure_bp = getbpStack(target_structure)[0]

	for lp in LP:
		if bp[lp] == LP[lp]: # if the base pair is within the current solution structure, re-add the basepair into the constraint structure.
			#tmp_target_structure[lp] = "("
			#tmp_target_structure[LP[lp]] = ")"
			tmp_target_structure_bp[lp] = LP[lp]
			tmp_target_structure_bp[LP[lp]] = lp
			
	# REMOVE BLOCK CONSTRAINT AND SUBSTITUTE IT WITH SINGLE STRAND INFORMATION repsective with brackets, if allowed base pairs occure
	# check for all allowed implicit constraint block declarators
	for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
		occurances = []
		for m in re.finditer(c, target_structure): # search for a declarator in the requested structure
			occurances.append(m.start()) # save the corresponding index

		# transform declarator into single stranded request
		for i in occurances: 
			#tmp_target_structure[i] = "."
			tmp_target_structure_bp[i] = i
		# infer a base pair within the block declarated positions, if the current structure provides it.
		for i in occurances:
			for j in occurances:
				if i < j:
					if bp[i] == j:
						#tmp_target_structure[i] = "("
						#tmp_target_structure[bp[i]] = ")"
						
						tmp_target_structure_bp[i] = bp[i]
						tmp_target_structure_bp[bp[i]] = i
						
	# CHECK FOR SEQUENCE CONSTRAINT WHICH INDUCES STRUCTURE CONSTRAINT IN THE MOMENTARY SITUATION
	#print "Checking Cseq influence and it's induced basepairs..."
	IUPACinducers, tmp_Cseq = getInducingSequencePositions(Cseq, degreeOfSequenceInducement)
	if len(Cseq.strip("N")) > 0:
    #print "Processing Cseq influence"
		# Iterate over all positions within the Base Pair stack
		for i in BP: # Check for each base index i 
			
			if i < bp[i]: # if the current index is samller that the affiliated in the basepair stack of the current solution

				bp_j = bp[i] # Actual j index of the current solution
				BP_j = BP[i][1][0] # j index of the requested structure
				if (i != bp_j and i == BP_j and BP[i][0] in IUPACinducers ): # if i pairs with some other base in the current structure, and i is requested single stranded and the Sequence constraint is allowed to induce...
					if (BP[bp_j][1][0] == bp_j and BP[bp_j][0] in IUPACinducers):# If position j is requested singlestranded and position j nucleotide can induce base pairs
						#if isCompatible(bp[i][0], bp[i][1][1], IUPAC_compatibles): # If both nucleotides, i and j are actually compatible
						#tmp_target_structure[i] = "("
						#tmp_target_structure[bp_j] = ")"
						
						tmp_target_structure_bp[i] = bp[i]
						tmp_target_structure_bp[bp_j] = i
						
	#tts = "".join(tmp_target_structure)
	dsreg = getBPDifferenceDistance(tmp_target_structure_bp, bp)
		
	# CHECK FOR ALL DETERMINED LONELY BASE PAIRS (i<j), if they are formed
	failLP = 0
	for lp in LP: 

		if bp[lp] != LP[lp]: 
	
			isComp = isCompatible(path[lp],path[LP[lp]], IUPAC_compatibles)
			isStru = isStructureCompatible(lp, LP[lp] ,bp)
			if not ( isStru and isStru ): # check if the bases at the specific positions are compatible and check if the 
		# basepair can be formed according to pseudoknot free restriction. If one fails, a penalty distance is raised for that base pair
				failLP += 1
				
	#print dsreg, failLP, float(len(tmp_target_structure_bp))
	dsLP = float(failLP)
  
	return (dsreg + dsLP) /float(len(tmp_target_structure_bp))
  
  
def getGC(sequence):
	"""
		Calculate GC content of a sequence
	"""
	GC = 0
	for nt in sequence:
		if nt == "G" or nt == "C":
			GC = GC + 1
	GC = GC/float(len(sequence))
	return GC
  

def getGCDistance(tGC, gc2, L):
  """
    Calculate the pseudo GC content distance 
  """
  nt_coeff = L * tGC
  pc_nt = (1/float(L))*100
  #     
  d = gc2 - tGC
  d = d * 100
  
  f = math.floor(nt_coeff)
  c = math.ceil(nt_coeff)

  if d < 0: # 
    #print "case x",(abs(nt_coeff - f)), pc_nt, (abs(nt_coeff - f)) * pc_nt,
    d = d + (abs(nt_coeff - f)) * pc_nt
  elif d > 0: # case y
    #print "case y", abs(nt_coeff - c), pc_nt, abs(nt_coeff - c) * pc_nt,
    d = d - abs(nt_coeff - c) * pc_nt
  elif d == 0:
    pass
  
  d = round(d, 7)
  #d = max(0, abs(d)- ( max ( abs( math.ceil(nt_coeff)-(nt_coeff)) , abs(math.floor(nt_coeff)-(nt_coeff)) )/L)*100 )  
  return abs(d)


def getSequenceEditDistance(SC, path):
  """
    Calculate sequence edit distance of a solution to the constraint
  """#
  IUPAC = {"A":"A", "C":"C", "G":"G", "U":"U", "R":"AG", "Y":"CU", "S":"GC", "W":"AU","K":"GU", "M":"AC", "B":"CGU", "D":"AGU", "H":"ACU", "V":"ACG", "N":"ACGU"}         
  edit = 0
  for i in xrange(len(SC)):
    if path[i] not in IUPAC[SC[i]]:
      edit += 1
  return edit/float(len(path))



def getTransitions(p):
	"""
		Retreive transitions of a specific path/sequence
	"""
	transitions = []
	for pos in xrange(len(p)):
		if pos == 0:
			transitions.append(str(pos) + "." + p[pos])

		else:
			insert = p[pos-1] + p[pos]
			transitions.append(str(pos) + "." + insert)

	return transitions

  
def evaporate(t, er): 
	"""
	Evaporate the terrain's pheromone trails
	"""
	terr, BP = t
	c = 1
	for key in terr:
		p,l,c = terr[key]
		p *= (1-er)
		terr[key] = (p, l, c)


def updateValue(distance, correction_term, omega):
  """
    Retrieves a distance dependend pheromone value
  """
  if correction_term == 0:
    return 0
  else:
    if distance == 0:
      return omega * correction_term
    else:
      return (1/float(distance)) * correction_term
 
 
def trailBlaze(p, c_s, s, ds, dgc, dseq, dn, t, correction_terms, BPstack, verbose):
	"""
		Pheromone Update function accorinding to the quality of the solution
	"""
	terr, BP = t
	bpstack, LP = getbpStack(c_s)

	struct_correction_term , GC_correction_term, seq_correction_term = correction_terms
	omega = 2.23

	bs = updateValue(ds, struct_correction_term, omega)
	bGC = updateValue(dgc, GC_correction_term, omega)
	if dseq != "n.a.":
		bSeq = updateValue(dseq, seq_correction_term, omega)
		d = bs + bGC + bSeq
	else:
		d = bs + bGC  
	transitions = getTransitions(p)

	for trans in xrange(len(transitions)): # for each transition in the path
		id1 = int(transitions[trans].split(".")[0])
		tar_id2 = int(BPstack[id1][1][0]) # getting requested  position id2
		curr_id2 = int(bpstack[id1]) # getting the current situation
		multiplicator = 0
		if tar_id2 == curr_id2 and id1 != tar_id2 and id1 != curr_id2: # case of a base pair, having both brackets on the correct position
			multiplicator = 1
		elif tar_id2 == curr_id2 and id1 == tar_id2 and id1 == curr_id2: # case of a single stranded base in both structures
			multiplicator = 1
		p, l, c = terr[transitions[trans]] # getting the pheromone and the length value of the single path transition
		p +=  d * multiplicator
		terr[transitions[trans]] = (p, l, c) # updating the values wihtin the terrain's
	t = (terr, BP)


def updateTerrain(p, c_s, s, ds, dgc, dseq, dn, t, er, correction_terms, BPstack, verbose, ant_count):
	"""
		General updating function
	"""
	evaporate(t,er)
	trailBlaze(p, c_s, s, ds, dgc, dseq, dn, t, correction_terms, BPstack, verbose)

  
def getUsedTime(start_time):
	"""
		Return the used time between -start time- and now.
	"""
	end_time = time.time()
	return end_time - start_time


def good2Go(SC, L, CC, STR):
	"""
		Check, if all input is correct and runnable
	"""
	if (SC == 1 and L == 1 and CC == 1 and STR == 1):
		return True
	else:
		print SC,L,CC,STR
		return False
  
  
def getPathFromSelection( aps, s, terrain, alpha, beta,  RNAfold, RNAfold_pattern, GC, SC, LP, verbose, IUPAC_compatibles, degreeOfSequenceInducement, IUPAC_reverseComplements, IUPAC, pseudoknots, strategy):
	"""
		Returns the winning path from a selection of pathes...
	"""
	terr, BPs = terrain
	win_path = 0
	for i in xrange(aps):
		# Generate Sequence
		path = getPath(s, terr, BPs, alpha, beta, IUPAC, IUPAC_reverseComplements)
		# Measure sequence features and transform them into singular distances
		distance_structural = float(getStructuralDistance(s, SC , path, RNAfold, verbose, LP, BPs, RNAfold_pattern, IUPAC_compatibles, degreeOfSequenceInducement, pseudoknots, strategy)) 
		distance_GC = float(getGCDistance(GC,getGC(path), len(path)))
		distance_seq = float(getSequenceEditDistance(SC, path))
		# Calculate Distance Score
		D = distance_structural + distance_GC + distance_seq
      
		# SELECT THE BEST-OUT-OF-k-SOLUTIONS according to distance score
		if i == 0:
			win_path = (path, D, distance_structural, distance_GC, distance_seq)
		else:
			if D < win_path[1]:
				win_path = (path, D, distance_structural, distance_GC, distance_seq)
	return win_path


def substr(x, string, subst):
	"""
		Classical substring function
	"""
	s1 = string[:x-1]
  
	s2 = string[x-1:x]
	s3 = string[x:]
  #s2 = s[x+len(string)-x-1:]
  
	return s1 + subst + s3
  
  
def inConvergenceCorridor(d_struct, d_gc, BS_d_struct, BS_d_gc):
	"""
		Check if a solutions qualities are within the convergence corridor
	"""
	struct_var = ((BS_d_struct/float(4)) + 3 ) * 4
	gc_var = (BS_d_gc + 1/float(100) * 5) + BS_d_gc + 1

	if d_struct <= struct_var and d_gc <= gc_var:
		return True
	else:
		return False
  
def getGCSamplingValue(GC, tGCmax, tGCvar):
	"""
	Returns a suitable GC value, dependend on the user input: Either returning the single GC value,
	which the user entered, or a smpled GC value
	from a designated distribution in it's interavals
	"""
	returnval = 0
	if tGCmax == -1.0 and tGCvar == -1.0: # regular plain tGC value as requested 
		return GC
	elif tGCmax != -1.0 and tGCvar == -1.0: # uniform distribution tGC value sampling
		if GC < tGCmax:
			tmp_GC = tGCmax
			tGCmax = GC
			GC = tmp_GC
		while returnval <= 0:
			returnval = float(numpy.random.uniform(low=GC, high=tGCmax, size=1))
		return returnval
	elif tGCmax == -1.0 and tGCvar != -1.0: # normal distribution tGC value sampling
		while returnval <= 0:
			returnval = float(numpy.random.normal(GC, tGCvar, 1))
		return returnval
	
	
def reachableGC(C_struct):
	"""
		Checks if a demanded GC target content is reachable in dependence with the given sequence constraint.
	"""
	AU = 0
	for i in C_struct:
		if i == "A" or i == "U":
			AU += 1
	maxGC = 1 - (AU / float(len(C_struct))) # 1 - min_GC
	return maxGC
  
 
def runColony(s, SC, objective_to_target_distance, GC, alpha, beta, evaporation_rate, correction_terms, verbose, IUPAC, IUPAC_compatibles, degreeOfSequenceInducement, IUPAC_reverseComplements, termination_convergence, convergence_count, reset_limit, improve, temperature, paramFile, pseudoknots, strategy):
	"""
		Execution function of a single ant colony finding one solution sequence
	"""
	retString = ""
	retString2 = []
	BPstack, LP = getBPStack(s, SC)
   
	rGC = reachableGC(SC)
	GC_message = ""
	if GC > rGC:
		print >> sys.stderr, "WARNING: Chosen target GC %s content is not reachable due to sequence constraint! Sequence Constraint GC-content is: %s" % (GC, rGC) 
		GC = rGC
    
	# Initial Constraint Checks prior to execution
	STR = isValidStructure(s)
	START_SC , SC = checkSequenceConstraint(str(SC))
	START_LENGTH = checkSimilarLength(str(s), str(SC)) 
	START_constraint_compatibility , CompReport = checkConstaintCompatibility(BPstack, SC, IUPAC_compatibles)
  
	g2g = good2Go(START_SC, START_LENGTH, START_constraint_compatibility, STR)
	if (g2g == 1):
		start_time = time.time()
		max_time = 600 # seconds
		
		
		
			
		#### 
		# INITIALIZATION OF THE RNA TOOLs
		#
		RNAfold = init_RNAfold(213, temperature, paramFile)
		#RNAdistance = init_RNAdistance()
		RNAfold_pattern = re.compile('.+\n([.()]+)\s.+')
		#RNAdist_pattern = re.compile('.*\s([\d]+)')
		#
		####
		
		terrain = initTerrain(s) 
		#print len(terrain), 
		terrain = applyTerrainModification(terrain, s, GC, SC, BPstack, IUPAC, IUPAC_compatibles, IUPAC_reverseComplements)
		#print len(terrain[0])
		#printTerrain(terrain)
		#exit(0)
		global_ant_count = 0
		global_best_ants = 0
		criterion = False
		met = True  
		ant_no = 1
		prev_res = 0
		seq = ""

		counter = 0
		
		dstruct_log = []
		dGC_log = []
		
		
		distance_structural = 1000
		distance_GC = 1000
		distance_seq = 1000
		
		convergence = convergence_count
		convergence_counter = 0
		
		resets = 0
		
		path = ""
		curr_structure = ""

		Dscore = 100000
		distance_structural = 10000
		distance_GC = 10000
		distance_seq = 10000
		best_solution = (path, curr_structure, Dscore, distance_structural, distance_GC, distance_seq)
		best_solution_local = (path, curr_structure, Dscore, distance_structural, distance_GC, distance_seq)
		
		best_solution_since = 0
		
		ants_per_selection = 10
		if len(LP) > 0 :
			for lp in LP:
				s = substr(lp + 1, s, ".")
				s = substr(LP[lp] + 1, s, ".")
    
		init = 1
		while criterion != met and getUsedTime(start_time) < max_time:
			iteration_start = time.time()
			global_ant_count += 1
			global_best_ants += 1

			path_info = getPathFromSelection(ants_per_selection, s, terrain, alpha, beta, RNAfold, RNAfold_pattern, GC, SC, LP, verbose, IUPAC_compatibles, degreeOfSequenceInducement, IUPAC_reverseComplements, IUPAC, pseudoknots, strategy)

			distance_structural_prev = distance_structural
			distance_GC_prev = distance_GC
			distance_seq_prev = distance_seq

			path, Dscore , distance_structural, distance_GC, distance_seq = path_info
			curr_structure = ""
			if pseudoknots:
				curr_structure = getPKStructure(path, strategy)
			else:
				curr_structure = getRNAfoldStructure(path, RNAfold)
				
			curr_solution = (path,curr_structure, Dscore, distance_structural, distance_GC, distance_seq)
			# BEST SOLUTION PICKING
			if improve == "h": # hierarchical check
				# for the global best solution
				if distance_structural < best_solution[3] or (distance_structural == best_solution[3] and distance_GC < best_solution[4]):
					best_solution = curr_solution
					ant_no = 1
				# for the local (reset) best solution
				if distance_structural < best_solution_local[3] or (distance_structural == best_solution_local[3] and distance_GC < best_solution_local[4]):
					best_solution_local = curr_solution
					
			elif improve == "s": #score based check
				# store best global solution
				if Dscore < best_solution[2]:
					best_solution = curr_solution
					ant_no = 1
				# store best local solution for this reset
				if Dscore < best_solution_local[2]:
					best_solution_local = curr_solution

# OLD ' BEST SOLUTION ' PICKING
#			if Dscore < best_solution[2]:
#				best_solution = (path,curr_structure, Dscore, distance_structural, distance_GC, distance_seq)
#      
#			if Dscore < best_solution_local[2]:
#				best_solution_local = (path,curr_structure, Dscore, distance_structural, distance_GC, distance_seq)


			distance_DN = 0
      
			if verbose:
				print "SCORE " + str(Dscore) + " Resets " + str(resets) + " #Ant " + str(global_ant_count) + " out of " + str(ants_per_selection)  + " cc " + str(convergence_counter)

				print s, " <- target struct" 
				print best_solution[0] , " <- BS since ", str(best_solution_since), "Size of Terrrain:", len(terrain[0])
				print best_solution[1] , " <- BS Dscore " + str(best_solution[2]) + " ds " + str(best_solution[3]) + " dGC " + str(best_solution[4]) + " dseq " + str(best_solution[5])+ " LP " + str(len(LP)) + " <- best solution stats"
				print curr_structure, " <- CS"
				print path,
				print " <- CS", "Dscore", str(Dscore), "ds", distance_structural, "dGC", distance_GC, "GC", getGC(path)*100, "Dseq", distance_seq
	
			#### UPDATING THE TERRAIN ACCORDING TO THE QUALITY OF THE CURRENT BESTO-OUT-OF-k SOLUTION
			updateTerrain(path, curr_structure, s, distance_structural,distance_GC, distance_seq, distance_DN, terrain, evaporation_rate, correction_terms, BPstack, verbose, global_ant_count) 
			####
			if verbose: print "Used time for one iteration", time.time() - iteration_start
				
				
			# CONVERGENCE AND TERMINATION CRITERION MANAGEMENT
			#print distance_structural, distance_GC, best_solution_local[3], best_solution_local[4]
			if inConvergenceCorridor(curr_solution[3], curr_solution[4], best_solution_local[3], best_solution_local[4]):
				convergence_counter += 1
			if distance_structural_prev == distance_structural and distance_GC_prev == distance_GC:
				convergence_counter += 1
	
			if best_solution[3] == objective_to_target_distance:
				if best_solution[4] == 0.0:
					break
				ant_no = ant_no + 1
				convergence_counter -= 1
			else:
				ant_no = 1

	  
			if ant_no == termination_convergence or resets >= reset_limit or global_ant_count >= 100000 or best_solution_since == 5:
				break

			# RESET
			if ant_no < termination_convergence and convergence_counter >= convergence:
	
				terrain = initTerrain(s)
				terrain = applyTerrainModification(terrain, s, GC, SC, BPstack, IUPAC, IUPAC_compatibles, IUPAC_reverseComplements)
				criterion = False
				met = True  
				ant_no = 1
				prev_res = 0
				pre_path = "_" * len(s)
				path = ""
				curr_structure = ""
				counter = 0
				Dscore = 100000
				distance_structural = 1000
				distance_GC = 1000
				distance_seq = 1000
				best_solution_local = (path, curr_structure, Dscore, distance_structural, distance_GC, distance_seq)
				convergence = convergence_count
				convergence_counter = 0

				if resets == 0:
					sentinel_solution = best_solution
					best_solution_since += 1
				else:
					if best_solution[2] < sentinel_solution[2]:
						sentinel_solution = best_solution
						best_solution_since = 0
					else:
						best_solution_since += 1

				resets += 1

		duration  = getUsedTime(start_time)

		retString += "|Ants:" + str(global_ant_count)
		retString += "|Resets:" + str(resets) + "/" + str(reset_limit)
		retString += "|AntsTC:" + str(termination_convergence) 
		retString += "|CC:" + str(convergence_count) 
		retString += "|IP:" + str(improve) 
		retString += "|BSS:" + str(best_solution_since)
		#if GC_message != "":
		#	retString += GC_message + "\n"
      
		sequence = best_solution[0]
		struct = best_solution[1]

		retString += "|LP:" + str(len(LP))
		retString += "|ds:" + str(getStructuralDistance(s,SC, sequence, RNAfold, verbose, LP, BPstack, RNAfold_pattern, IUPAC_compatibles, degreeOfSequenceInducement, pseudoknots, strategy))
		retString += "|dGC:" + str(best_solution[4])
		retString += "|GC:" + str(getGC(sequence)*100)
		retString += "|dseq:" + str(getSequenceEditDistance(SC, sequence))
		retString += "|L:" + str(len(sequence))
		retString += "|Time:" + str(duration)
		
		retString2.append(struct)
		retString2.append(sequence)
    
		# CLOSING THE PIPES TO THE PROGRAMS
		RNAfold.communicate()
		#RNAdistance.communicate()

	else: # Structural premisses are not met, htherefore the program will halt with a failure message
		retString +=  "\nSome mistake detected\n"
		retString +=  "SequenceConstraintCheck: " + str(START_SC) + "\nSequenceConstraint: " +  str(SC) + "\nLengthCheck: " + str(START_LENGTH) + "\nConstraintCompatibility: " + str(START_constraint_compatibility)+ "\n" + CompReport + "\n"

	return (retString, retString2) 

def findSequence(structure, Cseq, tGC, colonies, name, alpha, beta, evaporation_rate, struct_correction_term, GC_correction_term, seq_correction_term, degreeOfSequenceInducement, file_id, verbose, output_verbose, tGCmax, tGCvar, termination_convergence, convergence_count, reset_limit, improve, seed, temperature, paramFile, pseudoknots, strategy, useGU, return_mod = False):
	"""
		MAIN antaRNA - ant assembled RNA
	"""

	if seed != "none":
		random.seed(seed)
	
	if Cseq == "":
		sequenceconstraint = "N" * len(structure)
	else:
		sequenceconstraint = str(Cseq)
  
	alpha = float(alpha)
	beta = float(beta)
	tGC = float(tGC)
	evaporation_rate = float(evaporation_rate)
	struct_correction_term = float(struct_correction_term)
	GC_correction_term = float(GC_correction_term)
	seq_correction_term = float(seq_correction_term)
	colonies = int(colonies)
	file_id = str(file_id)
	tmp_verbose = verbose
	tmp_output_verbose = output_verbose
	verbose = tmp_output_verbose # Due to later change, this is a twistaround and a switching of purpose
	output_verbose = tmp_verbose # Due to later change, this is a twistaround and a switching of purpose
	correction_terms = struct_correction_term, GC_correction_term, seq_correction_term
	temperature = float(temperature)
	print_to_STDOUT = (file_id == "STDOUT")
	
	useGU = useGU
	
	if return_mod == False:
		if print_to_STDOUT == False:
			outfolder = '/'.join(file_id.strip().split("/")[:-1])
			curr_dir = os.getcwd()
			if not os.path.exists(outfolder):
				os.makedirs(outfolder)
			os.chdir(outfolder)  
		

	sequenceconstraint = transform(sequenceconstraint)
	###############################################
  
	# Allowed deviation from the structural target:
	objective_to_target_distance = 0.0

	# Loading the IUPAC copatibilities of nuleotides and their abstract representing symbols
	IUPAC = {"A":"A", "C":"C", "G":"G", "U":"U", "R":"AG", "Y":"CU", "S":"GC", "W":"AU","K":"GU", "M":"AC", "B":"CGU", "D":"AGU", "H":"ACU", "V":"ACG", "N":"ACGU"}         
	IUPAC_compatibles = loadIUPACcompatibilities(IUPAC, useGU)

	IUPAC_reverseComplements = {}
	if useGU == False: ## Without the GU basepair
		IUPAC_reverseComplements = {"A":"U", "C":"G", "G":"C", "U":"A", "R":"UC", "Y":"AG", "S":"GC", "W":"UA","K":"CA", "M":"UG", "B":"AGC", "D":"ACU", "H":"UGA", "V":"UGC", "N":"ACGU"}         
	else: ## allowing the GU basepair
		IUPAC_reverseComplements = {"A":"U", "C":"G", "G":"UC", "U":"AG", "R":"UC", "Y":"AG", "S":"UGC", "W":"UAG","K":"UCAG", "M":"UG", "B":"AGCU", "D":"AGCU", "H":"UGA", "V":"UGC", "N":"ACGU"}         
	
	result = []
	for col in xrange(colonies):
		# Checking the kind of taget GC value should be used
		GC = getGCSamplingValue(tGC, tGCmax, tGCvar)

		# Actual execution of a ant colony procesdure
		output_v, output_w  =  runColony(structure, sequenceconstraint, objective_to_target_distance, GC, alpha, beta, evaporation_rate, correction_terms, verbose, IUPAC, IUPAC_compatibles, degreeOfSequenceInducement, IUPAC_reverseComplements, termination_convergence, convergence_count, reset_limit, improve, temperature, paramFile, pseudoknots, strategy)

		# Post-Processing the output of a ant colony procedure
		line = ">" + name + str(col)
		if output_verbose:
			line += "|Cstr:" + structure + "|Cseq:" + sequenceconstraint + "|Alpha:" + str(alpha) + "|Beta:" + str(beta) + "|tGC:" +  str(GC)  + "|ER:" + str(evaporation_rate) + "|Struct_CT:" + str(struct_correction_term) + "|GC_CT:" + str(GC_correction_term) + "|Seq_CT:" + str(seq_correction_term) + output_v + "\n" + "\n".join(output_w)  
		else:
			line += "\n" + output_w[1]
		if return_mod == False:
			if print_to_STDOUT:
				print line
			else:
				if col == 0:
					print2file(file_id, line, 'w')
				else:
					print2file(file_id, line, 'a')
		else:
			result.append(line)

	if return_mod == True:
		return result
	if print_to_STDOUT == False:    
		os.chdir(curr_dir)
  
def execute(args):
	"""
		CHECK AND PARSE THE COMMAND LINE STUFF
	"""
  
	# Checking the arguments, parsed from the shell
	###############################################
	name = args.name
	structure = args.Cstr
	
	if args.Cseq == "":
		sequenceconstraint = "N" * len(structure)
	else:
		sequenceconstraint = args.Cseq
	
	seed = args.seed

		
	alpha = args.alpha
	beta = args.beta
	tGC = args.tGC
	if tGC < 0 or tGC > 1:
		print "Error: Chosen tGC not in range [0,1]"
		exit(1)
	evaporation_rate = args.ER
	struct_correction_term = args.Cstrweight
	GC_correction_term = args.Cgcweight
	seq_correction_term = args.Cseqweight
	colonies = args.noOfColonies
	degreeOfSequenceInducement = args.level
	file_id = args.output_file
	verbose = args.verbose
	output_verbose = args.output_verbose
	
	tGCmax = args.tGCmax
	tGCvar = args.tGCvar
  
	termination_convergence = args.antsTerConv
	convergence_count = args.ConvergenceCount
	temperature = args.temperature
	reset_limit = args.Resets
	
	improve = args.improve_procedure

	### RNAfold parameterfile
	paramFile = args.paramFile
	
	# Using the pkiss program under user changeable parameters
	pseudoknots = args.pseudoknots
	
	# Loading the optimized parameters for pseudoknots and ignore user input
	if args.pseudoknot_parameters:
		alpha = 1.0
		beta = 0.1
		evaporation_rate = 0.2 
		struct_correction_term = 0.1 
		GC_correction_term = 1.0 
		seq_correction_term = 0.5 
		termination_convergence = 50 
		convergence_count = 100


	strategy = args.strategy
	useGU = args.useGUBasePair

	checkForViennaTools()
	if pseudoknots:
		checkForpKiss()
	findSequence(structure, sequenceconstraint, tGC, colonies, name, alpha, beta, evaporation_rate, struct_correction_term, GC_correction_term, seq_correction_term, degreeOfSequenceInducement, file_id, verbose, output_verbose, tGCmax, tGCvar, termination_convergence, convergence_count, reset_limit, improve, seed, temperature, paramFile, pseudoknots,  strategy, useGU)
  
  
def exe():
	"""
	MAIN EXECUTABLE WHICH PARSES THE INPUT LINE
	"""

	argument_parser = argparse.ArgumentParser(
	description = """
    
	#########################################################################
	#       antaRNA - ant assembled RNA                                     #
	#       -> Ant Colony Optimized RNA Sequence Design                     #
	#       ------------------------------------------------------------    #
	#       Robert Kleinkauf (c) 2015                                       #
	#       Bioinformatics, Albert-Ludwigs University Freiburg, Germany     #
	#########################################################################
  
	- For antaRNA only the VIENNNA RNA Package must be installed on your linux system.
	  antaRNA will only check, if the executables of RNAfold and RNAdistance of the ViennaRNA package can be found. If those programs are 
	  not installed correctly, no output will be generated, an also no warning will be prompted.
	  So the binary path of the Vienna Tools must be set up correctly in your system's PATH variable in order to run antaRNA correctly!
   
    - antaRNA was only tested under Linux.
    
    - For questions and remarks please feel free to contact us at http://www.bioinf.uni-freiburg.de/

	""",
	
	epilog = """   
	Example calls:
		python antaRNA.py --Cstr "...(((...)))..." --tGC 0.5 -n 2
		python antaRNA.py --Cstr ".........AAA(((...)))AAA........." --tGC 0.5 -n 10 --output_file /path/to/antaRNA_TESTRUN -ov
		python antaRNA.py --Cstr "BBBBB....AAA(((...)))AAA....BBBBB" --Cseq "NNNNANNNNNCNNNNNNNNNNNGNNNNNNUNNN" --tGC 0.5 -n 10

	#########################################################################
	#       --- Hail to the Queen!!! All power to the swarm!!! ---          #
	#########################################################################
		""",
	#formatter_class=RawTextHelpFormatter
	)
  
	# mandatorys
	argument_parser.add_argument("-Cstr", "--Cstr", help="Structure constraint using RNA dotbracket notation with fuzzy block constraint. \n(TYPE: %(type)s)\n\n", type=str, required=True)
	argument_parser.add_argument("-tGC", "--tGC", help="Objective target GC content in [0,1].\n(TYPE: %(type)s)\n\n", type=float, required=True)
	argument_parser.add_argument("-n", "--noOfColonies", help="Number of sequences which shall be produced. \n(TYPE: %(type)s)\n\n\n\n", type=int,  default=1)
	argument_parser.add_argument("-GU", "--useGUBasePair", help="Allowing GU base pairs. \n\n", action="store_true")
	
	argument_parser.add_argument("-s", "--seed", help = "Provides a seed value for the used pseudo random number generator.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=str, default="none")
	argument_parser.add_argument("-ip", "--improve_procedure", help = "Select the improving method.  h=hierarchical, s=score_based.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=str, default="s")  
	argument_parser.add_argument("-r", "--Resets", help = "Amount of maximal terrain resets, until the best solution is retuned as solution.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=int, default=5)  
	argument_parser.add_argument("-CC", "--ConvergenceCount", help = "Delimits the convergence count criterion for a reset.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=int, default=130)  
	argument_parser.add_argument("-aTC", "--antsTerConv", help = "Delimits the amount of internal ants for termination convergence criterion for a reset.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=int, default=50)
	
	argument_parser.add_argument("-p", "--pseudoknots", help = "Switch to pseudoknot based prediction using pKiss. Check the pseudoknot parameter usage!!!\n\n", action="store_true")
	argument_parser.add_argument("-pkPar", "--pseudoknot_parameters", help = "Enable optimized parameters for the usage of pseudo knots (Further parameter input ignored).\n\n", action="store_true")
	argument_parser.add_argument("--strategy", help = "Defining the pKiss folding strategy.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=str, default="A")
	
	# Mutual Exclusiv target GC distribution variables
	#tGCgroup = argument_parser.add_mutually_exclusive_group()
	argument_parser.add_argument("-tGCmax", "--tGCmax", help = "Provides a maximum tGC value [0,1] for the case of uniform distribution sampling. The regular tGC value serves as minimum value.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=-1.0)
	argument_parser.add_argument("-tGCvar", "--tGCvar", help = "Provides a tGC variance (sigma square) for the case of normal distribution sampling. The regular tGC value serves as expectation value (mu).\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=-1.0)
	
	argument_parser.add_argument("-t", "--temperature", help = "Provides a temperature for the folding algorithms.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=37.0)
	argument_parser.add_argument("-P", "--paramFile", help = "Changes the energy parameterfile of RNAfold. If using this explicitly, please provide a suitable energy file delivered by RNAfold. \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=str, default="")
	argument_parser.add_argument("-of","--output_file", help="Provide a path and an output file, e.g. \"/path/to/the/target_file\". \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=str, default="STDOUT")
	argument_parser.add_argument("-Cseq", "--Cseq", help="Sequence constraint using RNA nucleotide alphabet {A,C,G,U} and wild-card \"N\". \n(TYPE: %(type)s)\n\n", type=str, default = "")  
	argument_parser.add_argument("-l", "--level", help="Sets the level of allowed influence of sequence constraint on the structure constraint [0:no influence; 3:extensive influence].\n(TYPE: %(type)s)\n\n", type=int, default = 1)
	argument_parser.add_argument("--name", help="Defines a name which is used in the sequence output. \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=str, default="antaRNA_")
	argument_parser.add_argument("-a", "--alpha", help="Sets alpha, probability weight for terrain pheromone influence. [0,1] \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=1.0)
	argument_parser.add_argument("-b", "--beta", help="Sets beta, probability weight for terrain path influence. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=1.0)
	argument_parser.add_argument("-er", "--ER", help="Pheromone evaporation rate. \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=0.2)
	argument_parser.add_argument("-Cstrw", "--Cstrweight", help="Structure constraint quality weighting factor. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=0.5)
	argument_parser.add_argument("-Cgcw", "--Cgcweight", help="GC content constraint quality weighting factor. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", type=float, default=5.0)
	argument_parser.add_argument("-Cseqw", "--Cseqweight", help="Sequence constraint quality weighting factor. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n\n", type=float, default=1.0)
	argument_parser.add_argument("-ov", "--output_verbose", help="Displayes intermediate output.\n\n", action="store_true") 
	argument_parser.add_argument("-v", "--verbose", help="Prints additional features and stats to the headers of the produced sequences. Also adds the structure of the sequence.\n\n", action="store_true")
	
	args = argument_parser.parse_args()

	execute(args)
  
def checkForViennaTools():
	"""
	Checking for the presence of the Vienna tools in the system by which'ing for RNAfold and RNAdistance
	"""
	RNAfold_output = subprocess.Popen(["which", "RNAfold"], stdout=subprocess.PIPE).communicate()[0].strip()
	if len(RNAfold_output) > 0 and RNAfold_output.find("found") == -1 and RNAfold_output.find(" no ") == -1:
		return True
	else:
		print "It seems the Vienna RNA Package is not installed on your machine. Please do so!"
		print "You can get it at http://www.tbi.univie.ac.at/"
		exit(0)

		
def checkForpKiss():
	"""
		Checking for the presence of pKiss
	"""
  	pKiss_output = subprocess.Popen(["which", "pKiss"], stdout=subprocess.PIPE).communicate()[0].strip()
	if len(pKiss_output) > 0 and pKiss_output.find("found") == -1 and pKiss_output.find(" no ") == -1:
		return True
	else:
		print "It seems that pKiss is not installed on your machine. Please do so!"
		print "You can get it at http://bibiserv2.cebitec.uni-bielefeld.de/pkiss"
		exit(0)
		
		
		
if __name__ == "__main__":

	exe()
    
