#!/usr/bin/python
# # Gangman Yi
# Ver 0.1 - modified on 2017.04.04 
# Modified for use in ReGSP by Ejigu Girum
# Modified on September 2021

import os
import subprocess
import shutil
import time
#from cPlot7 import *
#max_allowed_hits=0
class DNA: 
    """Class representing DNA as a string sequence.""" 
 
    basecomplement = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n':'n'} 

    standard = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
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
                 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
                 }
 
    def __init__(self, s = "", name = ""): 
        """Create DNA instance initialized to string s.""" 
        self.seq = s.lower()                                                      # change to lower case
        self.seq = self.cleandna(self.seq)                                        # clean dna is below
        self.len = len(self.seq)
        self.name = name

    def cleandna(self, s):
        """Return dna only composed by letters ['a', 'c', 'g', 't', 'n']."""
        new_sequence = ""
        nucleotides = ['a', 'c', 'g', 't', 'n']
        for c in s:
            if c not in nucleotides:                    # If letter is not of type nucleotide continie to next
                continue
            new_sequence += c
        return new_sequence

    def getname(self):
        """Return the name of the sequence."""
        if self.name:
            return self.name
        else:
            return 'unknown'
     
    def getsequence(self):
        """Return the dna sequence."""
        return self.seq

    def setname(self, name):
        """Set the name of the sequence."""
        self.name = name

    def setsequence(self, s):
        """Set the sequence content."""
        self.seq = s.lower()
        self.seq = self.cleandna(self.seq)
        self.len = len(self.seq)
     
    def transcribe(self): 
        """Return as rna string.""" 
        return self.seq.replace('t', 'u') 
     
    def reverse(self): 
        """Return dna string in reverse order.""" 
        letters = list(self.seq) 
        letters.reverse() 
        return ''.join(letters) 
     
    def complement(self): 
        """Return the complementary dna string.""" 
        comp = ''
        letters = list(self.seq) 
        for base in letters:
            comp += self.basecomplement[base]
        return comp 
     
    def reversecomplement(self): 
        """Return the reverse complement of the dna string.""" 
        revcomp = ''
        letters = list(self.seq) 
        letters.reverse() 
        for base in letters:
            revcomp += self.basecomplement[base] 
        return revcomp 
     
    def gc_percentage(self): 
        """Return the percentage of dna composed of G+C.""" 
        s = self.seq 
        gc = s.count('g') + s.count('c') 
        return gc * 100.0 / len(s) 
 
    def translate(self, frame = 1):
        """ translate a DNA like cDNA sequence to a protein 
            possible frames 1,2,3,-1,-2,-3 
        """
        possibleframe=(1, 2, 3, -1, -2, -3)
        if frame not in possibleframe:
            frame = 1 # First frame
        if frame < 0 :
            cdna = self.reversecomplement()
            frame = abs(frame) - 1
        else:
            cdna = self.seq
            frame = frame - 1
        code = self.standard
        prot = ""
        i = frame # Starting frame
        while i <= len(cdna) - 3: # While there are at least 3 letters
            prot += code.get(cdna[i:i+3], "W")
            i += 3
        return prot

def READ_SEQ(fileName):

	new_db = {}
	db_idx = {}
	seq_id = ""

	fp = open(fileName, "r")

	for lineStr in fp:

		if lineStr.strip() == "" :
			continue

		if lineStr[0] == ">":
			seq_id = lineStr[1:].strip()
			new_db[seq_id] = ""
			db_idx[len(db_idx)] = seq_id
			continue

		new_db[seq_id] += lineStr		

	fp.close()

	return new_db, db_idx




def SPLIT_SEQ(seqDB, dbIDX, numSplit, outputFileName):

	totalNumSeq   = len(seqDB)
	seqNumPerFile = totalNumSeq / numSplit

	for i in range(0, numSplit):

		startIDX = i * seqNumPerFile 
		endIDX   = startIDX + seqNumPerFile

		if i == numSplit - 1:
			endIDX += totalNumSeq % numSplit
			
		fp_out = open(outputFileName+"."+str(i+1),"w")
		for j in range(startIDX, endIDX):
			fp_out.write(seqDB[dbIDX[j]]+"\n")
		
		fp_out.close()


def READ_BLAST_SEQ_ALIGNMENT(fp):
	# this function reads sequence alignments and returns them

	result = ""
	frame  = 1
	alg_length = 0
	start_pos  = -1
	end_pos    = -1  

	while 1:

		lineStr = fp.readline()

		if lineStr.strip() == "":
			continue

		if "Length=" == lineStr.strip()[0:7]:
			alg_length = int(lineStr.strip().split("=")[1])
			continue

		if "Score = " == lineStr.strip()[0:8]:
			result += lineStr.strip()+"\n\t\t       "
			lineStr = fp.readline()
			result += lineStr.strip()+", 	"
			lineStr = fp.readline()
			result += lineStr.strip()+"\n\n"
			frame   = lineStr.strip().split()[2]

			count = 0
			while 1:
				lineStr = fp.readline()
				if lineStr.strip() == "":
					count += 1					
					if count == 2:
						return frame, int(start_pos), int(end_pos), result
					else:
						continue
				if lineStr.strip()[0:6] == "Query ":
					tmp_arr = lineStr.strip().split()
					if start_pos == -1 :
						start_pos = int(tmp_arr[1])
					end_pos   = int(tmp_arr[-1])

				result += "\t\t\t"+lineStr
				count   = 0


def SHOW_ALIGNMENTS_BY_FRAME_ORDER(frame, dna_seq):

	print ("\tFrame = %s" % frame)
	print ("\tSeq = %s" % (dna_seq.translate(frame)))


def SPLIT_BY_LEN(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]


def GET_NUM_CORES():

	cmd = "uname -a"
	cmd_out = subprocess.getoutput(cmd).strip()
	
	if "Darwin" in cmd_out:		
		# Mac OS
		cmd = "sysctl -a | grep machdep.cpu.core_count"
		cmd_out = subprocess.getoutput(cmd).split()[1].strip()
	elif "Linux" in cmd_out:
		# Linux 
		cmd = "grep -c \"processor\" /proc/cpuinfo"
		cmd_out = subprocess.getoutput(cmd).strip()
	else:
		cmd_out = 1
	
	return int(cmd_out)
	
	
def SHOW_SEQUENCE_WTH_POS_NUM(seq_id, seq, tab=0):

	new_seq = ""	
	output  = ""
	gap     = 10
	indent  = 5

	for ch in seq:
		if ch.lower() in ['a', 'c', 'g', 't', 'n']:
			new_seq += ch

	if seq_id != "":
		output += ("Query = "+seq_id+"\n")
	bar = "\t"*tab+" "*indent
	for i in range(0, 11):
		if i == 0 :
			bar += ("."+" "*(gap-2))
		else:
			if i == 5 :
				bar += ("*"+" "*(gap-1))
			elif i % 2 == 0 :
				bar += ("."+" "*(gap-1))
			else:
				bar += (":"+" "*(gap-1))
	
	output += (bar+"\n")
	tmp_arr = SPLIT_BY_LEN(new_seq, 100)
	for i in range(0, len(tmp_arr)):
		output += "\t"*tab+("%-5s%s\n" % (i,tmp_arr[i] ))

	return output


def CHECK_GROUPS(groups, start_pos, end_pos,max_allowed_hits):

	if start_pos > end_pos :
		tmp = start_pos
		start_pos = end_pos
		end_pos = tmp

	flag = 0
	for g_num in groups:
		g_start = int(groups[g_num][0].split()[0])
		g_end   = int(groups[g_num][0].split()[1])
		if start_pos <= g_start and end_pos >= g_start:
			flag = 1
			break
		if start_pos >= g_start and end_pos <= g_end :
			flag = 1
			break
		if end_pos >= g_end and start_pos <= g_end :
			flag = 1
			break

	if flag == 1:
		if max_allowed_hits <= groups[g_num][1] :
			return groups, 0
		else:
			groups[g_num][1] += 1
			return groups, 1

	max_groups = len(groups)
	groups[max_groups] = {}
	groups[max_groups][0] = "%s %s" % ( start_pos, end_pos )
	groups[max_groups][1] = 1

	return groups, 1


def REFINE_POS(seq, start_pos, end_pos):

	if start_pos > end_pos:
		tmp       = start_pos
		start_pos = end_pos
		end_pos   = tmp

	if start_pos == 0:
		start_pos = 1

	if end_pos > len(seq):
		end_pos = len(seq)
		
	return start_pos, end_pos
	
	
def GET_FRAGMENT(seq, start_pos, end_pos):

	start_pos, end_pos = REFINE_POS(seq, start_pos, end_pos)

	return SHOW_SEQUENCE_WTH_POS_NUM("", seq[start_pos-1:end_pos].upper(), 2)


def REPLACE_FRAGMENT(seq, start_pos, end_pos):

	start_pos, end_pos = REFINE_POS(seq, start_pos, end_pos)

	return seq[:start_pos-1]+("N"*(end_pos-start_pos+1))+seq[end_pos:]


def READ_BLAST_OUT(org_seq, seq, max_allowed_hits, blastFileName):
	# this function reads a blast output file and parses queries and their matching positions.
	# 6 frame translation is applied for checking different matching positions.
	#print(os.getcwd())
	fp = open(blastFileName, "r")

	database = ""
	start_pos = 0
	end_pos   = 0
	groups    = {}
	results   = ""
	nohits    = 0

	while 1:

		lineStr = fp.readline()

		if not lineStr:
			break

		if lineStr.strip() == "":
			continue

		if lineStr[0:6] == "Query=":
			if results != "":
				groups = {}
			continue

		if lineStr[0] == ">":
			nohits += 1
			database = lineStr[1:].strip()
			frame, start_pos, end_pos, seq_alignment = READ_BLAST_SEQ_ALIGNMENT(fp)

			groups, flag = CHECK_GROUPS(groups, start_pos, end_pos,max_allowed_hits)
			if flag == 1 :
				results += "\tDB = %s \n" % (database)
				results +=  GET_FRAGMENT(org_seq, start_pos, end_pos)				
				seq = REPLACE_FRAGMENT(seq, start_pos, end_pos)
				results += "\t\t%sPos = [%s,%s], Len = %s, %s\n"  % (frame[0], start_pos, end_pos, abs(end_pos-start_pos)+1, seq_alignment)

	fp.close()

	return seq, results, nohits


def RUN_BLAST_OUT(cores, seq, query, dbFileName, blastFileName):

	seqFileName = blastFileName+".seq"
	fp = open(seqFileName, "w")
	fp.write(">%s\n%s" % (query, seq.upper()))
	fp.close()
	cmd = "blastx -db %s -query %s -num_descriptions 250 -num_alignments 250 -num_threads %s -out %s " % ( dbFileName, seqFileName, cores, blastFileName)
	cmd_out = subprocess.getoutput(cmd)	
	cmd = "rm -rf "+seqFileName
	cmd_out = subprocess.getoutput(cmd)


def PARSE_BLAST_OUT(cores, seq, query, fp_out, dbFileName, outputFileName,max_allowed_hits):
	num_hits   = 1
	hit_counts = 0

	blastFileName = outputFileName+".blastx+"
	org_seq = seq 

	result_title = SHOW_SEQUENCE_WTH_POS_NUM(query, seq.upper())+"\n"
	results = ""
	fp_out.write(result_title)
	
	while 1:
		RUN_BLAST_OUT(cores, seq, query, dbFileName, blastFileName)
		seq, tmp_result, num_hits = READ_BLAST_OUT(org_seq, seq, max_allowed_hits, blastFileName)
		results += tmp_result
		fp_out.write(tmp_result)
		fp_out.flush()
		
		if num_hits == 0 :
			break
		
		hit_counts += num_hits
	
	

	cmd = "rm -rf "+blastFileName
	cmd_out = subprocess.getoutput(cmd)
	fp_out.write(tmp_result)	
	
	result_title = result_title.strip()+"\n[Total DB counts = %d]\n\n" % hit_counts
	results = result_title+results

	return hit_counts, results

def FIND_SUBREGION_CTG(seq_db, seq_db_idx, dbFileName, outputFileName,max_allowed_hits): ## dbFileName chnged to a list for mult ref ['ref1', 'ref2']
	
	# find blastx hits on ctgs until no more hits.
    fp_out = open(outputFileName, "w")
    cores = GET_NUM_CORES()

    raw_data = {}
    raw_hits = {}

	#widgets = ['   working: ', Percentage(), ' ', Bar(marker='=',left='[',right=']'), ' ', ETA(), ' ']
	#pbar = ProgressBar(widgets=widgets, maxval=len(seq_db_idx))
	#pbar.start()
    for ridx in range(0,len(dbFileName)):
        cmd = "makeblastdb -dbtype prot -in %s " % (dbFileName[ridx])
        cmd_out = subprocess.getoutput(cmd)	
        raw_data[ridx] = {}
        raw_hits[ridx] = {}
        for idx in range(0,len(seq_db_idx)):
            query   = seq_db_idx[idx]
            dna_seq = DNA(seq_db[query], query)
            num_hits, tmp_result = PARSE_BLAST_OUT(cores, dna_seq.getsequence(), query, fp_out, dbFileName[ridx], outputFileName, max_allowed_hits)
            if num_hits not in raw_hits[ridx]:
                raw_hits[ridx][num_hits] = {}
            raw_hits[ridx][num_hits][idx] = 1
            raw_data[ridx][idx] = tmp_result 	

            #pbar.update(idx)
        #print(raw_hits)
        #pbar.finish()

        cmd ="rm -rf %s.p?? " % (dbFileName[ridx])
        cmd_out = subprocess.getoutput(cmd)	
    fp_out.close()
    return raw_hits, raw_data

	
def HANDLE_RAW_DATA(seq_db, seq_db_idx, raw_hits, raw_data, refFileName, outputFileName, minDBcount,directory): # refFileName becomes list

	fp_out = open(os.path.join(directory, outputFileName+".sorted"), "w")
	fp_seq_out = open(os.path.join(directory,outputFileName+".seq"), "w")
	plot_seq_file = os.path.join(directory,outputFileName+".plot.seq")
	fp_adj_seq_out = open(plot_seq_file, "w")

	regsp=open(os.path.join(directory,outputFileName+".regsp"),"w")


	for ridx in raw_hits.keys():
		for i in sorted(raw_hits[ridx].keys(),  reverse=True):
			for idx in raw_hits[ridx][i]:
				fp_out.write(raw_data[ridx][idx])
				query = seq_db_idx[idx]
				dna_seq = DNA(seq_db[query], query)
				fp_seq_out.write(">%s\n%s\n\n" % ( query.split("_")[-1], dna_seq.getsequence() ))
				if(minDBcount<i):
					regsp.write(raw_data[ridx][idx])
	fp_out.close()
	fp_seq_out.close()
	regsp.close()
    #for ridx in raw_hits.keys(): ##we are here
	for i in sorted(raw_hits[0].keys(),  reverse=True):

		if i < minDBcount:
			break	

		for idx in raw_hits[0][i]:
			query = seq_db_idx[idx]
			dna_seq = DNA(seq_db[query], query)
			fp_adj_seq_out.write(">%s\n%s\n\n" % ( query, dna_seq.getsequence().upper()))
	
	fp_adj_seq_out.close()
	refFileName2=joinfiles(refFileName)
	refFileName2='Ref_file.fasta'
	cmd = "promer -p %s %s %s" % (outputFileName+".mummer", refFileName2, plot_seq_file)
	cmd_out = subprocess.getoutput(cmd)		
	cmd = "mummerplot -postscript -p %s %s " % ( outputFileName+".plot", outputFileName+".mummer.delta")
	cmd_out = subprocess.getoutput(cmd)		
	cmd = "ps2pdf %s %s" % ( outputFileName+".plot.ps", os.path.join(directory,outputFileName+".plot.pdf"))
	cmd_out = subprocess.getoutput(cmd)		
	cmd = "rm -rf %s.delta %s.pl  ot.ps %s.plot.fplot %s.plot.rplot %s.plot.gp %s.plot.ps %s" % ( outputFileName+".mummer", outputFileName, outputFileName, outputFileName, outputFileName, outputFileName, outputFileName)
	cmd_out = subprocess.getoutput(cmd)
    
    	
	
def joinfiles(filelist):
    with open('Ref_file.fasta','wb') as wfd:
        for f in filelist:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return wfd



	
	
