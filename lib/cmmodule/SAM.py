#!/usr/bin/env python
'''manipulate SAM file.'''

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
import random

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
import pysam
from cmmodule import mystat
from cmmodule import fasta
from cmmodule import cigar
from cmmodule import BED
#changes to the paths


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2010, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production


class ParseSAM:
	'''This class provides fuctions to parsing/processing/transforming SAM format file
	Format of SAM file see: http://samtools.sourceforge.net/SAM1.pdf'''
	_reExpr1=re.compile(r'\s+')				#>=1 spaces
	_reExpr2=re.compile(r'^\s*$')			#blank line
	_splicedHit_pat = re.compile(r'(\d+)[M|N]',re.IGNORECASE)	#regular expression for spliced mapped reads
	_monoHit_pat = re.compile(r'^(\d+)M$',re.IGNORECASE)				#regular expresion for Non-spliced mapped reads
	_insertionHit_pat =re.compile(r'\d+I',re.IGNORECASE)
	_deletionHit_pat =re.compile(r'\d+D',re.IGNORECASE)
	_softClipHi_pat =re.compile(r'\d+S',re.IGNORECASE)
	_hardClipHi_pat =re.compile(r'\d+H',re.IGNORECASE)
	_padHit_pat =re.compile(r'\d+P',re.IGNORECASE)
	_uniqueHit_pat = re.compile(r'[NI]H:i:1\b')

	def __init__(self,samFile):
		'''constructor'''
		if samFile == '-':
			self.fileName = "STDIN"
			self.f = sys.stdin
		else:
			self.fileName = os.path.basename(samFile)
			self.f=open(samFile,'r')

	def stat (self):
		'''Calculate mapping statistics'''		
		total_read=0
		
		pcr_duplicate =0
		low_qual =0
		secondary_hit =0
		
		unmapped_read1=0
		mapped_read1=0
		reverse_read1=0
		forward_read1=0
		
		unmapped_read2=0
		mapped_read2=0
		reverse_read2=0
		forward_read2=0
		
		_numSplitHit =0
		_numMonoHit =0
		_numInsertion =0
		_numDeletion =0

		minus_minus=0
		minus_plus =0
		plus_minus=0
		plus_plus=0
		paired=True
		
		unmap_SE=0
		map_SE=0
		reverse_SE=0
		forward_SE=0
		for line in self.f:
			line=line.rstrip()
			if line.startswith('@'):continue				#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			
			total_read +=1
			field = line.split()
			flagCode=string.atoi(field[1])
			if  (flagCode & 0x0400 !=0):							#PCR or optical duplicate
				pcr_duplicate +=1
				continue
			if  (flagCode & 0x0200 !=0):							#Low quality
				low_qual +=1
				continue
			if  (flagCode & 0x0200 !=0):							#Not primary alignment
				secondary_hit +=1
				continue
			if (len(ParseSAM._splicedHit_pat.findall(field[5]))>1):_numSplitHit +=1			#Splicing mapped reads							
			if (len(ParseSAM._splicedHit_pat.findall(field[5]))==1):_numMonoHit +=1			#mono mapped reads					
			if (ParseSAM._insertionHit_pat.search(field[5])):_numInsertion +=1				#insertion in reads
			if (ParseSAM._deletionHit_pat.search(field[5])):_numDeletion +=1				#deletion in reads
			


			if (flagCode & 0x0001 !=0):								#This is paired end sequencing
				if (flagCode & 0x0040 != 0):						#1st read
					if (flagCode & 0x0004 != 0):unmapped_read1 +=1
					if (flagCode & 0x0004 == 0):mapped_read1 +=1
					if (flagCode & 0x0010 != 0):reverse_read1 +=1
					if (flagCode & 0x0010 == 0):forward_read1 +=1
						
				if (flagCode & 0x0080 != 0):						#2nd read
					if (flagCode & 0x0004 != 0):unmapped_read2 +=1
					if (flagCode & 0x0004 == 0):mapped_read2 +=1
					if (flagCode & 0x0010 != 0):reverse_read2 +=1
					if (flagCode & 0x0010 == 0):forward_read2 +=1						
				if 	(flagCode & 0x0010 != 0 and flagCode & 0x0020 != 0):
					minus_minus +=1
				if 	(flagCode & 0x0010 != 0 and flagCode & 0x0020 == 0):
					minus_plus +=1
				if 	(flagCode & 0x0010 == 0 and flagCode & 0x0020 != 0):
					plus_minus +=1
				if 	(flagCode & 0x0010 == 0 and flagCode & 0x0020 == 0):
					plus_plus +=1				
			if (flagCode & 0x0001 ==0):								#This is single end sequencing
				paired=False
				if (flagCode & 0x0004 != 0):
					unmap_SE +=1
				if (flagCode & 0x0004 == 0):
					map_SE +=1
				if (flagCode & 0x0010 != 0):
					reverse_SE +=1
				if (flagCode & 0x0010 == 0):
					forward_SE +=1
				
		if paired:				
			print >>sys.stderr,"\n#=================================================="
			print >>sys.stderr,"#================Report (pair-end)================="
			print >>sys.stderr, "%-25s%d" % ("Total Reads:",total_read)
			print >>sys.stderr, "%-25s%d" % ("Total Mapped Reads:", (mapped_read1 + mapped_read2))
			print >>sys.stderr, "%-25s%d" % ("Total Unmapped Reads:",(unmapped_read1 + unmapped_read2))
			print >>sys.stderr, "%-25s%d" % ("PCR duplicate:",pcr_duplicate)
			print >>sys.stderr, "%-25s%d" % ("QC-failed:",low_qual)
			print >>sys.stderr, "%-25s%d" % ("Not primary mapping:",secondary_hit)
			print >>sys.stderr, "\n",
			print >>sys.stderr, "%-25s%d" % ("Unmapped Read-1:",unmapped_read1)
			print >>sys.stderr, "%-25s%d" % ("Mapped Read-1:",mapped_read1)
			print >>sys.stderr, "%-25s%d" % ("  Forward (+):",forward_read1)
			print >>sys.stderr, "%-25s%d" % ("  Reverse (-):",reverse_read1)
			
			print >>sys.stderr, "\n",
			print >>sys.stderr, "%-25s%d" % ("Unmapped Read-2:",unmapped_read2)
			print >>sys.stderr, "%-25s%d" % ("Mapped Read-2:",mapped_read2)
			print >>sys.stderr, "%-25s%d" % ("  Forward (+):",forward_read2)
			print >>sys.stderr, "%-25s%d" % ("  Reverse (-):",reverse_read2)
			
			print >>sys.stderr, "\n",
			print >>sys.stderr, "%-25s%d" % ("Mapped to (+/-):",plus_minus)
			print >>sys.stderr, "%-25s%d" % ("Mapped to (-/+):",minus_plus)
			print >>sys.stderr, "%-25s%d" % ("Mapped to (+/+):",plus_plus)
			print >>sys.stderr, "%-25s%d" % ("Mapped to (-/-):",minus_minus)
			print >>sys.stderr, "\n",
			print >>sys.stderr, "%-25s%d" % ("Spliced Hits:",_numSplitHit)
			print >>sys.stderr, "%-25s%d" % ("Non-spliced Hits:",_numMonoHit)
			print >>sys.stderr, "%-25s%d" % ("Reads have insertion:",_numInsertion)
			print >>sys.stderr, "%-25s%d" % ("Reads have deletion:",_numDeletion)
		else:
			print >>sys.stderr,"\n#===================================================="
			print >>sys.stderr,"#================Report (single-end)================="
			print >>sys.stderr, "%-25s%d" % ("Total Reads:",total_read)
			print >>sys.stderr, "%-25s%d" % ("Total Mapped Reads:", map_SE)
			print >>sys.stderr, "%-25s%d" % ("Total Unmapped Reads:",unmap_SE)
			print >>sys.stderr, "%-25s%d" % ("PCR duplicate:",pcr_duplicate)
			print >>sys.stderr, "%-25s%d" % ("QC-failed:",low_qual)
			print >>sys.stderr, "%-25s%d" % ("Not primary mapping:",secondary_hit)
			print >>sys.stderr, "%-25s%d" % ("froward (+):",forward_SE)
			print >>sys.stderr, "%-25s%d" % ("reverse (-):",reverse_SE)
			print >>sys.stderr, "\n",
			print >>sys.stderr, "%-25s%d" % ("Spliced Hits:",_numSplitHit)
			print >>sys.stderr, "%-25s%d" % ("Non-spliced Hits:",_numMonoHit)
			print >>sys.stderr, "%-25s%d" % ("Reads have insertion:",_numInsertion)
			print >>sys.stderr, "%-25s%d" % ("Reads have deletion:",_numDeletion)			
			
	def samTobed(self,outfile=None,mergePE=False):
		"""Convert SAM file to BED file. BED file will be saved as xxx.sam.bed unless otherwise specified.
		 If mergePE=False, each read will be one bed entry. If mergePE=True, pair-end (if there are and on
		 the same chr) reads will be displayed in bed entry."""
		if outfile is None:
			outfile=self.fileName + ".bed"

		print >>sys.stderr,"\tWriting bed entries to\"",outfile,"\"...",
		FO=open(outfile,'w')
		for line in self.f:
			if line.startswith(('@','track')):continue		#skip head lines
			if ParseSAM._reExpr2.match(line):continue	#skip blank line
			field=line.rstrip().split()
			if (string.atoi(field[1]) & 0x0004)!=0: continue	#skip unmapped line
			if (string.atoi(field[1]) & 0x0040)!=0:
				mate="/1"
			else:
				mate="/2"
			flag = string.atoi(field[1])
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			chrom = field[2]
			chromStart = string.atoi(field[3])-1
			chromEnd=chromStart
			for i in comb:
				chromEnd += i		
			name = field[0] + mate
			score = field[4]
			if(flag & 0x0010)==0:
				strand = '+'
			else:
				strand = '-'
			thickStart = chromStart
			thickEnd = chromEnd
			itemRgb = "0,255,0"
			blockCount = (len(comb) +1) /2
			blockSize = []
			for i in range(0,len(comb),2):
				blockSize.append(str(comb[i]))
			blockSizes = ','.join(blockSize)
			blockStart=[]
			for i in range(0,len(comb),2):
				blockStart.append(str(sum(comb[:i])))
			blockStarts = ','.join(blockStart)
			print >>FO, string.join((str(i) for i in [chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts]),sep="\t")			
		print >>sys.stderr, "Done"
		FO.close()
		self.f.seek(0)

		if mergePE:
		#creat another bed file. pair-end reads will be merged into single bed entry
			print >>sys.stderr, "Writing consoidated bed file ...",
			bedfile = open(outfile,'r')
			outfile_2 = outfile + ".consolidate.bed"
			outfile_3 = outfile + '.filter'
			FO = open(outfile_2,'w')
			FOF = open(outfile_3,'w')
			count={}
			chr=collections.defaultdict(set)
			txSt={}
			txEnd={}
			strand=collections.defaultdict(list)
			blocks={}
			sizes=collections.defaultdict(list)
			starts=collections.defaultdict(list)
			for line in bedfile:
				field=line.strip().split()
				if field[3] not in count:
					count[field[3]] = 1
					chr[field[3]].add(field[0])
					txSt[field[3]] = string.atoi(field[1])
					txEnd[field[3]] = string.atoi(field[2])
					strand[field[3]].append(field[5])
					blocks[field[3]] = string.atoi(field[9])
					sizes[field[3]].extend( field[10].split(',') )
					starts[field[3]].extend([string.atoi(i) + string.atoi(field[1]) for i in field[11].split(',') ])
				else:
					count[field[3]] += 1
					chr[field[3]].add(field[0])
					if string.atoi(field[1]) < txSt[field[3]]:
						txSt[field[3]] = string.atoi(field[1]) 
					if string.atoi(field[2]) > txEnd[field[3]]:
						txEnd[field[3]] = string.atoi(field[2]) 
					blocks[field[3]] += string.atoi(field[9])
					strand[field[3]].append(field[5])
					sizes[field[3]].extend( field[10].split(',') )
					starts[field[3]].extend([string.atoi(i) + string.atoi(field[1]) for i in field[11].split(',') ])				
			#befile.seek(0)
		
			for key in count:
				st=[]	#real sorted starts
				sz=[]	#real sorted sizes
				if(count[key] ==1):			#single-end read
					if(blocks[key] ==1):	#single end, single hit
						st = [i - txSt[key] for i in starts[key]]
						st = string.join([str(i) for i in st],',')
						print >>FO, chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key,"\t","11\t",strand[key][0],"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sizes[key],','),"\t",st
					else:
						st = [i - txSt[key] for i in starts[key]]	#single end, spliced hit
						st = string.join([str(i) for i in st],',')
						print >>FO, chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key,"\t","12\t",strand[key][0],"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sizes[key],','),"\t",st						
								
				elif(count[key]==2):	#pair-end read
					direction = string.join(strand[key],'/')	
					for i,j in sorted (zip(starts[key],sizes[key])):
						st.append(i-txSt[key])
						sz.append(j)
					#st=[string.atoi(i) for i in st]
					if(len(chr[key])==1):	#pair-end reads mapped to same chromosome
						if blocks[key] ==2:	#pair end, single hits
							print >>FO, chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key + "|strand=" + direction + "|chrom=same","\t","21\t",'.',"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sz,','),"\t",string.join([str(i) for i in st],',')
						elif blocks[key] >2:	#
							print >>FO, chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key + "|strand=" + direction + "|chrom=same","\t","22\t",'.',"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sz,','),"\t",string.join([str(i) for i in st],',')
					else:
						print >>FOF,key,"\t","pair-end mapped, but two ends mapped to different chromosome"
				elif(count[key] >2):	#reads occur more than 2 times
					print >>FOF,key,"\t","occurs more than 2 times in sam file"
					continue
			FO.close()
			FOF.close()
			print >>sys.stderr, "Done"


	def samTowig(self,outfile=None,log2scale=False,header=False,strandSpecific=False):
		"""Convert SAM file to wig file. WIG file will be saved as xxx.sam.wig unless otherwise specified.
		log2scale has no effect if strandSpecific was set"""
		
		if outfile is None:
			outfile = self.fileName + ".wig"
		FO=open(outfile,'w')
		print >>sys.stderr, "Writing wig file to\"",outfile,"\"..."
		
		headline="track type=wiggle_0 name=" + outfile + " track_label description='' visibility=full color=255,0,0"
		wig=collections.defaultdict(dict)
		Pwig=collections.defaultdict(dict)
		Nwig=collections.defaultdict(dict)
		
		#strand_rule={'1+':'+','1-':'-','2+':'-','2-':'+'}	
		strand_rule={'1+':'-','1-':'+','2+':'+','2-':'-'}
		
		for line in self.f:
			hits=[]
			if line.startswith('@'):continue							#skip head lines
			if ParseSAM._reExpr2.match(line):continue					#skip blank lines
			field=line.rstrip().split()
			flagCode = string.atoi(field[1])
			if (flagCode & 0x0004)!=0: continue			#skip unmapped line										
			if not ParseSAM._uniqueHit_pat.search(line):continue		#skip multiple mapped reads
			if (flagCode & 0x0100 !=0): continue						#skip non primary hit
			if (flagCode & 0x0200 !=0): continue						#skip QC-failed
			if (flagCode & 0x0400 !=0): continue						#skip PCR artifact
			
			chrom=field[2]
			txStart=string.atoi(field[3])
			if (flagCode & 0x0010 != 0):
				strand='-'
			else:
				strand='+'
			if (flagCode & 0x0040 !=0): read_type='1'
			if (flagCode & 0x0080 !=0): read_type='2'
			
			blocks = cigar.fetch_exon(chrom,txStart,field[5])	
			for block in blocks:
				hits.extend(range(block[1]+1,block[2]+1))
				
			if strandSpecific is not True:
				for i in hits:
					if wig[chrom].has_key(i):
						wig[chrom][i] +=1
					else:
						wig[chrom][i]=1
			else:
				if strand_rule[read_type + strand] == '-':
					for i in hits:
						if Nwig[chrom].has_key(i):
							Nwig[chrom][i] += 1
						else:
							Nwig[chrom][i] = 1
				if strand_rule[read_type + strand] == '+':
					for i in hits:
						if Pwig[chrom].has_key(i):
							Pwig[chrom][i] +=1
						else:
							Pwig[chrom][i]=1					
		
		if header:FO.write(headline + "\n")
		
		if strandSpecific is not True:
			for chr in sorted(wig.keys()):
				print >>sys.stderr, "Writing ",chr, " ..."
				FO.write('variableStep chrom='+chr+'\n')
				for coord in sorted(wig[chr]):
					if log2scale:FO.write("%d\t%5.3f\n" % (coord,math.log(wig[chr][coord],2)))
					else:FO.write("%d\t%d\n" % (coord,wig[chr][coord]))
		else:
			chroms=set(Pwig.keys() + Nwig.keys())
			for chr in sorted(chroms):
				print >>sys.stderr, "Writing ",chr, " ..."
				FO.write('variableStep chrom='+chr+'\n')
				coords=sorted(set(Pwig[chr].keys() + Nwig[chr].keys()))
				for coord in coords:
					if ((coord in Pwig[chr]) and (coord not in Nwig[chr])):
						FO.write("%d\t%d\n" % (coord,Pwig[chr][coord]))
					elif ((coord in Nwig[chr]) and (coord not in Pwig[chr])):
						FO.write("%d\t%d\n" % (coord,-Nwig[chr][coord]))
					elif ((coord in Nwig[chr]) and (coord in Pwig[chr])):
						if (Pwig[chr][coord] >= Nwig[chr][coord]):
							FO.write("%d\t%d\n" % (coord,Pwig[chr][coord]))
						else:
							FO.write("%d\t%d\n" % (coord,-Nwig[chr][coord]))				
					
		self.f.seek(0)
		FO.close()
		

	def getUnmap(self, outfile=None,fastq=True):
		'''Extract unmapped reads from SAM file and write to fastq [when fastq=True] or fasta [when fastq=False] file'''
		if outfile is None:
			if fastq: outfile = self.fileName + ".unmap.fq"
			else: outfile = self.fileName + ".unmap.fa"
		FO=open(outfile,'w')
		unmapCount=0
		print >>sys.stderr, "Writing unmapped reads to\"",outfile,"\"... ",
		
		for line in self.f:
			hits=[]
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip().split()	
			flagCode=string.atoi(field[1])
			seq=field[9]
			qual=field[10]
			if (flagCode & 0x0004) != 0:							#read unmap
				unmapCount +=1
				if (flagCode & 0x0001) != 0:							#paried in sequencing
					if (flagCode & 0x0040)!=0:seqID=field[0] + '/1'			#first read	
					if (flagCode & 0x0080)!=0:seqID=field[0] + '/2'			#second read
				else: seqID=field[0]
				
				if fastq: FO.write('@' + seqID + '\n' + seq +'\n' + '+' +'\n' + qual+'\n')
				else: FO.write('>' + seqID + '\n' + seq +'\n')

		print >>sys.stderr, str(unmapCount) + " reads saved!\n"
		FO.close()
		self.f.seek(0)
	
	
	def getProperPair(self,outfile=None):
		'''Extract proper paried mapped reads.'''
		if outfile is None:
			outfile = self.fileName + ".PP.sam"
		FO=open(outfile,'w')
		PPcount=0
		print >>sys.stderr, "Writing proper paired reads to\"",outfile,"\"... ",
		for line in self.f:
			hits=[]
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			if ((flagCode & 0x0001) != 0) and ((flagCode & 0x0002)!=0):
				PPcount +=1
				FO.write(line)
		FO.close()
		print >>sys.stderr, str(PPcount) + " reads were saved!\n",
		self.f.seek(0)

	def samNVC(self,outfile=None):
		'''for each read, calculate nucleotide frequency vs position'''
		if outfile is None:
			outfile1 = self.fileName + ".NVC.xls"
			outfile2 = self.fileName +".NVC_plot.r"
		else:
			outfile1 = outfile + ".NVC.xls"
			outfile2 = outfile +".NVC_plot.r"
		FO=open(outfile1,'w')
		RS=open(outfile2,'w')
		PPcount=0
		
		transtab = string.maketrans("ACGTNX","TGCANX")
		base_freq=collections.defaultdict(int)
		a_count=[]
		c_count=[]
		g_count=[]
		t_count=[]
		print >>sys.stderr, "reading sam file ... "
		for line in self.f:
			if line.startswith('@'):continue									#skip head lines	
			if ParseSAM._reExpr2.match(line):continue					#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			
			if flagCode & 0x0010 ==0:	#plus strand
				RNA_read = field[9].upper()
			else:
				RNA_read = field[9].upper().translate(transtab)[::-1]
			for i in xrange(len(RNA_read)):
				key = str(i) + RNA_read[i]
				base_freq[key] += 1
		
		print >>sys.stderr, "generating data matrix ..."
		print >>FO, "Position\tA\tC\tG\tT\tN\tX"
		for i in xrange(len(RNA_read)):
			print  >>FO, str(i) + '\t',
			print  >>FO, str(base_freq[str(i) + "A"]) + '\t',
			a_count.append(str(base_freq[str(i) + "A"]))
			print  >>FO, str(base_freq[str(i) + "C"]) + '\t',
			c_count.append(str(base_freq[str(i) + "C"]))
			print  >>FO, str(base_freq[str(i) + "G"]) + '\t',
			g_count.append(str(base_freq[str(i) + "G"]))
			print  >>FO, str(base_freq[str(i) + "T"]) + '\t',
			t_count.append(str(base_freq[str(i) + "T"]))
			print  >>FO, str(base_freq[str(i) + "N"]) + '\t',
			print  >>FO, str(base_freq[str(i) + "X"]) + '\t'
		FO.close()
		
		#generating R scripts
		print >>sys.stderr, "generating R script  ..."
		print >>RS, "position=c(" + ','.join([str(i) for i in xrange(len(RNA_read))]) + ')'
		print >>RS, "A_count=c(" + ','.join(a_count) + ')'
		print >>RS, "C_count=c(" + ','.join(c_count) + ')'
		print >>RS, "G_count=c(" + ','.join(g_count) + ')'
		print >>RS, "T_count=c(" + ','.join(t_count) + ')'
		print >>RS, "total= A_count + C_count + G_count + T_count"
		print >>RS, "ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05"
		print >>RS, "yn=min(A_count/total,C_count/total,G_count/total,T_count/total)"
		
		print >>RS, 'pdf("NVC_plot.pdf")'
		print >>RS, 'plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")'
		print >>RS, 'lines(position,T_count/total,type="o",pch=20,col="red")'
		print >>RS, 'lines(position,G_count/total,type="o",pch=20,col="blue")'
		print >>RS, 'lines(position,C_count/total,type="o",pch=20,col="cyan")'
		print >>RS, 'legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))'
		print >>RS, "dev.off()"
		
		RS.close()
		#self.f.seek(0)

	def samGC(self,outfile=None):
		'''GC content distribution of reads'''
		if outfile is None:
			outfile1 = self.fileName + ".GC.xls"
			outfile2 = self.fileName +".GC_plot.r"
		else:
			outfile1 = outfile + ".GC.xls"
			outfile2 = outfile +".GC_plot.r"
		FO=open(outfile1,'w')
		RS=open(outfile2,'w')
		
		gc_hist=collections.defaultdict(int)	#key is GC percent, value is count of reads
		print >>sys.stderr, "reading sam file ... "
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			gc_percent = "%4.2f" % ((field[9].upper().count('C') + field[9].upper().count('G'))/(len(field[9])+0.0)*100)
			#print gc_percent
			gc_hist[gc_percent] += 1
		
		print >>sys.stderr, "writing GC content ..."
		
		print >>FO, "GC%\tread_count"
		for i in gc_hist.keys():
			print >>FO, i + '\t' + str(gc_hist[i])
			
		print >>sys.stderr, "writing R script ..."
		print >>RS, "pdf('GC_content.pdf')"
		print >>RS, 'gc=rep(c(' + ','.join([i for i in gc_hist.keys()]) + '),' + 'times=c(' + ','.join([str(i) for i in gc_hist.values()]) + '))'
		print >>RS, 'hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100
		#print >>RS, "lines(density(gc),col='red')"
		print >>RS ,"dev.off()"		
		#self.f.seek(0)
		
	def samDupRate(self,outfile=None,up_bound=500):
		'''Calculate reads's duplicate rates'''
		if outfile is None:
			outfile1 = self.fileName + ".seq.DupRate.xls"
			outfile2 = self.fileName + ".pos.DupRate.xls"
			outfile3 = self.fileName +".DupRate_plot.r"
		else:
			outfile1 = outfile + ".seq.DupRate.xls"
			outfile2 = outfile + ".pos.DupRate.xls"
			outfile3 = outfile +".DupRate_plot.r"
		SEQ=open(outfile1,'w')
		POS=open(outfile2,'w')
		RS=open(outfile3,'w')
		
		seqDup=collections.defaultdict(int)
		posDup=collections.defaultdict(int)
		
		seqDup_count=collections.defaultdict(int)
		posDup_count=collections.defaultdict(int)
		print >>sys.stderr, "reading sam file ... "
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			if (flagCode & 0x0004) == 1: 
				continue									#skip unmapped reads
			seqDup[field[9]] +=1			#key is read sequence
			
			#calculte duplicate read based on coordinates
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			chrom = field[2]
			chromStart = string.atoi(field[3])-1
			chromEnd=chromStart + sum(map(int,comb))
			blockSize = []
			for i in range(0,len(comb),2):
				blockSize.append(str(comb[i]))
			blockSizes = ','.join(blockSize)
			blockStart=[]
			for i in range(0,len(comb),2):
				blockStart.append(str(sum(comb[:i])))
			blockStarts = ','.join(blockStart)
			
			coord = chrom + ":" + str(chromStart) + "-" + str(chromEnd) + ":" + blockSizes + ":" + blockStarts
			posDup[coord] +=1
			
		print >>sys.stderr, "report duplicte rate based on sequence ..."
		print >>SEQ, "Occurrence\tUniqReadNumber"
		for i in seqDup.values():			#key is occurence, value is uniq reads number (based on seq)
			seqDup_count[i] +=1
		for k in sorted(seqDup_count.iterkeys()):	
			print >>SEQ, str(k) +'\t'+ str(seqDup_count[k])
		SEQ.close()
		
		print >>sys.stderr, "report duplicte rate based on mapping  ..."
		print >>POS, "Occurrence\tUniqReadNumber"
		for i in posDup.values():			#key is occurence, value is uniq reads number (based on coord)
			posDup_count[i] +=1
		for k in sorted(posDup_count.iterkeys()):	
			print >>POS, str(k) +'\t'+ str(posDup_count[k])
		POS.close()
		
		
		print >>sys.stderr, "generate R script ..."
		print >>RS, "pdf('duplicateRead.pdf')"
		print >>RS, "par(mar=c(5,4,4,5),las=0)"
		print >>RS, "seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
		print >>RS, "seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
		print >>RS, "pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.iterkeys()) ]) + ')'
		print >>RS, "pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.iterkeys()) ]) + ')'
		print >>RS, "plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Frequency',pch=4,cex=0.8,col='blue',xlim=c(1,%d),yaxt='n')" % up_bound
		print >>RS, "points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')"
		print >>RS, 'ym=floor(max(log10(pos_uniqRead)))'
		print >>RS, "legend(%d,ym,legend=c('Sequence-base','Mapping-base'),col=c('red','blue'),pch=c(4,20))" % max(up_bound-200,1)
		print >>RS, 'axis(side=2,at=0:ym,labels=0:ym)'
		print >>RS, 'axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead)),round(pos_uniqRead[2]*100/sum(pos_uniqRead)),round(pos_uniqRead[3]*100/sum(pos_uniqRead)),round(pos_uniqRead[4]*100/sum(pos_uniqRead))))'
		print >>RS, 'mtext(4, text = "Reads %", line = 2)'
		#self.f.seek(0)
		
	def getUniqMapRead(self,outfile=None):
		'''Extract uniquely mapped reads.'''
		if outfile is None:
			outfile = self.fileName + ".uniq.sam"
		FO=open(outfile,'w')
		Uniqcount=0
		print >>sys.stderr, "Writing uniquely mapped reads to\"",outfile,"\"... ",
		for line in self.f:
			hits=[]
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			if (flagCode & 0x0004) == 1: 
				continue			#skip unmapped reads
			#else:
				#print >>sys.stderr,line,
			if (ParseSAM._uniqueHit_pat.search(line)):
				print >>sys.stderr,line,
				Uniqcount +=1
				FO.write(line)
		FO.close()
		print >>sys.stderr, str(Uniqcount) + " reads were saved!\n",
		self.f.seek(0)

	def getWrongStrand(self,outfile=None):
		'''Extract pair-end reads mapped in incorrectly strand, such +/+ or -/-'''
		if outfile is None:
			outfile = self.fileName + ".wrongStrand.sam"
		FO=open(outfile,'w')
		wrongStrand=0
		print >>sys.stderr, "Writing incorrectly stranded reads to\"",outfile,"\"... ",
		for line in self.f:
			hits=[]
			if line.startswith('@'):continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			if (flagCode & 0x0004) != 0: 
				continue			#skipped if read itself is unmapped
			if (flagCode & 0x0008) !=0: 
				continue			#skipped if mate  is unmapped
			if (flagCode & 0x0001) == 0: 
				continue			#skip single end sequencing'
			if  ((flagCode & 0x0010) ==0) and ((flagCode & 0x0020)==0 ):
				FO.write(line)
				wrongStrand+=1
			if  ((flagCode & 0x0010) !=0) and ((flagCode & 0x0020)!=0 ):
				FO.write(line)
				wrongStrand+=1

		FO.close()
		print >>sys.stderr, str(wrongStrand) + " reads were saved!\n",
		self.f.seek(0)
		
	def filterSpliceRead(self,outfile=None,min_overhang=8,min_gap=50,max_gap=1000000):
		'''filter spiced mapped reads from sam file. min_overhang is used to determine the reliability of splice sites
		splice reads with overhang size <8 will also be reported if the same splice sites has been suported by
		at least 1 read with overhang size >8. Multiple spliced reads (belong to the same splice junction) 
		will always be reported. min_overhang, min_gap and max_gap only applied to one-time splice read'''
		
		if outfile is None:
			outfile = self.fileName + ".SR.sam"
			#outfile2 = self.fileName + ".SR.filter.sam"			
		splice_sites=collections.defaultdict(set)
		print >>sys.stderr, "\tDetermine splice sites with proper overhang, intron size ... ",
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			if not (ParseSAM._uniqueHit_pat.search(line)):	#skip non unique mapped read
				continue
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			map_st = int(field[3])
			chrom = field[2]
			if (flagCode & 0x0004) == 1:continue			#skip unmapped reads
			
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			if (len(comb)==1):	#skip non-spliced
				#print line,
				continue
			if (len(comb)>3):	#skip multiple spliced
				continue
			else:				#one-time spliced
				if (comb[1] < min_gap or comb[1] > max_gap):
					continue
				else:
					if (comb[0] >= min_overhang):
						splice_sites[chrom].add(map_st + comb[0])
					if (comb[2] >= min_overhang):
						splice_sites[chrom].add(map_st + comb[0] + comb[1])			
		self.f.seek(0)	
		print >>sys.stderr, "Done"
		
		
		FO=open(outfile,'w')
		#FO2=open(outfile2,'w')
		print >>sys.stderr, "\tExtracting splicing reads  ... ",
		total_SR =0
		extract_SR =0
		total_read =0
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			map_st = int(field[3])
			chrom = field[2]
			if (flagCode & 0x0004) == 1:continue			#skip unmapped reads
			total_read +=1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			if (len(comb)==1):	#skip non-spliced
				continue
			total_SR +=1
			
			
			if (len(comb)>3):	#multipel splice read. report directly
				FO.write(line)
				extract_SR +=1
			else:				#one-time splice read
				if (comb[1] < min_gap or comb[1] > max_gap):
					continue
				if (chrom in splice_sites) and ((map_st + comb[0]) in splice_sites[chrom]) and ((map_st + comb[0] + comb[1]) in splice_sites[chrom]):
					FO.write(line)
					#print line
					extract_SR +=1
				else:
					#FO2.write(line)
					continue
		print >>sys.stderr, "Done"
		print >>sys.stderr, "\tTotal mapped Read: " + str(total_read)
		print >>sys.stderr, "\tTotal Splicing Read: " + str(total_SR)
		print >>sys.stderr, "\Usable Splicing Read: " + str(extract_SR)
		FO.close()
		#FO2.close()
		self.f.seek(0)	

	def getSpliceRead(self,outfile=None):
		'''Extract spiced mapped reads from sam file'''
		
		if outfile is None:
			outfile = self.fileName + ".SR.sam"
		FO=open(outfile,'w')
		print >>sys.stderr, "\tExtract splicing reads without any filter ...",
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines
			field=line.rstrip('\n').split()	
			flagCode=string.atoi(field[1])
			if (flagCode & 0x0004) == 1:continue			#skip unmapped reads			
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			if (len(comb)>=3):
				FO.write(line)
		
		print >>sys.stderr, "Done"
		self.f.seek(0)	
		FO.close()

	def collapseSAM(self, outfile=None,collevel=10):
		'''At most collevel[default=10] identical reads will be retained in outputting SAM file
		The original SAM file must be sorted before hand. if not, using linux command like "sort -k3,3 -k4,4n myfile.sam >myfile.sorted.sam" '''
		if outfile is None:
			outfile = self.fileName + ".collapsed.sam"
		print >>sys.stderr, "Writing collapsed SAM file to\"",outfile,"\"... "
		FO=open(outfile,'w')		
		flag=""
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines	
			field=line.rstrip('\n').split()	
			if (string.atoi(field[1]) & 0x0004)!=0: continue	#skip unmapped line	
			id=field[2] + field[3] + field[5]
			if (id != flag):
				FO.write(line)
				flag=id
				skipTrigger=0
			else:
				skipTrigger+=1
				if skipTrigger < collevel:
					FO.write(line)
				else:continue
		FO.close()
		self.f.seek(0)

	def qualSAM(self,read_len,outfile=None):
		'''calculate phred quality score for each base in read (5->3)'''
		if outfile is None:
			outfile = self.fileName + ".qual.plot.r"
		else:
			outfile = outfile + ".qual.plot.r"
		FO=open(outfile,'w')
		print >>sys.stderr, "\tcalculating quality score ... "
		qual_min={}
		qual_max={}
		qual_sum={}
		total_read=0
		for i in range(0,read_len):
			qual_min.setdefault(i,1000)
			qual_max.setdefault(i,-1)
			qual_sum.setdefault(i,0.0)
		
		for line in self.f:
			if line[0] == '@':continue						#skip head lines	
			if ParseSAM._reExpr2.match(line):continue		#skip blank lines	
			field=line.rstrip('\n').split()
			#if (string.atoi(field[1]) & 0x0004)!=0: continue	#skip unmapped line
			
			if (len(field[10]) != read_len):
				continue
			if (string.atoi(field[1]) & 0x0010)==0:	#query map to +
				qual_str=field[10]
			else:
				qual_str=field[10][::-1]
			total_read +=1
			for i in range(0,read_len):
				#print ord(qual_str[i])-33,
				qual_sum[i] += ord(qual_str[i])-33
				if(qual_min[i] > (ord(qual_str[i])-33)):
					qual_min[i] = ord(qual_str[i])-33
				if(qual_max[i] < (ord(qual_str[i])-33)):
					qual_max[i] = ord(qual_str[i])-33
			#print '\n',
		min_qualities = [str(qual_min[i]) for i in range(0,read_len)]
		max_qualities =[str(qual_max[i]) for i in range(0,read_len)]
		avg_qualities = [str(qual_sum[i]/total_read) for i in range(0,read_len)]
		nt_pos = [str(i) for i in range(0,read_len)]
		print >>FO, "nt_pos=c(" + ','.join(nt_pos)  + ')'
		print >>FO, "max_qual=c(" + ','.join(max_qualities)  + ')'
		print >>FO, "min_qual=c(" + ','.join(min_qualities)  + ')'
		print >>FO, "avg_qual=c(" + ','.join(avg_qualities)  + ')'
		print >>FO, "pdf('phred_qual.pdf')"
		print >>FO, "plot(nt_pos,avg_qual, xlab=\"Nucleotide Position (5'->3')\", ylab='Phred Quality',ylim=c(0,97),lwd=2,type='s')"
		print >>FO, 'lines(nt_pos,max_qual,type="s",lwd=2,col="red")'
		print >>FO, 'lines(nt_pos,min_qual,type="s",lwd=2,col="blue")'
		print >>FO, 'legend(0,100,legend=c("Max","Average","Min"),col=c("red","black","blue"),lwd=2)'
		print >>FO, 'dev.off()'
		#for i in range(0,read_len):
		#	print >>sys.stderr, str(i) + '\t' + str(qual_max[i]) + '\t' + str(qual_min[i]) + '\t' + str(qual_sum[i]/total_read)
		#self.f.seek(0)


	def samToBinnedArray(self):
		"""Convert SAM file to BinnedArray."""
		
		lines=0
		for line in self.f:
			if line.startswith('@'):continue						#skip head lines
			if ParseSAM._reExpr2.match(line):continue				#skip blank lines
			field=line.rstrip().split()	
			if (string.atoi(field[1]) & 0x0004)!=0: continue		#skip unmapped line		
			txStart=string.atoi(field[3])
			#if (string.atoi(field[1]) & 0x0010 != 0):
			#	strand='-'
			#else:
			#	strand='+'
			lines +=1
			scores={}
			chrom = field[2]
			
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			if not chrom in scores:scores[chrom] = BinnedArray()
			
			for i in range(0,len(comb),2):
				for pos in range(txStart + sum(comb[:i]), txStart + sum(comb[:i]) + comb[i]):
					tmp = scores[chrom][pos]
					if isNaN(tmp):				
						scores[chrom][pos] =1
					else:
						scores[chrom][pos] +=1			
			if lines % 10000 == 0: print >>sys.stderr, "%i lines loaded \r" % lines
		return scores
		self.f.seek(0)

class QCSAM:
	'''Perform basic quality control. Useful for RNA-seq experiment'''
	

	def __init__(self,samFile):
		'''constructor'''
		if samFile == '-':
			self.fileName = "STDIN"
			self.f = sys.stdin
		else:
			self.fileName = os.path.basename(samFile)
			self.f=open(samFile,'r')

		
	def distribSAM(self,refbed,outfile=None):
		'''calculate reads distribution over genome features (Exon reads, Inton reads, Intergenic 
		reads, spliced reads). A bed file representing the gene model (i.e. refseq) must be provided
		Two bed format files will be generated: outfile_exon.count.bed and outfile_intron.count.bed.
		The 5th column is number of reads fallen into the region defined by the first 3 columns'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if outfile is None:
			exon_count = self.fileName + "_exon.count.bed"
			intron_count = self.fileName + "_intron.count.bed"
			rscript=self.fileName + ".piechart.r"
			rpdf=self.fileName + ".piechart.pdf"
		else:
			exon_count = outfile + "_exon.count.bed"
			intron_count = outfile +  "_intron.count.bed"
			rscript=outfile + ".piechart.r"
			rpdf=outfile + ".piechart.pdf"
		
		EXON_OUT = open(exon_count,'w')
		INTRON_OUT =open(intron_count,'w')
		R_OUT = open(rscript,'w')
		
		ranges={}
		intronReads=0
		exonReads=0
		intergenicReads=0
		totalReads=0
		splicedReads=0
		
		#read SAM 
		print >>sys.stderr, "reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue		#skip unmap reads
			totalReads +=1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			if( len(comb)>1):
				splicedReads +=1	
				continue
			else:
				chrom=fields[2].upper()
				#st=int(fields[3])-1
				#end= st +len(fields[9])
				mid = int(fields[3]) + int(len(fields[9])/2)
				if chrom not in ranges:
					ranges[chrom] = Intersecter()
				else:
					ranges[chrom].add_interval( Interval( mid, mid ) )
    			
		self.f.seek(0)
		print >>sys.stderr, "Done"
		
		#read refbed file
		print >>sys.stderr, "Assign reads to "+ refbed + '...',
		for line in open(refbed,'r'):
			try:
				if line.startswith('#'):continue
				if line.startswith('track'):continue
				if line.startswith('browser'):continue   
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
				intron_starts = exon_ends[:-1]
				intron_ends=exon_starts[1:]
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue

				# assign reads to intron				
			if(strand == '-'):
				intronNum=len(intron_starts)
				exonNum=len(exon_starts)
				for st,end in zip(intron_starts,intron_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						intronReads += hits
						INTRON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\n')
						intronNum -= 1
							
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						exonReads += hits
						EXON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\n')
						exonNum -= 1
			elif(strand == '+'):
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						intronReads += hits
						INTRON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\n')
						intronNum += 1    
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						exonReads += hits
						EXON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\n')
						exonNum += 1		
		intergenicReads=totalReads-exonReads-intronReads-splicedReads		
		print >>sys.stderr, "Done." + '\n'
		print >>sys.stderr, "Total reads:\t" + str(totalReads)
		print >>sys.stderr, "Exonic reads:\t" + str(exonReads) 
		print >>sys.stderr, "Intronic reads:\t" + str(intronReads) 
		print >>sys.stderr, "Splicing reads:\t" + str(splicedReads)
		print >>sys.stderr, "Intergenic reads:\t" + str(intergenicReads)
		
		print >>sys.stderr,"writing R script ...",
		totalReads=float(totalReads)
		print >>R_OUT, "pdf('%s')" % rpdf
		print >>R_OUT, "dat=c(%d,%d,%d,%d)" % (exonReads,splicedReads,intronReads,intergenicReads)
		print >>R_OUT, "lb=c('exon(%.2f)','junction(%.2f)','intron(%.2f)','intergenic(%.2f)')" % (exonReads/totalReads,splicedReads/totalReads,intronReads/totalReads,intergenicReads/totalReads)
		print >>R_OUT, "pie(dat,labels=lb,col=rainbow(4),clockwise=TRUE,main='Total reads = %d')" % int(totalReads)
		print >>R_OUT, "dev.off()"
		print >>sys.stderr, "Done."
	
	
	def coverageGeneBody(self,refbed,outfile=None):
		'''Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divied
		into 100 regsions'''
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if outfile is None:
			outfile1 = self.fileName + ".geneBodyCoverage_plot.r"
			outfile2 =  self.fileName + ".geneBodyCoverage.txt"
		else:
			outfile1 = outfile + ".geneBodyCoverage_plot.r"
			outfile2 = outfile + ".geneBodyCoverage.txt"
		OUT1 = open(outfile1,'w')
		OUT2 = open(outfile2,'w')

		ranges={}
		totalReads=0
		fragment_num=0		#splice reads will counted twice
		rpkm={}
		
		#read SAM 
		print >>sys.stderr, "reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue		#skip unmap reads
			totalReads +=1
			
			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			fragment_num += (len(comb) +1)/2
			blockStart=[]
			blockSize=[]
			
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			for st,size in zip(blockStart,blockSize):
				if chrom not in ranges:
					ranges[chrom] = Intersecter()
				else:
					ranges[chrom].add_interval( Interval( st, st+size ) )
		print >>sys.stderr, "Done"   	

		print >>sys.stderr, "calculating coverage over gene body ..."
		coverage=collections.defaultdict(int)
		flag=0
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5]
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue
			gene_all_base=[]
			percentile_base=[]
			mRNA_len =0
			flag=0
			for st,end in zip(exon_starts,exon_ends):
				gene_all_base.extend(range(st+1,end+1))		#0-based coordinates on genome
				mRNA_len = len(gene_all_base)
				if mRNA_len <100:
					flag=1
					break
			if flag==1: continue
			if strand == '-':
				gene_all_base.sort(reverse=True)			#deal with gene on minus stand
			else:
				gene_all_base.sort(reverse=False)
			percentile_base = mystat.percentile_list (gene_all_base)	#get 101 points from each gene's coordinates
			
			for i in range(0,len(percentile_base)):
				if chrom in ranges:
					coverage[i] += len(ranges[chrom].find(percentile_base[i], percentile_base[i]+1))
		x_coord=[]
		y_coord=[]
		print >>OUT2, "Total reads: " + str(totalReads)
		print >>OUT2, "Fragment number: " + str(fragment_num)
		print >>OUT2, "percentile\tcount"
		for i in coverage:
			x_coord.append(str(i))
			y_coord.append(str(coverage[i]))
			print >>OUT2, str(i) + '\t' + str(coverage[i])
		print >>OUT1, "pdf('geneBody_coverage.pdf')"
		print >>OUT1, "x=0:100"
		print >>OUT1, "y=c(" + ','.join(y_coord) + ')'
		print >>OUT1, "plot(x,y,xlab=\"percentile of gene body (5'->3')\",ylab='read number',type='s')"
		print >>OUT1, "dev.off()"
			
	def calculateRPKM(self,refbed,outfile=None):
		'''calculate RPKM values for each gene in refbed. Only uniquely aligned reads are used. 
		Spilced reads are split. output raw read connt and eRPKM (eRPKM = exon Represented times Per Kb 
		exon per Million mapped reads) for each exon, intron and mRNA'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if outfile is None:
			rpkm_file = self.fileName + ".rpkm.xls"
		else:
			rpkm_file = outfile + ".rpkm.xls"
		RPKM_OUT=open(rpkm_file,'w')
		
		ranges={}
		totalReads=0
		cUR=0
		sR=0
		multiMapReads=0
		rpkm={}
		
		#read SAM 
		print >>sys.stderr, "reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue		#skip unmap reads
			totalReads +=1
			if not ParseSAM._uniqueHit_pat.search(line):		#skip multiple mapped reads
				multiMapReads +=1
				continue

			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			cUR += (len(comb) +1)/2
			if(len(comb)>1):
				sR+=1
			blockStart=[]
			blockSize=[]
			
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			for st,size in zip(blockStart,blockSize):
				mid = int(st) + (size/2)
				if chrom not in ranges:
					ranges[chrom] = Intersecter()
				else:
					ranges[chrom].add_interval( Interval( mid, mid ) )
    			
		self.f.seek(0)
		print >>sys.stderr, "Done"
		print >>RPKM_OUT, "Total mapped reads (TR): " + str(totalReads) 
		print >>RPKM_OUT, "Multiple mapped reads (MR): " + str(multiMapReads)
		print >>RPKM_OUT, "Uniquely mapped reads (UR): " + str(totalReads - multiMapReads)
		print >>RPKM_OUT, "Spliced  mapped reads (SR): " + str(sR)
		print >>RPKM_OUT, "Corrected uniquely mapped reads (cUR): " + str(cUR)
		if totalReads ==0:
			sys.exit(1)
		
		#read refbed file
		print >>sys.stderr, "Assign reads to "+ refbed + '...',
		for line in open(refbed,'r'):
			try:
				if line.startswith('#'):continue
				if line.startswith('track'):continue
				if line.startswith('browser'):continue   
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
				exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
				intron_starts = exon_ends[:-1]
				intron_ends=exon_starts[1:]
				key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue

			# assign reads to intron				
			mRNA_count=0
			mRNA_len=sum(exon_sizes)
			if(strand == '-'):
				intronNum=len(intron_starts)
				exonNum=len(exon_starts)
								
				for st,end in zip(intron_starts,intron_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						intronNum -= 1
							
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))					
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						exonNum -= 1
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
					rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
					rpkm[key] = 0
			elif(strand == '+'):
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						intronNum += 1    
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						exonNum += 1		
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
					rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
					rpkm[key] = 0
		print >>sys.stderr, "Done"
		return rpkm
		self.f.seek(0)

	def calculateRPKM2(self,refbed,outfile=None):
		'''calculate RPKM values for each gene in refbed. Only uniquely aligned reads are used. 
		Spilced reads are split. output raw read connt and eRPKM (eRPKM = exon Represented times Per Kb 
		exon per Million mapped reads) for each exon, intron and mRNA
		NOTE: intronic reads are not counted as part of total reads'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if outfile is None:
			rpkm_file = self.fileName + ".rpkm.xls"
		else:
			rpkm_file = outfile + ".rpkm.xls"
		RPKM_OUT=open(rpkm_file,'w')
		
		ranges={}
		exon_ranges={}
		totalReads=0
		#intronic=0
		cUR=0
		sR=0
		multiMapReads=0
		rpkm={}
		
		#read gene model file, the purpose is to remove intronic reads
		print >>sys.stderr, "Reading reference gene model "+ refbed + '...'
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue
 
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue	

			for st,end in zip(exon_starts,exon_ends):				
				if chrom not in exon_ranges:
					exon_ranges[chrom] = Intersecter()
				else:
					exon_ranges[chrom].add_interval( Interval( st, end ) )		

		#read SAM 
		print >>sys.stderr, "reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue		#skip unmap reads
			totalReads +=1
			if not ParseSAM._uniqueHit_pat.search(line):		#skip multiple mapped reads
				multiMapReads +=1
				continue

			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			#cUR += (len(comb) +1)/2
			if(len(comb)>1):
				sR+=1
			blockStart=[]
			blockSize=[]
			
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			#build bitset only for exonic reads
			for st,size in zip(blockStart,blockSize):
				if (chrom in exon_ranges) and (len(exon_ranges[chrom].find(st,st+size)) >0):	#if we found this fragment is overlapped with exon
					cUR += 1
					mid = int(st) + (size/2)
					if chrom not in ranges:
						ranges[chrom] = Intersecter()
					else:
						ranges[chrom].add_interval( Interval( mid, mid ) )
    			else:																			#if this framgnet is intronic, skip it.
    				#intronic +=1
    				continue	
		self.f.seek(0)
		print >>sys.stderr, "Done"
		print >>RPKM_OUT, "Total mapped reads (TR): " + str(totalReads) 
		print >>RPKM_OUT, "Multiple mapped reads (MR): " + str(multiMapReads)
		print >>RPKM_OUT, "Uniquely mapped reads (UR): " + str(totalReads - multiMapReads)
		print >>RPKM_OUT, "Spliced  mapped reads (SR): " + str(sR)
		print >>RPKM_OUT, "Corrected uniquely mapped reads (cUR, non-intronic fragments): " + str(cUR)
		#print >>RPKM_OUT, "Intronic Fragments (IF): " + str(intronic)
		if totalReads ==0:
			sys.exit(1)
		
		#read refbed file
		print >>sys.stderr, "Assign reads to "+ refbed + '...',
		for line in open(refbed,'r'):
			try:
				if line.startswith('#'):continue
				if line.startswith('track'):continue
				if line.startswith('browser'):continue   
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
				exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
				intron_starts = exon_ends[:-1]
				intron_ends=exon_starts[1:]
				key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue

			# assign reads to intron				
			mRNA_count=0
			mRNA_len=sum(exon_sizes)
			if(strand == '-'):
				intronNum=len(intron_starts)
				exonNum=len(exon_starts)
								
				for st,end in zip(intron_starts,intron_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						intronNum -= 1
							
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))					
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						exonNum -= 1
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
					rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
					rpkm[key] = 0
			elif(strand == '+'):
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						intronNum += 1    
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						hits= len(ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
						exonNum += 1		
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
					rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
					rpkm[key] = 0
		print >>sys.stderr, "Done"
		return rpkm
		self.f.seek(0)

	def filterKnownReads(self,refbed,outfile=None):
		'''Compare SAM files with reference gene model, all reads mapped to gene model will be filted
		out. The remainning unknown reads will be writtern to a new SAM file'''
		
		totalReads=0	#total mapped reads
		unknownReads=0
		ranges={}
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		
		if outfile is None:
			out_file = self.fileName + ".unknownReads.SAM"
		else:
			out_file = outfile + ".unknownReads.SAM"	
		OUT=open(out_file,'w')
		
		print >>sys.stderr, "Reading reference gene model "+ refbed + '...'
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue
 
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue	

			for st,end in zip(exon_starts,exon_ends):				
				if chrom not in ranges:
					ranges[chrom] = Intersecter()
				else:
					ranges[chrom].add_interval( Interval( st, end ) )		

		print >>sys.stderr, "Processing SAM file "+ self.fileName + '...'
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue			#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):	#skip multiple mapped reads
				continue
			
			blockStart=[]
			blockSize=[]
			totalReads +=1
			
			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']			
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])

			for st,size in zip(blockStart,blockSize):
				if (chrom in ranges) and (len(ranges[chrom].find(st,st+size)) >0):	#if we found this read is overlapped with known gene
					break
			else:
				OUT.write(line)
				unknownReads +=1
		OUT.close()
		print >>sys.stderr, "Total reads mapped to genome: " + str(totalReads)
		print >>sys.stderr, "Total reads not overlapped with any exon: " + str(unknownReads)
		self.f.seek(0)

	def genomicFragSize(self,outfile=None,low_bound=0,up_bound=1000,step=10):
		'''estimate the genomic fragment size of mRNA experiment. fragment size = insert_size + 2 x read_length'''
		
		
		if outfile is None:
			out_file1 = self.fileName + ".fragSize.txt"
			out_file2 = self.fileName + ".fragSize.Freq.txt"
			out_file3 = self.fileName + ".fragSize_plot.r"
		else:
			out_file1 = outfile + ".fragSize.txt"	
			out_file2 = outfile + ".fragSize.Freq.txt"
			out_file3 = outfile + ".fragSize_plot.r"
		
		FO=open(out_file1,'w')
		FQ=open(out_file2,'w')
		RS=open(out_file3,'w')
		
		chrom="chr100"		#this is the fake chromosome
		ranges={}
		ranges[chrom]=Intersecter()
		
		window_left_bound = range(low_bound,up_bound,step)
 		frag_size=0

 		pair_num=0.0
 		ultra_low=0.0
 		ultra_high=0.0
 		size=[]
 		counts=[]
 		count=0
		print >>sys.stderr, "Reading SAM file "+ self.fileName + ' ... ',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			#if fields[0] in pairRead_info:
			#	continue
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0001) ==0:
				print >>sys.stderr,"NOT pair-end sequencing"
				sys.exit(1)
			if (flagCode & 0x0004) != 0: continue			#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):	#skip multiple mapped reads
				continue
			if (flagCode & 0x0008 !=0):						#skip single-end mapped reads
				continue
			if (fields[7] =='0'):
				continue
			if (int(fields[3]) > int(fields[7])):			#left < right
				continue
			pair_num +=1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			read_len = len(fields[9])
			if (len(comb)==1):		# this read is NOT spliced
				frag_size = (int(fields[7]) - int(fields[3]) +1) + read_len
			elif (len(comb) >1):	# this read is spliced
				frag_size = (int(fields[7]) - int(fields[3]) +1) + read_len - sum(comb[1::2])
			FO.write(fields[0] + '\t' + str(frag_size) + '\n')
			if frag_size <= low_bound:
				ultra_low+=1
				continue
			elif frag_size > up_bound:
				ultra_high +=1
				continue
			ranges[chrom].add_interval( Interval( frag_size-1, frag_size ) )
		print >>sys.stderr, "Done"
		
		if pair_num==0:
			print >>sys.stderr, "Cannot find paired reads"
			sys.exit(0)
		print >>FQ, "Total paired read " + str(pair_num)
		print >>FQ, "<=" + str(low_bound) + "\t"+ str(ultra_low)
		for st in window_left_bound:
			size.append(str(st + step/2))
			count = str(len(ranges[chrom].find(st,st + step)))
			counts.append(count)
			print >>FQ, str(st) + '\t' + str(st+step) +'\t' + count		
		print >>FQ, ">" + str(up_bound) + "\t"+ str(ultra_high)
		
		print >>RS, "pdf('gFragSize.pdf')"
		print >>RS, "par(mfrow=c(2,1),cex.main=0.8,cex.lab=0.8,cex.axis=0.8,mar=c(4,4,4,1))"
		print >>RS, 'pie(c(%d,%d,%d),col=rainbow(3),cex=0.5,radius=1,main="Total %d fragments",labels=c("fraSize <= %d\\n(%4.2f%%)","fragSize > %d\\n(%4.2f%%)","%d < fragSize <= %d\\n(%4.2f%%)"), density=rep(80,80,80),angle=c(90,140,170))' % (ultra_low, ultra_high, pair_num -ultra_low -ultra_high, pair_num, low_bound, ultra_low*100/pair_num, up_bound, ultra_high*100/pair_num, low_bound, up_bound, 100-ultra_low*100/pair_num - ultra_high*100/pair_num)
		print >>RS, 'fragsize=rep(c(' + ','.join(size) + '),' + 'times=c(' + ','.join(counts) + '))'
		print >>RS, 'frag_sd = round(sd(fragsize))'
		print >>RS, 'frag_mean = round(mean(fragsize))'
		print >>RS, 'hist(fragsize,probability=T,breaks=%d,xlab="Fragment size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(window_left_bound)
		print >>RS, "lines(density(fragsize,bw=%d),col='red')" % (2*step)
		print >>RS ,"dev.off()"
		FO.close()
		FQ.close()
		RS.close()
		#self.f.seek(0)
		
		
	def saturation_RPKM(self,refbed,outfile=None,sample_start=5,sample_step=5,sample_end=100):
		'''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if outfile is None:
			rpkm_file = self.fileName + ".eRPKM.xls"
			raw_file = self.fileName + ".rawCount.xls"
		else:
			rpkm_file = outfile + ".eRPKM.xls"
			raw_file = outfile + ".rawCount.xls"
		
		RPKM_OUT = open(rpkm_file,'w')
		RAW_OUT = open(raw_file ,'w')
		
		ranges={}
		totalReads=0
		cUR_num = 0	#number
		block_list=[]	#non-spliced read AS IS, splicing reads were counted multiple times
				
		#read SAM 
		my_pat = re.compile(r'NH:i:(\d+)\b')
		NH_tag=0
		print >>sys.stderr, "Reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue		#skip unmap reads
			totalReads +=1
			hitNum =[int(i) for i in my_pat.findall(line)]
			if len(hitNum) ==0:
				NH_tag=1								#cannot determine uniqness without NH tag
			elif len(hitNum) ==1:
				if int(hitNum[0])>1: continue				#skip multiple mapped reads
			else:
				print >>sys.stderr, "More than 1 NH tag found within a single line. Incorrect SAM format!"
				sys.exit(1)

			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			cUR_num += (len(comb) +1)/2
			blockStart=[]
			blockSize=[]
			
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			for st,size in zip(blockStart,blockSize):
				mid = int(st) + (size/2)
				block_list.append(chrom + ":" + str(mid))
		
		if NH_tag==1:
			print >>sys.stderr, "Warn: NO NH tag found. Cannot determine uniqueness of alignment. All alignments will be used"
		print >>sys.stderr, "Done"
		
		print >>sys.stderr, "shuffling alignments ...",
		random.shuffle(block_list)
		print >>sys.stderr, "Done"
		
		
		ranges={}
		sample_size=0
		frag_total = cUR_num
		RPKM_table=collections.defaultdict(list)
		rawCount_table=collections.defaultdict(list)
		RPKM_head=['chr','start','end','name','score','strand']

		tmp=range(sample_start,sample_end,sample_step)
		tmp.append(100)
		#=========================sampling uniquely mapped reads from population
		for pertl in tmp:	#[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
			index_st = int(frag_total * (pertl-sample_step)/100.0)
			index_end = int(frag_total * pertl/100.0)
			if index_st < 0: index_st = 0
			sample_size += index_end -index_st
			
			RPKM_head.append(str(pertl) + '%')
			print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(sample_size) + ") fragments ...",
			for i in range(index_st, index_end):
				(chr,coord) = block_list[i].split(':')
				if chr not in ranges:
					ranges[chr] = Intersecter()
				else:
					ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )				
			#========================= calculating RPKM based on sub-population
			#print >>sys.stderr, "assign reads to "+ refbed + '...',
			for line in open(refbed,'r'):
				try:
					if line.startswith(('#','track','browser')):continue  
            		# Parse fields from gene tabls
					fields = line.split()
					chrom     = fields[0].upper()
					tx_start  = int( fields[1] )
					tx_end    = int( fields[2] )
					geneName      = fields[3]
					strand    = fields[5].replace(" ","_")
					exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
					exon_starts = map((lambda x: x + tx_start ), exon_starts)
					exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
					exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
					exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
					key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
				except:
					print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
					continue
				mRNA_count=0	#we need to initializ it to 0 for each gene
				mRNA_len=sum(exon_sizes)
				for st,end in zip(exon_starts,exon_ends):
					if chrom in ranges:
						mRNA_count += len(ranges[chrom].find(st,end))		
				if mRNA_len ==0:
					print >>sys.stderr, geneName + " has 0 nucleotides. Exit!"
					sys.exit(1)
				if sample_size == 0:
					print >>sys.stderr, "Too few reads to sample. Exit!"
					sys.exit(1)
				mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * sample_size)
				RPKM_table[key].append(str(mRNA_RPKM))
				rawCount_table[key].append(str(mRNA_count))
			print >>sys.stderr, "Done"

		#self.f.seek(0)
		print >>RPKM_OUT, '\t'.join(RPKM_head)
		print >>RAW_OUT, '\t'.join(RPKM_head)
		for key in RPKM_table:
			print >>RPKM_OUT, key + '\t',
			print >>RPKM_OUT, '\t'.join(RPKM_table[key])
			print >>RAW_OUT, key + '\t',
			print >>RAW_OUT, '\t'.join(rawCount_table[key])			

	def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1):
		'''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction'''
		
		if outfile is None:
			out_file = self.fileName + ".junctionSaturation_plot.r"
		else:
			out_file = outfile + ".junctionSaturation_plot.r"
		if refgene is None:
			print >>sys.stderr, "You must provide reference gene model in bed format."
			sys.exit(1)
		
		OUT = open(out_file,'w')


		#reading reference gene 
		knownSpliceSites= set()
		print >>sys.stderr, "reading reference bed file: ",refgene, " ... ",
		for line in open(refgene,'r'):
			if line.startswith(('#','track','browser')):continue  
			fields = line.split()
			if(len(fields)<12):
				print >>sys.stderr, "Invalid bed line (skipped):",line,
				continue
			chrom     = fields[0].upper()
			tx_start = int( fields[1] )
			tx_end   = int( fields[2] )
			if int(fields[9] ==1):
				continue    	
			
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			for st,end in zip (intron_start, intron_end):
				knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))
		print >>sys.stderr,"Done! Total "+str(len(knownSpliceSites)) + " known splicing sites"


		#read SAM file
		samSpliceSites=[]
		intron_start=[]
		intron_end=[]
		uniqSpliceSites=collections.defaultdict(int)
		print >>sys.stderr, "Reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			
			if (flagCode & 0x0004) != 0: continue				#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):		#skip multiple mapped reads
				continue
			if (len(ParseSAM._splicedHit_pat.findall(fields[5]))==1):	#skip non-spilced reads
				continue

			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			blockStart=[]
			blockSize=[]
			blockEnd=[]		
			#if intron size < min_intron, skip. 
			flag=0
			for i in range(1,len(comb),2):
				if comb[i] < min_intron:
					flag=1
					break
			if flag ==1:
				continue		
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			for st,size in zip(blockStart,blockSize):
				end = st + size
				blockEnd.append(end)
			intron_st = blockEnd[:-1]
			intron_end = blockStart[1:]
			
			for st,end in zip(intron_st, intron_end):
				samSpliceSites.append(chrom + ":" + str(st) + "-" + str(end))			
		#self.f.seek(0)
		print >>sys.stderr, "Done"
		


		print >>sys.stderr, "shuffling alignments ...",
		random.shuffle(samSpliceSites)
		print >>sys.stderr, "Done"
				
		#resampling
		SR_num = len(samSpliceSites)
		sample_size=0
		knownSpliceSites_num = 0
		known_junc=[]
		all_junc=[]
		#=========================sampling uniquely mapped reads from population
		tmp=range(sample_start,sample_end,sample_step)
		tmp.append(100)
		for pertl in tmp:	#[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
			knownSpliceSites_num = 0
			index_st = int(SR_num * ((pertl - sample_step)/100.0))
			index_end = int(SR_num * (pertl/100.0))
			if index_st < 0: index_st = 0
			sample_size += index_end -index_st
			
			print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(sample_size) + ") unique splicing alignments ...",
			for i in range(index_st, index_end):
				uniqSpliceSites[samSpliceSites[i]] +=1	
			all_junc.append(str(len(uniqSpliceSites.keys())))
			for sj in uniqSpliceSites:
				if sj in knownSpliceSites and uniqSpliceSites[sj] >= recur:
					knownSpliceSites_num +=1
			print >>sys.stderr, str(knownSpliceSites_num) + " known splicing junctions"
			known_junc.append(str(knownSpliceSites_num))
			
		#for j in uniq_SJ:
			#print >>OUT, j + "\t" + str(uniq_SJ[j])
		print >>OUT, "pdf('junction_saturation.pdf')"
		print >>OUT, "x=c(" + ','.join([str(i) for i in tmp]) + ')'
		print >>OUT, "y=c(" + ','.join(known_junc) + ')'
		print >>OUT, "z=c(" + ','.join(all_junc) + ')'
		print >>OUT, "plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(%d,%d))" % (int(int(known_junc[0])/1000), int(int(all_junc[-1])/1000))
		print >>OUT, "points(x,y/1000,type='o',col='red')"
		print >>OUT, 'legend(5,%d, legend=c("All detected junction","Annotated junction"),col=c("blue","red"),lwd=1,pch=1)' % int(int(all_junc[-1])/1000)
		print >>OUT, "dev.off()"

	
	def annotate_junction(self,refgene,outfile=None,min_intron=50):
		'''Annotate splicing junctions in SAM file. Note that a (long) read might have multiple splicing
		events  (splice multiple times), and the same splicing events can be consolidated into a single
		junction'''
		
		if outfile is None:
			out_file = self.fileName + ".junction.xls"
			out_file2 = self.fileName + ".junction_plot.r"
		else:
			out_file = outfile + ".junction.xls"
			out_file2 = outfile + ".junction_plot.r"
		if refgene is None:
			print >>sys.stderr, "You must provide reference gene model in bed format."
			sys.exit(1)
		OUT = open(out_file,'w')
		ROUT = open(out_file2,'w')
		
		#reading reference gene model
		refIntronStarts=collections.defaultdict(dict)
		refIntronEnds=collections.defaultdict(dict)	
		total_junc =0
		novel35_junc =0
		novel3or5_junc =0
		known_junc =0
		splicing_events=collections.defaultdict(int)	
		
		print >>sys.stderr, "\treading reference bed file: ",refgene, " ... ",
		for line in open(refgene,'r'):
			if line.startswith(('#','track','browser')):continue  
           	# Parse fields from gene tabls
			fields = line.split()
			if(len(fields)<12):
				print >>sys.stderr, "Invalid bed line (skipped):",line,
				continue
			chrom     = fields[0].upper()
			tx_start = int( fields[1] )
			tx_end   = int( fields[2] )
			if int(fields[9] ==1):
				continue    	
			
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			for i_st,i_end in zip (intron_start, intron_end):
				refIntronStarts[chrom][i_st] =i_st
				refIntronEnds[chrom][i_end] =i_end			
		print >>sys.stderr,"Done"
		
		#reading input SAM file
		print >>sys.stderr, "\tProcessing "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			
			if (flagCode & 0x0004) != 0: continue				#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):		#skip multiple mapped reads
				continue
			if (len(ParseSAM._splicedHit_pat.findall(fields[5]))==1):	#skip non-spilced reads
				continue

			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			blockStart=[]
			blockSize=[]
			blockEnd=[]		
			#if intron size < min_intron, skip. 
			flag=0
			for i in range(1,len(comb),2):
				if comb[i] < min_intron:
					flag=1
					break
			if flag ==1:
				continue		
			
			total_junc += (len(comb) -1)/2
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			for st,size in zip(blockStart,blockSize):
				end = st + size
				blockEnd.append(end)
			intron_st = blockEnd[:-1]
			intron_end = blockStart[1:]
			for i_st,i_end in zip(intron_st, intron_end):
				splicing_events[chrom + ":" + str(i_st) + ":" + str(i_end)] += 1
				if (refIntronStarts[chrom].has_key(i_st) and refIntronEnds[chrom].has_key(i_end)):
					known_junc +=1																		#known both
				elif (not refIntronStarts[chrom].has_key(i_st) and not refIntronEnds[chrom].has_key(i_end)):
					novel35_junc +=1																
				else:
					novel3or5_junc +=1
		#self.f.seek(0)
		print >>sys.stderr, "Done"
		
		print >>ROUT, 'pdf("splicing_events_pie.pdf")'
		print >>ROUT, "events=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc)])+ ')'
		print >>ROUT, 'pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing events",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc))
		print >>ROUT, "dev.off()"
		
		print >>sys.stderr, "\n==================================================================="
		print >>sys.stderr, "Total splicing  Events:\t" + str(total_junc)
		print >>sys.stderr, "Known Splicing Events:\t" + str(known_junc)
		print >>sys.stderr, "Partial Novel Splicing Events:\t" + str(novel3or5_junc)
		print >>sys.stderr, "Novel Splicing Events:\t" + str(novel35_junc)
		
		#reset variables
		total_junc =0
		novel35_junc =0
		novel3or5_junc =0
		known_junc =0
		
		print >>OUT, "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation"
		for i in splicing_events:
			total_junc += 1
			(chrom, i_st, i_end) = i.split(":")
			print >>OUT, '\t'.join([chrom.replace("CHR","chr"),i_st,i_end]) + '\t' + str(splicing_events[i]) + '\t',
			i_st = int(i_st)
			i_end = int(i_end)
			if (refIntronStarts[chrom].has_key(i_st) and refIntronEnds[chrom].has_key(i_end)):
				print >>OUT, "annotated"
				known_junc +=1
			elif (not refIntronStarts[chrom].has_key(i_st) and not refIntronEnds[chrom].has_key(i_end)):
				print >>OUT, 'complete_novel'
				novel35_junc +=1
			else:
				print >>OUT, 'partial_novel'
				novel3or5_junc +=1
		
		
		print >>sys.stderr, "\nTotal splicing  Junctions:\t" + str(total_junc)
		print >>sys.stderr, "Known Splicing Junctions:\t" + str(known_junc)
		print >>sys.stderr, "Partial Novel Splicing Junctions:\t" + str(novel3or5_junc)
		print >>sys.stderr, "Novel Splicing Junctions:\t" + str(novel35_junc)
		print >>sys.stderr, "\n==================================================================="
		
		print >>ROUT, 'pdf("splicing_junction_pie.pdf")'
		print >>ROUT, "junction=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc,)])+ ')'
		print >>ROUT, 'pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing junctions",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc))
		print >>ROUT, "dev.off()"
		#print >>ROUT, "mat=matrix(c(events,junction),byrow=T,ncol=3)"
		#print >>ROUT, 'barplot(mat,beside=T,ylim=c(0,100),names=c("known","partial\nnovel","complete\nnovel"),legend.text=c("splicing events","splicing junction"),ylab="Percent")'

	def mRNA_RPKM(self,refbed,outfile=None):
		'''calculate mRNA's RPKM value'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if outfile is None:
			rpkm_file = self.fileName + ".RPKM.xls"
		else:
			rpkm_file = outfile + ".RPKM.xls"		
		RPKM_OUT = open(rpkm_file,'w')                  
		
		ranges={}
		totalReads=0
		cUR_num = 0	#number

		RPKM_table={}
		rawCount_table={}
		mRNAlen_table={}
		RPKM_head=['chr','start','end','name','score','strand','length','rawCount','RPKM']		
	
		#read SAM 
		print >>sys.stderr, "Reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue		#skip unmap reads
			totalReads +=1
			if not ParseSAM._uniqueHit_pat.search(line):		#skip multiple mapped reads
				continue

			chrom = fields[2].upper()
			chromStart = string.atoi(fields[3])-1
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			cUR_num += (len(comb) +1)/2
			blockStart=[]
			blockSize=[]
			
			for i in range(0,len(comb),2):
				blockStart.append(chromStart + sum(comb[:i]) )
				
			for i in range(0,len(comb),2):
				blockSize.append(comb[i])
			
			for st,size in zip(blockStart,blockSize):
				mid = int(st) + (size/2)
				if chrom not in ranges:
					ranges[chrom] = Intersecter()
				else:
					ranges[chrom].add_interval( Interval( mid, mid ) )
					
		print >>sys.stderr, "Done"

		
		print >>sys.stderr, "Calculating RPKM ...",
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
				# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
				exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
				key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
				continue
			mRNA_count=0	#we need to initializ it to 0 for each gene
			mRNA_len=sum(exon_sizes)
			for st,end in zip(exon_starts,exon_ends):
				if chrom in ranges:
					mRNA_count += len(ranges[chrom].find(st,end))					
			mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * cUR_num)
			
			mRNAlen_table[key] = mRNA_len
			RPKM_table[key] = str(mRNA_RPKM)
			rawCount_table[key] = str(mRNA_count)
		print >>sys.stderr, "Done"
		
		print >>RPKM_OUT, '\t'.join(RPKM_head)
		for k in RPKM_table:
			print >>RPKM_OUT, k + '\t',
			print >>RPKM_OUT, str(mRNAlen_table[k]) + '\t',
			print >>RPKM_OUT, str(rawCount_table[k]) + '\t',
			print >>RPKM_OUT, str(RPKM_table[k]) + '\t'
		return RPKM_table
		self.f.seek(0)	

	def strand_specificity(self,refbed,genome,outfile=None,sp="GTAG,GCAG,ATAC"):
		'''Measure the strand specificity of strand specific RNA-seq protocol. For non-splice read,
		use the parental gene as standard, for spliced read, use the splicing motif as strandard'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		if genome is None:
			print >>sys.stderr,"You must specify genome sequence in fasta format\n"
			exit(0)
		
		if outfile is None:
			strand_file = self.fileName + ".strand.infor"
		else:
			strand_file = outfile + ".strand.infor"		
		OUT = open(strand_file,'w')
		print >>OUT,"read_type\tread_id\tread_seq\tchr\tStart\tCigar\tprotocol_strand\tgene_strand"	
		
		transtab = string.maketrans("ACGTNX","TGCANX")
		motif=sp.upper().split(',')
		motif_rev = [m.translate(transtab)[::-1] for m in motif]
		
		#load genome
		print >>sys.stderr, "\tloading "+genome+'...'
		tmp=fasta.Fasta(genome)
		
		#load reference gene model
		gene_ranges={}
		print >>sys.stderr, "reading reference gene model ...",
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
				# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0]
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5]
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
				continue
			if chrom not in gene_ranges:
				gene_ranges[chrom]=Intersecter()
			gene_ranges[chrom].insert(tx_start,tx_end,strand)							
		print >>sys.stderr, "Done"		

		#read SAM 
		
		read_type="unknown"
		strand_from_protocol = 'unknown'
		strand_from_gene='unknown'
		strand_stat=collections.defaultdict(int)
		print >>sys.stderr, "Reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue						#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):continue		#skip multiple mapped reads
			if (flagCode & 0x0100 !=0): continue						#skip non primary hit
			if (flagCode & 0x0200 !=0): continue						#skip QC-failed
			if (flagCode & 0x0400 !=0): continue						#skip PCR artifact
			if (flagCode & 0x0010 !=0): strand_from_protocol = '-'
			if (flagCode & 0x0010 ==0): strand_from_protocol = '+'
			if (flagCode & 0x0040 !=0): read_type="read_1"
			if (flagCode & 0x0080 !=0): read_type="read_2"
			chrom = fields[2]
			comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
			readStart = string.atoi(fields[3])-1
			
			
			#for non spliced read
			if len(comb)==1:
				readEnd = readStart + len(fields[9])
				if chrom in gene_ranges:
					if len(set(gene_ranges[chrom].find(readStart,readEnd)))>1:    
						strand_from_gene="overlap"
					elif len(set(gene_ranges[chrom].find(readStart,readEnd)))==1:
						strand_from_gene = set(gene_ranges[chrom].find(readStart,readEnd)).pop()
					else:
						strand_from_gene="intergenic"
				
				print >>OUT,read_type + '\t' + fields[0] + '\t' + fields[9] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[5] +'\t',
				print >>OUT,strand_from_protocol + '\t' + strand_from_gene		  
		 		strand_stat[read_type + '\t' + strand_from_protocol +'\t' + strand_from_gene] +=1  


			#for spliced read
			if len(comb)>=3:
				splice_strand=[]
				blockStart=[]
				blockSize=[]
				blockEnd =[]
				for i in range(0,len(comb),2):
					blockStart.append(readStart + sum(comb[:i]) )
				for i in range(0,len(comb),2):
					blockSize.append(comb[i])
				blockEnd=map((lambda x,y:x+y),blockStart,blockSize)
				intron_start=blockEnd[:-1]
				intron_end=blockStart[1:]
				for st,end in zip(intron_start,intron_end):
					try:
						splice_motif = str(tmp.fetchSeq(chrom, st, st+2)) + str(tmp.fetchSeq(chrom, end-2,end))
					except:
							print line
					if splice_motif in motif:
						splice_strand.append('+')
					elif splice_motif in motif_rev:
						splice_strand.append('-')
					else:
						splice_strand.append('unknown motif')


				if len(set(splice_strand))>1:
						strand_from_splice = 'unknown motif'
				else:
						strand_from_splice = set(splice_strand).pop()
				print >>OUT,read_type + '\t' + fields[0] + '\t' + fields[9] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[5] +'\t',
				print >>OUT,strand_from_protocol + '\t' + strand_from_splice
				
				strand_stat[read_type + '\t' + strand_from_protocol +'\t' + strand_from_splice] +=1
							
		print >>sys.stderr, "Done"
		
		print "read_type\tstrand_expected\tstrand_observed\tcount"
		for i in sorted(strand_stat):
				print str(i) +'\t' + str(strand_stat[i])
		
	def clipping_profile(self,outfile=None):
		'''calculate profile of soft clipping'''
		if outfile is None:
			out_file1 = self.fileName + ".clipping_profile.xls"
			out_file2 = self.fileName + ".clipping_profile.r"
		else:
			out_file1 = outfile + ".clipping_profile.xls"
			out_file2 = outfile + ".clipping_profile.r"
		
		OUT=open(out_file1,'w')
		ROUT=open(out_file2,'w')
		print >>OUT, "Position\tRead_Total\tRead_clipped"
		soft_p = re.compile(r'(.*?)(\d+)S')
		read_part = re.compile(r'(\d+)[MIS=X]')
		total_read =0
		skip_part_of_read =0
		soft_clip_profile=collections.defaultdict(int)
		
		read_pos=[]
		clip_count=[]
		print >>sys.stderr, "Reading "+ self.fileName + '...'
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue						#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):continue		#skip multiple mapped reads
			if (flagCode & 0x0100 !=0): continue						#skip non primary hit
			if (flagCode & 0x0200 !=0): continue						#skip QC-failed
			if (flagCode & 0x0400 !=0): continue						#skip PCR artifact
			total_read +=1
			m = soft_p.findall(fields[5])
			
			skip_part_of_read =0
			if len(m)==0:					#NO soft clip
				continue
			else:
				for j in m: 
					skip_part_of_read += sum([int(i) for i in read_part.findall(j[0])])
					for n in range(skip_part_of_read,(skip_part_of_read + int(j[1]))):
						soft_clip_profile[n]+=1
					skip_part_of_read += int(j[1])
		for i in soft_clip_profile:
			read_pos.append(str(i))
			clip_count.append(str(soft_clip_profile[i]))
			print >>OUT, str(i) + '\t' + str(total_read) + '\t' + str(soft_clip_profile[i])
		print >>ROUT, "pdf('clipping_profile.pdf')"
		print >>ROUT, "read_pos=c(" + ','.join(read_pos) + ')'
		print >>ROUT, "count=c(" + ','.join(clip_count) + ')'
		print >>ROUT, 'plot(read_pos,1-(count/%d),col="blue",main="clipping profile",xlab="Position of reads",ylab="Mappability",type="b")' % total_read
		print >>ROUT, "dev.off()"
		
	def insertion_profile(self,read_len,outfile=None):
		'''calculate profile of insertion (insertion means insertion to the reference)'''
		if outfile is None:
			out_file1 = self.fileName + ".insertion_profile.xls"
			out_file2 = self.fileName + ".insertion_profile.r"
		else:
			out_file1 = outfile + ".insertion_profile.xls"
			out_file2 = outfile + ".insertion_profile.r"
		
		OUT=open(out_file1,'w')
		ROUT=open(out_file2,'w')
		print >>OUT, "Position\tRead_Total\tRead_clipped"
		soft_p = re.compile(r'(.*?)(\d+)I')
		read_part = re.compile(r'(\d+)[MIS=X]')
		total_read =0
		skip_part_of_read =0
		soft_clip_profile=collections.defaultdict(int)
		print >>sys.stderr, "Reading "+ self.fileName + '...',
		for line in self.f:
			if line.startswith("@"):continue
			fields=line.rstrip('\n ').split()
			flagCode=string.atoi(fields[1])
			if (flagCode & 0x0004) != 0: continue						#skip unmap reads
			if not ParseSAM._uniqueHit_pat.search(line):continue		#skip multiple mapped reads
			if (flagCode & 0x0100 !=0): continue						#skip non primary hit
			if (flagCode & 0x0200 !=0): continue						#skip QC-failed
			if (flagCode & 0x0400 !=0): continue						#skip PCR artifact
			total_read +=1
			m = soft_p.findall(fields[5])
			
			skip_part_of_read =0
			if len(m)==0:					#NO soft clip
				continue
			else:
				for j in m: 
					skip_part_of_read += sum([int(i) for i in read_part.findall(j[0])])
					for n in range(skip_part_of_read,(skip_part_of_read + int(j[1]))):
						soft_clip_profile[n]+=1
					skip_part_of_read += int(j[1])
		for i in range(0,read_len):
			print >>OUT, str(i) + '\t' + str(total_read) + '\t' + str(soft_clip_profile[i])

class ParseBAM:
	'''This class provides fuctions to parsing/processing/transforming SAM or BAM files. The input
	file could be either SAM or BAM format file'''
	
	multi_hit_tags=['H0','H1','H2','IH','NH']
	def __init__(self,inputFile):
		'''constructor. input could be bam or sam'''
		try:
			self.samfile = pysam.Samfile(inputFile,'rb')
			if len(self.samfile.header) ==0:
				print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
				sys.exit(1)
			self.bam_format = True
		except:
			self.samfile = pysam.Samfile(inputFile,'r')
			if len(self.samfile.header) ==0:
				print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
				sys.exit(1)
			self.bam_format = False

	def stat (self):
		'''Calculate mapping statistics'''
		R_total=0
		R_qc_fail=0
		R_duplicate=0
		R_nonprimary=0
		R_unmap =0
		
		R_multipleHit=0
		R_uniqHit=0	#all the following count should be sum to uniqHit
		
		R_read1=0
		R_read2=0
		R_reverse =0
		R_forward=0
		R_nonSplice=0
		R_splice=0
		R_properPair =0 

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				R_total +=1
				if aligned_read.is_qcfail:			#skip QC fail read
					R_qc_fail +=1
					continue
				if aligned_read.is_duplicate:		#skip duplicate read
					R_duplicate +=1
					continue
				if aligned_read.is_secondary:		#skip non primary hit
					R_nonprimary +=1
					continue
				if aligned_read.is_unmapped:		#skip unmap read
					R_unmap +=1
					continue		
				
				
				if len(aligned_read.tags) > 0:
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:
					R_multipleHit +=1
					continue						#skip multiple map read				
				if flag==0:
					R_uniqHit +=1
					if aligned_read.is_read1:
						R_read1 +=1
					if aligned_read.is_read2:
						R_read2 +=1
					if aligned_read.is_reverse:
						R_reverse +=1
					else:
						R_forward +=1
					cigar_str = cigar.list2str(aligned_read.cigar)
					introns = cigar.fetch_intron('chr1', aligned_read.pos, cigar_str)
					if len(introns)==0:
						R_nonSplice +=1
					if len(introns)>=1:
						R_splice +=1
					if aligned_read.is_proper_pair:
						R_properPair +=1
				
		except StopIteration:
			print >>sys.stderr, "Done"		
		#self.samfile.seek(current_pos)
				
		print >>sys.stderr,"\n#=================================================="
		print >>sys.stderr, "%-30s%d" % ("Total Reads (Records):",R_total)
		print >>sys.stderr, "\n",
		print >>sys.stderr, "%-30s%d" % ("QC failed:",R_qc_fail)
		print >>sys.stderr, "%-30s%d" % ("Optical/PCR duplicate:", R_duplicate)
		print >>sys.stderr, "%-30s%d" % ("Non Primary Hits", R_nonprimary)
		print >>sys.stderr, "%-30s%d" % ("Unmapped reads:",R_unmap)
		print >>sys.stderr, "%-30s%d" % ("Multiple mapped reads:",R_multipleHit)
		print >>sys.stderr, "\n",
		print >>sys.stderr, "%-30s%d" % ("Uniquely mapped:",R_uniqHit)
		print >>sys.stderr, "%-30s%d" % ("Read-1:",R_read1)
		print >>sys.stderr, "%-30s%d" % ("Read-2:",R_read2)
		print >>sys.stderr, "%-30s%d" % ("Reads map to '+':",R_forward)
		print >>sys.stderr, "%-30s%d" % ("Reads map to '-':",R_reverse)
		print >>sys.stderr, "%-30s%d" % ("Non-splice reads:",R_nonSplice)
		print >>sys.stderr, "%-30s%d" % ("Splice reads:",R_splice)	
		print >>sys.stderr, "%-30s%d" % ("Reads mapped in proper pairs:",R_properPair)
	
	def configure_experiment(self,refbed,sample_size = 200000):
		'''Given a BAM/SAM file, this function will try to guess the RNA-seq experiment:
			1) single-end or pair-end
			2) strand_specific or not
			3) if it is strand-specific, what's the strand_ness of the protocol
		'''
		
			#how many reads you want to sample
		count =0
		p_strandness=collections.defaultdict(int)
		s_strandness=collections.defaultdict(int)
		#load reference gene model
		gene_ranges={}
		print >>sys.stderr, "Reading reference gene model " + refbed + ' ...',
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
				# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0]
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5]
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
				continue
			if chrom not in gene_ranges:
				gene_ranges[chrom]=Intersecter()
			gene_ranges[chrom].insert(tx_start,tx_end,strand)							
		print >>sys.stderr, "Done"		
		
		#read SAM/BAM file
		#current_pos = self.samfile.tell()
		print >>sys.stderr, "Loading SAM/BAM file ... ",
		try:
			while(1):
				if count >= sample_size:
					break
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:			#skip low quanlity
					continue
				if aligned_read.is_duplicate:		#skip duplicate read
					continue
				if aligned_read.is_secondary:		#skip non primary hit
					continue
				if aligned_read.is_unmapped:		#skip unmap read
					continue		
				if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:continue						#skip multiple map read		
														
				
				chrom = self.samfile.getrname(aligned_read.tid)
				if aligned_read.is_paired:
					if aligned_read.is_read1:
						read_id = '1'
					if aligned_read.is_read2:
						read_id = '2'
					if aligned_read.is_reverse:
						map_strand = '-'
					else:
						map_strand = '+'
					readStart = aligned_read.pos
					readEnd = readStart + aligned_read.qlen
					if chrom in gene_ranges:
						if len(set(gene_ranges[chrom].find(readStart,readEnd)))==1:
							strand_from_gene = set(gene_ranges[chrom].find(readStart,readEnd)).pop()
							p_strandness[read_id + map_strand + strand_from_gene]+=1
							count += 1
				else:
					if aligned_read.is_reverse:
						map_strand = '-'
					else:
						map_strand = '+'					
					readStart = aligned_read.pos
					readEnd = readStart + aligned_read.qlen
					if chrom in gene_ranges:
						if len(set(gene_ranges[chrom].find(readStart,readEnd)))==1:
							strand_from_gene = set(gene_ranges[chrom].find(readStart,readEnd)).pop()
							s_strandness[map_strand + strand_from_gene]+=1
							count += 1

		except StopIteration:
			print >>sys.stderr, "Finished"		
		#self.samfile.seek(current_pos)
		
		print >>sys.stderr, "Total " + str(count) + " usable reads were sampled"
		protocol="unknown"
		strandness=None
		spec1=0.0
		spec2=0.0
		other=0.0
		if len(p_strandness) >0 and len(s_strandness) ==0 :
			protocol="PairEnd"
			#for k,v in p_strandness.items():
			#	print >>sys.stderr, k + '\t' + str(v)
			spec1= (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
			spec2= (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
			other = 1-spec1-spec2
			
		elif len(s_strandness) >0 and len(p_strandness) ==0 :
			protocol="SingleEnd"
			#for k,v in s_strandness.items():
			#	print  >>sys.stderr, k + '\t' + str(v)
			spec1 = (s_strandness['++'] + s_strandness['--'])/float(sum(s_strandness.values()))
			spec2 = (s_strandness['+-'] + s_strandness['-+'])/float(sum(s_strandness.values()))
			other = 1-spec1-spec2
		else:
			protocol="Mixture"
			spec1 = "NA"
			spec2 = "NA"
			other = "NA"
		return [protocol,spec1,spec2,other]

	def bamTowig(self,outfile,chrom_sizes, skip_multi=True,strand_rule=None):
		"""Convert BAM/SAM file to wig file. chrom_size is dict with chrom as key and chrom_size as value
		strandRule should be determined from \"infer_experiment\". such as \"1++,1--,2+-,2-+\""""
		
		#strand_rule={'1+':'-','1-':'+','2+':'+','2-':'-'}
		strandRule={}
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of option :'strand_rule' " + strand_rule
			sys.exit(1)
		if len(strandRule) == 0:
			FWO = open(outfile + '.wig','w')
		else:
			FWO = open(outfile + '.Forward.wig','w')
			RVO = open(outfile + '.Reverse.wig','w')
		
		
		read_id=''
		for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
			try:
				self.samfile.fetch(chr_name,0,chr_size)
			except:
				print >>sys.stderr, "No alignments for " + chr_name + '. skipped'
				continue
			print >>sys.stderr, "Processing " + chr_name + " ..."
			if len(strandRule) == 0: FWO.write('variableStep chrom='+chr_name+'\n')
			else:
				FWO.write('variableStep chrom='+chr_name+'\n')
				RVO.write('variableStep chrom='+chr_name+'\n')
			Fwig = collections.defaultdict(int)
			Rwig = collections.defaultdict(int)
			alignedReads = self.samfile.fetch(chr_name,0,chr_size)		
			for aligned_read in alignedReads:
				flag=0
				if aligned_read.is_qcfail:continue			#skip low quanlity
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				
				if skip_multi:
					if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
						for i in aligned_read.tags:
							if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
								flag=1						#multiple hit read
								break
					if flag==1:continue						#skip multiple map read		
				
				if aligned_read.is_paired:
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'
				
				key = read_id + map_strand
				
				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)
				for block in cigar.fetch_exon(chr_name, hit_st, cigar_str): 
					for pos in range(block[1]+1,block[2]+1):	
						if len(strandRule) == 0: Fwig[pos] +=1.0	#this is NOT strand specific. everythign into Fwig
						else:										#this is strand specific. separate Fwig and Rwig
							if strandRule[key] == '+':Fwig[pos] +=1.0
							if strandRule[key] == '-':Rwig[pos] -=1.0
					
			if len(strandRule) == 0:							#this is NOT strand specific.
				for pos in sorted (Fwig.keys()):
					print >>FWO, "%d\t%.2f" % (pos,Fwig[pos])
			else:
				for pos in sorted (Fwig.keys()):
					print >>FWO, "%d\t%.2f" % (pos,Fwig[pos])
				for pos in sorted (Rwig.keys()):
					print >>RVO, "%d\t%.2f" % (pos,Rwig[pos])
				
					
	def calculate_rpkm(self,geneFile,outfile,strand_rule=None):
		'''calculate RPKM vaues. For single end RNA-seq, if it is strand specific, we assume that
		read plus mapped indicates a gene on plus strand.(similar to minus). 
		Advantages: works for both SAM and BAM
					works for both sorted and unsorted BAM/SAM file
					works for both index or unindexed BAM/SAM file
					much faster than indexing bam file
		Disadvantage: random access BAM file was disabled, thus large mount of RAM is required
		
		strand_rule: could be the following values:
			'1++,1--,2+-,2-+
			'1+-,1-+,2++,2--
			'++,--'
			'+-,-+'
			None
		'''
		
		strandRule={}
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of option :'strand_rule' " + strand_rule
			sys.exit(1)
		
		uniq_read=0
		total_tags=0
		plus_ranges={}
		minus_ranges={}
		unstrand_ranges={}
		
		rpkm_value={}
		
		RPKM_OUT = open(outfile,'w')
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		
		#current_pos = self.samfile.tell()
		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read

				if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:continue						#skip multiple map read		
				
				uniq_read +=1
				
				if aligned_read.is_paired:
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				else:
					read_id = ''
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'
				
				strand_key = read_id + map_strand
				
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)		
				exon_blocks = cigar.fetch_exon(chrom, hit_st, cigar_str)
				total_tags += len(exon_blocks)
						
				#construct bitset
				if strand_rule is not None:	
					if strandRule[strand_key] == '+':
						for block in exon_blocks:
							mid = block[1] + int((block[2] - block[1])/2)
							if chrom not in plus_ranges:plus_ranges[chrom] = Intersecter()
							else:plus_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
					elif strandRule[strand_key] == '-':
						for block in exon_blocks:
							mid = block[1] + int((block[2] - block[1])/2)
							if chrom not in minus_ranges:minus_ranges[chrom] = Intersecter()	
							else:minus_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
				elif strand_rule is None:	
					for block in exon_blocks:
						mid = block[1] + int((block[2] - block[1])/2)
						if chrom not in unstrand_ranges:unstrand_ranges[chrom] = Intersecter()
						else:unstrand_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
					
		except StopIteration:
			print >>sys.stderr, "Done"
		#self.samfile.seek(current_pos)
		print >>RPKM_OUT, "#Total uniquely mapped reads = " + str(uniq_read)
		print >>RPKM_OUT, "#Total fragments = " + str(total_tags)
		print >>sys.stderr, "Assign reads to "+ geneFile + '...',
		for line in open(geneFile,'r'):
			try:
				if line.startswith('#'):continue
				if line.startswith('track'):continue
				if line.startswith('browser'):continue   
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
				exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
				intron_starts = exon_ends[:-1]
				intron_ends=exon_starts[1:]
				key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue

					
			mRNA_count=0
			mRNA_len=sum(exon_sizes)
				
			if (strand_rule is not None) and (strand == '-'):
				intronNum=len(intron_starts)
				exonNum=len(exon_starts)
				
				# assign reads to intron	
				for st,end in zip(intron_starts,intron_ends):
					if chrom in minus_ranges:
						hits= len(minus_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						intronNum -= 1
				# assign reads to exon				
				for st,end in zip(exon_starts,exon_ends):
					if chrom in minus_ranges:
						hits= len(minus_ranges[chrom].find(st,end))					
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						exonNum -= 1
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
			elif (strand_rule is not None) and (strand == '+'):
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in plus_ranges:
						hits= len(plus_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						intronNum += 1    
				for st,end in zip(exon_starts,exon_ends):
					if chrom in plus_ranges:
						hits= len(plus_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						exonNum += 1		
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
			elif strand_rule is None:
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in unstrand_ranges:
						hits= len(unstrand_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						intronNum += 1    
				for st,end in zip(exon_starts,exon_ends):
					if chrom in unstrand_ranges:
						hits= len(unstrand_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						exonNum += 1		
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
		print >>sys.stderr, "Done"

	def readsNVC(self,outfile=None,nx=True):
		'''for each read, calculate nucleotide frequency vs position'''
		if outfile is None:
			outfile1 = self.fileName + ".NVC.xls"
			outfile2 = self.fileName +".NVC_plot.r"
		else:
			outfile1 = outfile + ".NVC.xls"
			outfile2 = outfile +".NVC_plot.r"
		FO=open(outfile1,'w')
		RS=open(outfile2,'w')
		PPcount=0
		
		transtab = string.maketrans("ACGTNX","TGCANX")
		base_freq=collections.defaultdict(int)
		a_count=[]
		c_count=[]
		g_count=[]
		t_count=[]
		n_count=[]
		x_count=[]
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				aligned_read = self.samfile.next()
				#if aligned_read.is_unmapped:continue	#skip unmapped read
				#if aligned_read.is_qcfail:continue	#skip low quality
				RNA_read = aligned_read.seq.upper()		
				if aligned_read.is_reverse:
					RNA_read = RNA_read.translate(transtab)[::-1]
				for i,j in enumerate(RNA_read):
					key = str(i) + j
					base_freq[key] += 1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print >>sys.stderr, "generating data matrix ..."
		print >>FO, "Position\tA\tC\tG\tT\tN\tX"
		for i in xrange(len(RNA_read)):
			print  >>FO, str(i) + '\t',
			print  >>FO, str(base_freq[str(i) + "A"]) + '\t',
			a_count.append(str(base_freq[str(i) + "A"]))
			print  >>FO, str(base_freq[str(i) + "C"]) + '\t',
			c_count.append(str(base_freq[str(i) + "C"]))
			print  >>FO, str(base_freq[str(i) + "G"]) + '\t',
			g_count.append(str(base_freq[str(i) + "G"]))
			print  >>FO, str(base_freq[str(i) + "T"]) + '\t',
			t_count.append(str(base_freq[str(i) + "T"]))
			print  >>FO, str(base_freq[str(i) + "N"]) + '\t',
			n_count.append(str(base_freq[str(i) + "N"]))
			print  >>FO, str(base_freq[str(i) + "X"]) + '\t'
			x_count.append(str(base_freq[str(i) + "X"]))
		FO.close()
		
		#generating R scripts
		print >>sys.stderr, "generating R script  ..."
		print >>RS, "position=c(" + ','.join([str(i) for i in xrange(len(RNA_read))]) + ')'
		print >>RS, "A_count=c(" + ','.join(a_count) + ')'
		print >>RS, "C_count=c(" + ','.join(c_count) + ')'
		print >>RS, "G_count=c(" + ','.join(g_count) + ')'
		print >>RS, "T_count=c(" + ','.join(t_count) + ')'
		print >>RS, "N_count=c(" + ','.join(n_count) + ')'
		print >>RS, "X_count=c(" + ','.join(x_count) + ')'
		
		if nx:
			print >>RS, "total= A_count + C_count + G_count + T_count + N_count + X_count"
			print >>RS, "ym=max(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total) + 0.05"
			print >>RS, "yn=min(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total)"
			
			print >>RS, 'pdf(\"%s\")' % (outfile +".NVC_plot.pdf")
			print >>RS, 'plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")'
			print >>RS, 'lines(position,T_count/total,type="o",pch=20,col="red")'
			print >>RS, 'lines(position,G_count/total,type="o",pch=20,col="blue")'
			print >>RS, 'lines(position,C_count/total,type="o",pch=20,col="cyan")'
			print >>RS, 'lines(position,N_count/total,type="o",pch=20,col="black")'		
			print >>RS, 'lines(position,X_count/total,type="o",pch=20,col="grey")'	
			print >>RS, 'legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C","N","X"),col=c("dark green","red","blue","cyan","black","grey"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan","black","grey"))'
			print >>RS, "dev.off()"
		else:
			print >>RS, "total= A_count + C_count + G_count + T_count"
			print >>RS, "ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05"
			print >>RS, "yn=min(A_count/total,C_count/total,G_count/total,T_count/total)"
		
			print >>RS, 'pdf(\"%s\")' % (outfile +".NVC_plot.pdf")
			print >>RS, 'plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")'
			print >>RS, 'lines(position,T_count/total,type="o",pch=20,col="red")'
			print >>RS, 'lines(position,G_count/total,type="o",pch=20,col="blue")'
			print >>RS, 'lines(position,C_count/total,type="o",pch=20,col="cyan")'
			print >>RS, 'legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))'
			print >>RS, "dev.off()"
		
		RS.close()
		#self.f.seek(0)

	def readsQual_boxplot(self,outfile,shrink=1000):
		'''calculate phred quality score for each base in read (5->3)'''

		output = outfile + ".qual.r"
		FO=open(output,'w')

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		quality = collections.defaultdict(dict)	#read_pos=>quality score=>count
		q_max = -1
		q_min = 10000
		q_list=[]
		i_box={}	#key is read postion,value is 
		try:
			while(1):
				aligned_read = self.samfile.next()
				
				#if aligned_read.is_unmapped:continue	#skip unmapped read
				#if aligned_read.is_qcfail:continue		#skip low quality
				
				qual_str = aligned_read.qqual
				read_len = aligned_read.rlen
				if aligned_read.is_reverse:
					qual_str = qual_str[::-1]

				for i,j in enumerate(qual_str):
					q=ord(j)-33
					if q > q_max: q_max = q
					if q < q_min: q_min = q
					try:
						quality[i][q] += 1
					except:
						quality[i][q] = 1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		for p in range(0,read_len):
			#print str(p) + ':',
			val=[]
			occurrence=[]
			for q in range(q_min,q_max+1):
				if quality.has_key(p) and quality[p].has_key(q):
					val.append(str(q))				
					occurrence.append(str(quality[p][q]))	
					q_list.append(str(quality[p][q]))
				else:
					q_list.append(str(0))
			i_box[p] = 'rep(c(' + ','.join(val) + '),times=c(' + ','.join(occurrence) + ')/' + str(shrink)+ ')'
		
		
		#generate R script for boxplot
		print >>FO, "pdf(\'%s\')" % (outfile + ".qual.boxplot.pdf")
		for i in sorted(i_box):
			print >>FO,'p'+str(i) + '<-' + i_box[i]
		print >>FO, 'boxplot(' + ','.join(['p'+str(i) for i in i_box]) + ',xlab=\"Position of Read(5\'->3\')\",ylab=\"Phred Quality Score\",outline=F' + ')'
		print >>FO,"dev.off()"
		
		
		#generate R script for heatmap
		print >>FO, '\n'
		print >>FO, "pdf(\'%s\')" % (outfile + ".qual.heatmap.pdf")
		print >>FO, "qual=c(" + ','.join(q_list)  + ')'
		print >>FO, "mat=matrix(qual,ncol=%s,byrow=F)" % (read_len)
		print >>FO, 'Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c(\'spline\'))'
		print >>FO, "heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from=%s,to=%s),col = Lab.palette(256),scale=\"none\" )" % (q_min,q_max)
		print >>FO, 'dev.off()'
		
		
	def readGC(self,outfile=None):
		'''GC content distribution of reads'''
		if outfile is None:
			outfile1 = self.fileName + ".GC.xls"
			outfile2 = self.fileName +".GC_plot.r"
		else:
			outfile1 = outfile + ".GC.xls"
			outfile2 = outfile + ".GC_plot.r"
		FO=open(outfile1,'w')
		RS=open(outfile2,'w')
		
		gc_hist=collections.defaultdict(int)	#key is GC percent, value is count of reads

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.is_unmapped:continue	#skip unmapped read
				if aligned_read.is_qcfail:continue		#skip low quality
				RNA_read = aligned_read.seq.upper()		
				gc_percent = "%4.2f" % ((RNA_read.count('C') + RNA_read.count('G'))/(len(RNA_read)+0.0)*100)
				#print gc_percent
				gc_hist[gc_percent] += 1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print >>sys.stderr, "writing GC content ..."	
		print >>FO, "GC%\tread_count"
		for i in gc_hist.keys():
			print >>FO, i + '\t' + str(gc_hist[i])
			
		print >>sys.stderr, "writing R script ..."
		print >>RS, "pdf(\"%s\")" % (outfile +  ".GC_plot.pdf")
		print >>RS, 'gc=rep(c(' + ','.join([i for i in gc_hist.keys()]) + '),' + 'times=c(' + ','.join([str(i) for i in gc_hist.values()]) + '))'
		print >>RS, 'hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100
		#print >>RS, "lines(density(gc),col='red')"
		print >>RS ,"dev.off()"		
		#self.f.seek(0)


	def readDupRate(self,outfile=None,up_bound=500):
		'''Calculate reads's duplicate rates'''
		if outfile is None:
			outfile1 = self.fileName + ".seq.DupRate.xls"
			outfile2 = self.fileName + ".pos.DupRate.xls"
			outfile3 = self.fileName + ".DupRate_plot.r"
		else:
			outfile1 = outfile + ".seq.DupRate.xls"
			outfile2 = outfile + ".pos.DupRate.xls"
			outfile3 = outfile +".DupRate_plot.r"
		SEQ=open(outfile1,'w')
		POS=open(outfile2,'w')
		RS=open(outfile3,'w')
		
		seqDup=collections.defaultdict(int)
		posDup=collections.defaultdict(int)
		
		seqDup_count=collections.defaultdict(int)
		posDup_count=collections.defaultdict(int)

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				exon_boundary=""
				aligned_read = self.samfile.next()
				if aligned_read.is_unmapped:continue	#skip unmapped read
				if aligned_read.is_qcfail:continue		#skip low quality
				RNA_read = aligned_read.seq.upper()		
				seqDup[RNA_read] +=1					#key is read sequence

				chrom = self.samfile.getrname(aligned_read.tid)
				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)		
				exon_blocks = cigar.fetch_exon(chrom, hit_st, cigar_str)				
				for ex in exon_blocks:
					exon_boundary += str(ex[1]) + '-' + str(ex[2]) + ":"
				key = chrom + ":" + str(hit_st) + ":" + exon_boundary
				posDup[key] +=1

		except StopIteration:
			print >>sys.stderr, "Done"

		print >>sys.stderr, "report duplicte rate based on sequence ..."
		print >>SEQ, "Occurrence\tUniqReadNumber"
		for i in seqDup.values():			#key is occurence, value is uniq reads number (based on seq)
			seqDup_count[i] +=1
		for k in sorted(seqDup_count.iterkeys()):	
			print >>SEQ, str(k) +'\t'+ str(seqDup_count[k])
		SEQ.close()
		
		print >>sys.stderr, "report duplicte rate based on mapping  ..."
		print >>POS, "Occurrence\tUniqReadNumber"
		for i in posDup.values():			#key is occurence, value is uniq reads number (based on coord)
			posDup_count[i] +=1
		for k in sorted(posDup_count.iterkeys()):	
			print >>POS, str(k) +'\t'+ str(posDup_count[k])
		POS.close()
		
		
		print >>sys.stderr, "generate R script ..."
		print >>RS, "pdf(\'%s\')" % (outfile +".DupRate_plot.pdf")
		print >>RS, "par(mar=c(5,4,4,5),las=0)"
		print >>RS, "seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
		print >>RS, "seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
		print >>RS, "pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.iterkeys()) ]) + ')'
		print >>RS, "pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.iterkeys()) ]) + ')'
		print >>RS, "plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Frequency',pch=4,cex=0.8,col='blue',xlim=c(1,%d),yaxt='n')" % up_bound
		print >>RS, "points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')"
		print >>RS, 'ym=floor(max(log10(pos_uniqRead)))'
		print >>RS, "legend(%d,ym,legend=c('Sequence-base','Mapping-base'),col=c('red','blue'),pch=c(4,20))" % max(up_bound-200,1)
		print >>RS, 'axis(side=2,at=0:ym,labels=0:ym)'
		print >>RS, 'axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead)),round(pos_uniqRead[2]*100/sum(pos_uniqRead)),round(pos_uniqRead[3]*100/sum(pos_uniqRead)),round(pos_uniqRead[4]*100/sum(pos_uniqRead))))'
		print >>RS, 'mtext(4, text = "Reads %", line = 2)'
		print >>RS, 'dev.off()'
		#self.f.seek(0)

	def clipping_profile(self,outfile):
		'''calculate profile of soft clipping'''
		out_file1 = outfile + ".clipping_profile.xls"
		out_file2 = outfile + ".clipping_profile.r"
		
		OUT=open(out_file1,'w')
		ROUT=open(out_file2,'w')
		print >>OUT, "Position\tRead_Total\tRead_clipped"
		soft_p = re.compile(r'(.*?)(\d+)S')
		read_part = re.compile(r'(\d+)[MIS=X]')
		total_read =0
		skip_part_of_read =0
		soft_clip_profile=collections.defaultdict(int)
		
		read_pos=[]
		clip_count=[]

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				exon_boundary=""
				aligned_read = self.samfile.next()
				if aligned_read.is_unmapped:continue	#skip unmapped read
				if aligned_read.is_qcfail:continue		#skip low quality
				
				total_read +=1
				cigar_str = cigar.list2str(aligned_read.cigar)	
				m = soft_p.findall(cigar_str)
				skip_part_of_read =0
				if len(m)==0:continue					#NO soft clip
				else:
					for j in m: 
						skip_part_of_read += sum([int(i) for i in read_part.findall(j[0])])
						for n in range(skip_part_of_read,(skip_part_of_read + int(j[1]))):
							soft_clip_profile[n]+=1
						skip_part_of_read += int(j[1])
		except StopIteration:
			print >>sys.stderr, "Done"
			

		for i in soft_clip_profile:
			read_pos.append(str(i))
			clip_count.append(str(soft_clip_profile[i]))
			print >>OUT, str(i) + '\t' + str(total_read) + '\t' + str(soft_clip_profile[i])
		print >>ROUT, "pdf('clipping_profile.pdf')"
		print >>ROUT, "read_pos=c(" + ','.join(read_pos) + ')'
		print >>ROUT, "count=c(" + ','.join(clip_count) + ')'
		print >>ROUT, 'plot(read_pos,1-(count/%d),col="blue",main="clipping profile",xlab="Position of reads",ylab="Mappability",type="b")' % total_read
		print >>ROUT, "dev.off()"

	def coverageGeneBody(self,refbed,outfile):
		'''Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divided
		into 100 regsions'''
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		OUT1 = open(outfile + ".geneBodyCoverage_plot.r",'w')
		OUT2 = open(outfile + ".geneBodyCoverage.txt",'w')

		ranges={}
		totalReads=0
		fragment_num=0		#splice reads will counted twice
		rpkm={}
		
		#read SAM 
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				totalReads +=1

				chrom = self.samfile.getrname(aligned_read.tid).upper()
				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)		
				exon_blocks = cigar.fetch_exon(chrom, hit_st, cigar_str)			
				fragment_num += len(exon_blocks)
			
				for exon in exon_blocks:
					if chrom not in ranges:
						ranges[chrom] = Intersecter()
					else:
						ranges[chrom].add_interval( Interval( exon[1], exon[2] ) )
		except StopIteration:
			print >>sys.stderr, "Done"

		print >>sys.stderr, "calculating coverage over gene body ..."
		coverage=collections.defaultdict(int)
		flag=0
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
            	# Parse fields from gene tabls
				fields = line.split()
				chrom     = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end    = int( fields[2] )
				geneName      = fields[3]
				strand    = fields[5]
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue
			gene_all_base=[]
			percentile_base=[]
			mRNA_len =0
			flag=0
			for st,end in zip(exon_starts,exon_ends):
				gene_all_base.extend(range(st+1,end+1))		#0-based coordinates on genome
				mRNA_len = len(gene_all_base)
				if mRNA_len <100:
					flag=1
					break
			if flag==1: continue
			if strand == '-':
				gene_all_base.sort(reverse=True)			#deal with gene on minus stand
			else:
				gene_all_base.sort(reverse=False)
			percentile_base = mystat.percentile_list (gene_all_base)	#get 101 points from each gene's coordinates
			
			for i in range(0,len(percentile_base)):
				if chrom in ranges:
					coverage[i] += len(ranges[chrom].find(percentile_base[i], percentile_base[i]+1))
		x_coord=[]
		y_coord=[]
		print >>OUT2, "Total reads: " + str(totalReads)
		print >>OUT2, "Fragment number: " + str(fragment_num)
		print >>OUT2, "percentile\tcount"
		for i in coverage:
			x_coord.append(str(i))
			y_coord.append(str(coverage[i]))
			print >>OUT2, str(i) + '\t' + str(coverage[i])
		print >>OUT1, "pdf(\'%s\')" % (outfile + ".geneBodyCoverage.pdf")
		print >>OUT1, "x=0:100"
		print >>OUT1, "y=c(" + ','.join(y_coord) + ')'
		print >>OUT1, "plot(x,y,xlab=\"percentile of gene body (5'->3')\",ylab='read number',type='s')"
		print >>OUT1, "dev.off()"
					
	def mRNA_inner_distance(self,outfile,refbed,low_bound=0,up_bound=1000,step=10):
		'''estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length'''
		
		out_file1 = outfile + ".inner_distance.txt"	
		out_file2 = outfile + ".inner_distance_freq.txt"
		out_file3 = outfile + ".inner_distance_plot.r"
		
		FO=open(out_file1,'w')
		FQ=open(out_file2,'w')
		RS=open(out_file3,'w')
		
		fchrom="chr100"		#this is the fake chromosome
		ranges={}
		ranges[fchrom]=Intersecter()
		
		window_left_bound = range(low_bound,up_bound,step)
 		frag_size=0
		
		inner_distance_bitsets=BinnedBitSet()
		tmp = BinnedBitSet()
		tmp.set_range(0,0)
 		pair_num=0.0
 		sizes=[]
 		counts=[]
 		count=0
 		
 		print >>sys.stderr, "Get intron regions from " + refbed + " ..."
 		bed_obj = BED.ParseBED(refbed)
 		ref_exons = []
 		
 		for exn in bed_obj.getExon():
 			ref_exons.append([exn[0].upper(), exn[1], exn[2]])
		exon_bitsets = binned_bitsets_from_list(ref_exons)		
		
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				splice_intron_size=0
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				if aligned_read.is_read2:continue	
				if not aligned_read.is_paired: continue		#skip single map read

				if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:continue						#skip multiple map read		

				pair_num +=1
				read1_len = aligned_read.qlen
				read1_start = aligned_read.pos
				read2_start = aligned_read.mpos
				
				if read1_start > read2_start:continue		#only consider read1_start < read2_start
				
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				cigar_str = cigar.list2str(aligned_read.cigar)		
				intron_blocks = cigar.fetch_intron(chrom, read1_start, cigar_str)				
				for intron in intron_blocks:
					splice_intron_size += intron[2] - intron[1]
				
				read1_end = read1_start + read1_len + splice_intron_size		
				#if read1_end > read2_start: continue
				
				inner_distance = read2_start - read1_end +1
				if inner_distance > 0: 
					if chrom in exon_bitsets:
						size =0 
						inner_distance_bitsets.set_range(read1_end, read2_start-read1_end)
						inner_distance_bitsets.iand(exon_bitsets[chrom])
						end=0
						while 1:
							start = inner_distance_bitsets.next_set( end )
							if start == inner_distance_bitsets.size: break
							end = inner_distance_bitsets.next_clear( start )
							size += (end - start) +1
						inner_distance_bitsets.iand(tmp)											#clear BinnedBitSet
						if size == inner_distance:
							FO.write(aligned_read.qname + '\t' + str(size) + '\tPE_within_same_exon\n')
							ranges[fchrom].add_interval( Interval( size-1, size ) )
						elif size < inner_distance and size >1:
							FO.write(aligned_read.qname + '\t' + str(size) + '\tPE_within_diff_exon\n')
							ranges[fchrom].add_interval( Interval( size-1, size ) )						
						elif size == 1:
							FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tPE_not_in_gene_region\n')
							ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
						
					else:
						FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tPE_chrom_has_no_refgene\n')
						ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
				else:
					FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tPE_reads_overlap\n')
					ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
					
		except StopIteration:
			print >>sys.stderr, "Done"
		
		
		if pair_num==0:
			print >>sys.stderr, "Cannot find paired reads"
			sys.exit(0)
		#print >>FQ, "Total paired read " + str(pair_num)
		for st in window_left_bound:
			sizes.append(str(st + step/2))
			count = str(len(ranges[fchrom].find(st,st + step)))
			counts.append(count)
			print >>FQ, str(st) + '\t' + str(st+step) +'\t' + count		
		
		print >>RS, "pdf(\'%s\')" % (outfile + ".inner_distance_plot.pdf")
		#print >>RS, "par(mfrow=c(2,1),cex.main=0.8,cex.lab=0.8,cex.axis=0.8,mar=c(4,4,4,1))"
		#print >>RS, 'pie(c(%d,%d,%d),col=rainbow(3),cex=0.5,radius=1,main="Total %d fragments",labels=c("fraSize <= %d\\n(%4.2f%%)","fragSize > %d\\n(%4.2f%%)","%d < fragSize <= %d\\n(%4.2f%%)"), density=rep(80,80,80),angle=c(90,140,170))' % (ultra_low, ultra_high, pair_num -ultra_low -ultra_high, pair_num, low_bound, ultra_low*100/pair_num, up_bound, ultra_high*100/pair_num, low_bound, up_bound, 100-ultra_low*100/pair_num - ultra_high*100/pair_num)
		print >>RS, 'fragsize=rep(c(' + ','.join(sizes) + '),' + 'times=c(' + ','.join(counts) + '))'
		print >>RS, 'frag_sd = round(sd(fragsize))'
		print >>RS, 'frag_mean = round(mean(fragsize))'
		print >>RS, 'hist(fragsize,probability=T,breaks=%d,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(window_left_bound)
		print >>RS, "lines(density(fragsize,bw=%d),col='red')" % (2*step)
		print >>RS ,"dev.off()"
		FO.close()
		FQ.close()
		RS.close()
		#self.f.seek(0)

	def annotate_junction(self,refgene,outfile,min_intron=50):
		'''Annotate splicing junctions in BAM or SAM file. Note that a (long) read might have multiple splicing
		events  (splice multiple times), and the same splicing events can be consolidated into a single
		junction'''
		
		out_file = outfile + ".junction.xls"
		out_file2 = outfile + ".junction_plot.r"
		if refgene is None:
			print >>sys.stderr, "You must provide reference gene model in bed format."
			sys.exit(1)
		OUT = open(out_file,'w')
		ROUT = open(out_file2,'w')
		
		#reading reference gene model
		refIntronStarts=collections.defaultdict(dict)
		refIntronEnds=collections.defaultdict(dict)	
		total_junc =0
		novel35_junc =0
		novel3or5_junc =0
		known_junc =0
		splicing_events=collections.defaultdict(int)	
		
		print >>sys.stderr, "Reading reference bed file: ",refgene, " ... ",
		for line in open(refgene,'r'):
			if line.startswith(('#','track','browser')):continue  
           	# Parse fields from gene tabls
			fields = line.split()
			if(len(fields)<12):
				print >>sys.stderr, "Invalid bed line (skipped):",line,
				continue
			chrom     = fields[0].upper()
			tx_start = int( fields[1] )
			tx_end   = int( fields[2] )
			if int(fields[9] ==1):
				continue    	
			
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			for i_st,i_end in zip (intron_start, intron_end):
				refIntronStarts[chrom][i_st] =i_st
				refIntronEnds[chrom][i_end] =i_end			
		print >>sys.stderr,"Done"
		
		#reading input SAM file
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read

				if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:continue						#skip multiple map read		

				chrom = self.samfile.getrname(aligned_read.tid).upper()
				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)		
				intron_blocks = cigar.fetch_intron(chrom, hit_st, cigar_str)			
				if len(intron_blocks)==0:
					continue
				for intrn in intron_blocks:
					total_junc +=1
					if intrn[2] - intrn[1] < min_intron:continue
					splicing_events[intrn[0] + ":" + str(intrn[1]) + ":" + str(intrn[2])] += 1
					if (refIntronStarts[chrom].has_key(intrn[1]) and refIntronEnds[chrom].has_key(intrn[2])):
						known_junc +=1																		#known both
					elif (not refIntronStarts[chrom].has_key(intrn[1]) and not refIntronEnds[chrom].has_key(intrn[2])):
						novel35_junc +=1																
					else:
						novel3or5_junc +=1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print "total = " + str(total_junc)
		#self.f.seek(0)
		
		print >>ROUT, 'pdf(\"%s\")' % (outfile + ".junction_plot.pdf")
		print >>ROUT, "events=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc)])+ ')'
		print >>ROUT, 'pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing events",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc))
		print >>ROUT, "dev.off()"
		
		print >>sys.stderr, "\n==================================================================="
		print >>sys.stderr, "Total splicing  Events:\t" + str(total_junc)
		print >>sys.stderr, "Known Splicing Events:\t" + str(known_junc)
		print >>sys.stderr, "Partial Novel Splicing Events:\t" + str(novel3or5_junc)
		print >>sys.stderr, "Novel Splicing Events:\t" + str(novel35_junc)
		
		#reset variables
		total_junc =0
		novel35_junc =0
		novel3or5_junc =0
		known_junc =0
		
		print >>OUT, "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation"
		for i in splicing_events:
			total_junc += 1
			(chrom, i_st, i_end) = i.split(":")
			print >>OUT, '\t'.join([chrom.replace("CHR","chr"),i_st,i_end]) + '\t' + str(splicing_events[i]) + '\t',
			i_st = int(i_st)
			i_end = int(i_end)
			if (refIntronStarts[chrom].has_key(i_st) and refIntronEnds[chrom].has_key(i_end)):
				print >>OUT, "annotated"
				known_junc +=1
			elif (not refIntronStarts[chrom].has_key(i_st) and not refIntronEnds[chrom].has_key(i_end)):
				print >>OUT, 'complete_novel'
				novel35_junc +=1
			else:
				print >>OUT, 'partial_novel'
				novel3or5_junc +=1
		
		if total_junc ==0:
			print >>sys.stderr, "No splice read found"
			sys.exit(1)
		print >>sys.stderr, "\nTotal splicing  Junctions:\t" + str(total_junc)
		print >>sys.stderr, "Known Splicing Junctions:\t" + str(known_junc)
		print >>sys.stderr, "Partial Novel Splicing Junctions:\t" + str(novel3or5_junc)
		print >>sys.stderr, "Novel Splicing Junctions:\t" + str(novel35_junc)
		print >>sys.stderr, "\n==================================================================="
		
		print >>ROUT, 'pdf("splicing_junction_pie.pdf")'
		print >>ROUT, "junction=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc,)])+ ')'
		print >>ROUT, 'pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing junctions",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc))
		print >>ROUT, "dev.off()"
		#print >>ROUT, "mat=matrix(c(events,junction),byrow=T,ncol=3)"
		#print >>ROUT, 'barplot(mat,beside=T,ylim=c(0,100),names=c("known","partial\nnovel","complete\nnovel"),legend.text=c("splicing events","splicing junction"),ylab="Percent")'

	def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1):
		'''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction'''
		
		out_file = outfile + ".junctionSaturation_plot.r"
		if refgene is None:
			print >>sys.stderr, "You must provide reference gene model in bed format."
			sys.exit(1)
		
		OUT = open(out_file,'w')


		#reading reference gene 
		knownSpliceSites= set()
		print >>sys.stderr, "reading reference bed file: ",refgene, " ... ",
		for line in open(refgene,'r'):
			if line.startswith(('#','track','browser')):continue  
			fields = line.split()
			if(len(fields)<12):
				print >>sys.stderr, "Invalid bed line (skipped):",line,
				continue
			chrom     = fields[0].upper()
			tx_start = int( fields[1] )
			tx_end   = int( fields[2] )
			if int(fields[9] ==1):
				continue    	
			
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			for st,end in zip (intron_start, intron_end):
				knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))
		print >>sys.stderr,"Done! Total "+str(len(knownSpliceSites)) + " known splicing junctions."


		#read SAM file
		samSpliceSites=[]
		intron_start=[]
		intron_end=[]
		uniqSpliceSites=collections.defaultdict(int)

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read

				if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:continue						#skip multiple map read		

				chrom = self.samfile.getrname(aligned_read.tid).upper()
				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)		
				intron_blocks = cigar.fetch_intron(chrom, hit_st, cigar_str)			
				if len(intron_blocks)==0:
					continue
				for intrn in intron_blocks:
					if intrn[2] - intrn[1] < min_intron:continue
					samSpliceSites.append(intrn[0] + ":" + str(intrn[1]) + "-" + str(intrn[2]))
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print >>sys.stderr, "shuffling alignments ...",
		random.shuffle(samSpliceSites)
		print >>sys.stderr, "Done"
				
		#resampling
		SR_num = len(samSpliceSites)
		sample_size=0
		all_junctionNum = 0	
		known_junc=[]
		all_junc=[]
		unknown_junc=[]
		#=========================sampling uniquely mapped reads from population
		tmp=range(sample_start,sample_end,sample_step)
		tmp.append(100)
		for pertl in tmp:	#[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
			knownSpliceSites_num = 0
			index_st = int(SR_num * ((pertl - sample_step)/100.0))
			index_end = int(SR_num * (pertl/100.0))
			if index_st < 0: index_st = 0
			sample_size += index_end -index_st
			
			print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(sample_size) + ") splicing reads.",
			
			#all splice juntion
			for i in range(index_st, index_end):
				uniqSpliceSites[samSpliceSites[i]] +=1	
			all_junctionNum = len(uniqSpliceSites.keys())
			all_junc.append(str(all_junctionNum))
			print >>sys.stderr, str(all_junctionNum) + " splicing junctions.",
			
			#known splice junction
			known_junctionNum = 0
			for sj in uniqSpliceSites:
				if sj in knownSpliceSites and uniqSpliceSites[sj] >= recur:
					known_junctionNum +=1
			print >>sys.stderr, str(known_junctionNum) + " known splicing junctions.",
			known_junc.append(str(known_junctionNum))
			
			#unknown splice junction
			unknown_junctionNum = 0
			for sj in uniqSpliceSites:
				if sj not in knownSpliceSites:
					unknown_junctionNum +=1
			unknown_junc.append(str(unknown_junctionNum))
			print >>sys.stderr, str(unknown_junctionNum) + " novel splicing junctions."
			
		#for j in uniq_SJ:
			#print >>OUT, j + "\t" + str(uniq_SJ[j])
		print >>OUT, "pdf(\'%s\')" % (outfile + '.junctionSaturation_plot.pdf')
		print >>OUT, "x=c(" + ','.join([str(i) for i in tmp]) + ')'
		print >>OUT, "y=c(" + ','.join(known_junc) + ')'
		print >>OUT, "z=c(" + ','.join(all_junc) + ')'
		print >>OUT, "w=c(" + ','.join(unknown_junc) + ')'
		print >>OUT, "m=max(%d,%d,%d)" % (int(int(known_junc[-1])/1000), int(int(all_junc[-1])/1000),int(int(unknown_junc[-1])/1000))
		print >>OUT, "n=min(%d,%d,%d)" % (int(int(known_junc[0])/1000), int(int(all_junc[0])/1000),int(int(unknown_junc[0])/1000))
		print >>OUT, "plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))"
		print >>OUT, "points(x,y/1000,type='o',col='red')"
		print >>OUT, "points(x,w/1000,type='o',col='green')"
		print >>OUT, 'legend(5,%d, legend=c("All junctions","known junctions", "novel junctions"),col=c("blue","red","green"),lwd=1,pch=1)' % int(int(all_junc[-1])/1000)
		print >>OUT, "dev.off()"

	def saturation_RPKM(self,refbed,outfile,sample_start=5,sample_step=5,sample_end=100,skip_multi=True, strand_rule=None):
		'''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		rpkm_file = outfile + ".eRPKM.xls"
		raw_file = outfile + ".rawCount.xls"
		
		RPKM_OUT = open(rpkm_file,'w')
		RAW_OUT = open(raw_file ,'w')
		
		ranges={}
		totalReads=0
		cUR_num = 0	#number of fragements
		cUR_plus = 0
		cUR_minus = 0
		block_list_plus = []	#non-spliced read AS IS, splicing reads were counted multiple times
		block_list_minus = []
		block_list = []
		strandRule = {}
				
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of: 'strand_rule' " +  strand_rule
			sys.exit(1)	


		#read SAM or BAM
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				
				if skip_multi:
					if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
						for i in aligned_read.tags:
							if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
								flag=1						#multiple hit read
								break
				if flag==1:continue						#skip multiple map read		
				
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				
				#determine read_id and read_strand
				if aligned_read.is_paired:						#pair end
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				else:read_id = ''								#single end
			
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'				
				strand_key = read_id + map_strand				#used to determine if a read should assign to gene(+) or gene(-)

				hit_st = aligned_read.pos
				cigar_str = cigar.list2str(aligned_read.cigar)		
				exon_blocks = cigar.fetch_exon(chrom, hit_st, cigar_str)	
				cUR_num += len(exon_blocks)
				
				#strand specific
				if strand_rule is not None:
					if strandRule[strand_key] == '+': cUR_plus += len(exon_blocks)
					if strandRule[strand_key] == '-': cUR_minus += len(exon_blocks)
					for exn in exon_blocks:
						if strandRule[strand_key] == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
						if strandRule[strand_key] == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
				#Not strand specific
				else:			
					for exn in exon_blocks:
						block_list.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
		except StopIteration:
			print >>sys.stderr, "Done"
		
		
		print >>sys.stderr, "shuffling alignments ...",
		random.shuffle(block_list_plus)
		random.shuffle(block_list_minus)
		random.shuffle(block_list)
		print >>sys.stderr, "Done"
		
		
		ranges_plus={}
		ranges_minus={}
		ranges={}
		sample_size=0
		RPKM_table=collections.defaultdict(list)
		rawCount_table=collections.defaultdict(list)
		RPKM_head=['#chr','start','end','name','score','strand']

		tmp=range(sample_start,sample_end,sample_step)
		tmp.append(100)
		#=========================sampling uniquely mapped reads from population
		for pertl in tmp:	#[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
			percent_st = (pertl-sample_step)/100.0
			percent_end = pertl/100.0
			if percent_st < 0: percent_st = 0
			sample_size = cUR_num * percent_end
			RPKM_head.append(str(pertl) + '%')
			
			if strand_rule is not None:
				print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(int(cUR_plus * percent_end)) + ") forward strand fragments ..."
				for i in block_list_plus[int(cUR_plus*percent_st):int(cUR_plus*percent_end)]:
					(chr,coord) = i.split(':')
					if chr not in ranges_plus:ranges_plus[chr] = Intersecter()
					else:ranges_plus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )								
				
				print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(int(cUR_minus * percent_end)) + ") reverse strand fragments ..."			
				for i in block_list_minus[int(cUR_minus*percent_st):int(cUR_minus*percent_end)]:
					(chr,coord) = i.split(':')
					if chr not in ranges_minus:ranges_minus[chr] = Intersecter()				
					else:ranges_minus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )						
			
			else:
				print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(int(sample_size)) + ") fragments ..."
				for i in block_list[int(cUR_num*percent_st):int(cUR_num*percent_end)]:
					(chr,coord) = i.split(':')
					if chr not in ranges:ranges[chr] = Intersecter()						
					else:ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )														

			#========================= calculating RPKM based on sub-population
			print >>sys.stderr, "assign reads to transcripts in " + refbed + ' ...'
			for line in open(refbed,'r'):
				try:
					if line.startswith(('#','track','browser')):continue  
            		# Parse fields from gene tabls
					fields = line.split()
					chrom     = fields[0].upper()
					tx_start  = int( fields[1] )
					tx_end    = int( fields[2] )
					geneName      = fields[3]
					strand    = fields[5]
					exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
					exon_starts = map((lambda x: x + tx_start ), exon_starts)
					exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
					exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
					exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
					key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
				except:
					print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
					continue
				mRNA_count=0	#we need to initializ it to 0 for each gene
				mRNA_len=sum(exon_sizes)
				for st,end in zip(exon_starts,exon_ends):
					#if chrom in ranges:
					if strand_rule is not None:
						if (strand == '+') and (chrom in ranges_plus): mRNA_count += len(ranges_plus[chrom].find(st,end))	
						if (strand == '-') and (chrom in ranges_minus): mRNA_count += len(ranges_minus[chrom].find(st,end))
					else:
						if chrom in ranges:
							mRNA_count += len(ranges[chrom].find(st,end))
				if mRNA_len ==0:
					print >>sys.stderr, geneName + " has 0 nucleotides. Exit!"
					sys.exit(1)
				if sample_size == 0:
					print >>sys.stderr, "Too few reads to sample. Exit!"
					sys.exit(1)
				mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * sample_size)
				RPKM_table[key].append(str(mRNA_RPKM))
				rawCount_table[key].append(str(mRNA_count))
			print >>sys.stderr, ""

		#self.f.seek(0)
		print >>RPKM_OUT, '\t'.join(RPKM_head)
		print >>RAW_OUT, '\t'.join(RPKM_head)
		for key in RPKM_table:
			print >>RPKM_OUT, key + '\t',
			print >>RPKM_OUT, '\t'.join(RPKM_table[key])
			print >>RAW_OUT, key + '\t',
			print >>RAW_OUT, '\t'.join(rawCount_table[key])		
		
	def fetchAlignments(self,chr,st,end):
		'''fetch alignment from sorted BAM file based on chr, st, end
		Note: BAM file must be indexed'''
		try:
			a=self.samfile.fetch(chr,st,end)
			return a
		except:
			return None
		

def print_bits_as_bed( bits ):
	end = 0
	while 1:
		start = bits.next_set( end )
		if start == bits.size: break
		end = bits.next_clear( start )
		print "%d\t%d" % ( start, end )
