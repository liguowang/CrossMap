#!/usr/bin/env python
#Liguo Wang
#04/13/2011

#import built-in modules
import os,sys
import re
import string
#from optparse import OptionParser
import warnings
import string
import collections
import math

#import third-party modules
#from bx.bitset import *
#from bx.bitset_builders import *
#from bx.intervals import *
#import fasta
import bx.wiggle
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN


class ParseWig:
	'''provie methods to manipulate wiggle format file. For wiggle format see:
	http://genome.ucsc.edu/goldenPath/help/wiggle.html'''
	
	def __init__(self,wigFile):
		'''read wig file, creat wig obj'''
		self.scores = {}
		self.num_re=re.compile(r'[\d\.\-\+]+');
		fh=open(wigFile)
		#infile=open(wigFile,'r')
		for i, ( chrom, pos, val ) in enumerate( bx.wiggle.Reader( fh ) ):
			chrom=chrom.upper()
			if not chrom in self.scores: self.scores[ chrom ] = BinnedArray()
			self.scores[chrom][pos] = val
			if i % 100000 == 0: print("%i datapoints loaded \r" % i) 
		#print self.scores.keys()
		print("total " + str(i) + " points loaded")
		
	def fetch_all_scores(self,chr,st,end):
		'''fetch all wiggle scores defined by st and end.  NOTE:
		1)st and end are 0-based, half-open. (st,end]
		2)points without score are indicated as "nan"
		'''
		chr=chr.upper()
		return [ self.scores[chr][i] for i in range(st,end)]

	def fetch_max_scores(self,chr,st,end):
		''' fetch maximum score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''		
		
		chr=chr.upper()
		return max([ self.scores[chr][i] for i in range(st,end)])

	def fetch_min_scores(self,chr,st,end):
		''' fetch minimum score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''	
		
		chr=chr.upper()
		return min([ self.scores[chr][i] for i in range(st,end)])			
	
	def fetch_avg_scores(self,chr,st,end):
		''' fetch average score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''	
		
		chr=chr.upper()
		lst=[ float(self.scores[chr][i]) for i in range(st,end) if self.num_re.match(str(self.scores[chr][i]))]
		return sum(lst)/len(list(range(st,end)))
		
	def fetch_sum_scores(self,chr,st,end):
		''' fetch sum score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''	
		
		chr=chr.upper()
		lst=[ float(self.scores[chr][i]) for i in range(st,end) if self.num_re.match(str(self.scores[chr][i]))]
		return sum(lst)		
#if __name__ == "__main__": main()

class ParseWig2:
	'''provie methods to manipulate wiggle format file. For wiggle format see:
	http://genome.ucsc.edu/goldenPath/help/wiggle.html. The same coordinate could occur more than
	one time in wig file, and the scores will be sumed up. Slower than ParseWig'''
	
	def __init__(self,wigFile):
		'''read wig file, creat wig obj'''
		self.scores = {}
		self.num_re=re.compile(r'[\d\.\-\+]+');
		fh=open(wigFile)
		#infile=open(wigFile,'r')
		for i, ( chrom, pos, val ) in enumerate( bx.wiggle.Reader( fh ) ):
			chrom=chrom.upper()
			if not chrom in self.scores: self.scores[chrom] = BinnedArray()
			tmp=self.scores[chrom][pos]
			if isNaN(tmp):
				self.scores[chrom][pos] = val
			else:
				self.scores[chrom][pos] += val
			if i % 100000 == 0: print("%i datapoints loaded \r" % i) 
		#print self.scores.keys()
		print("total " + str(i) + " points loaded")
		
	def fetch_all_scores_by_range(self,chr,st,end):
		'''fetch all wiggle scores defined by st and end.  NOTE:
		1)st and end are 0-based, half-open. (st,end]
		2)points without score are indicated as "nan"'''
		chr=chr.upper()
		return [ self.scores[chr][i] for i in range(st,end)]
		
	def fetch_all_scores_by_positions(self,chr,lst):
		'''fetch all wiggle scores defined by st and end.  NOTE:
		2)points without score are indicated as "nan"'''
		chr=chr.upper()
		return [ self.scores[chr][i] for i in lst]
		
	def fetch_max_scores_by_range(self,chr,st,end):
		''' fetch maximum score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''			
		chr=chr.upper()
		return max([ self.scores[chr][i] for i in range(st,end)])

	def fetch_max_scores_by_positions(self,chr,lst):
		'''fetch maximum score defined by chr, st, end'''
		
		chr=chr.upper()
		return max([ self.scores[chr][i] for i in lst])

	def fetch_min_scores_by_range(self,chr,st,end):
		''' fetch minimum score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''			
		chr=chr.upper()
		return min([ self.scores[chr][i] for i in range(st,end)])			

	def fetch_min_scores_by_positions(self,chr,lst):
		''' fetch minimum score defined by chr, st, end
		'''			
		chr=chr.upper()
		return min([ self.scores[chr][i] for i in lst])
		
	
	def fetch_avg_scores_by_range(self,chr,st,end):
		''' fetch average score defined by chr, st, end
		1)st and end are 0-based, half-open. (st,end]
		'''			
		chr=chr.upper()
		lst=[ float(self.scores[chr][i]) for i in range(st,end) if self.num_re.match(str(self.scores[chr][i]))]
		return sum(lst)/len(list(range(st,end)))
		
	def fetch_avg_scores_by_positions(self,chr,lst):
		''' fetch average score defined by chr, st, end
		'''			
		chr=chr.upper()
		lst_score =[ float(self.scores[chr][i]) for i in lst if self.num_re.match(str(self.scores[chr][i]))]
		return sum(lst_score)/len(lst_score)
		
	def fetch_sum_scores_by_range(self,chr,st,end):
		''' fetch sum score defined by chr, st, end
		'''			
		chr=chr.upper()
		lst=[ float(self.scores[chr][i]) for i in range(st,end) if self.num_re.match(str(self.scores[chr][i]))]
		return sum(lst)		

	def fetch_sum_scores_by_positions(self,chr,lst):
		''' fetch sum score defined by chr, st, end
		'''			
		chr=chr.upper()
		lst_score=[ float(self.scores[chr][i]) for i in lst if self.num_re.match(str(self.scores[chr][i]))]
		return sum(lst_score)

	def distriub_wig(self,bed,till_count=100):
		'''calculate coverage over bed file (only consider exon regions). The mRNA sequences in input
		bed file will be cut into 100 tills of equal size'''
		
		print("Reading " + bed + " ...", file=sys.stderr)
		for line in open(bed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue
				fields=line.rstrip('\r\n').split()
				txStart=int(fields[1])
				chrom=fields[0]
				strand=fields[5]
				geneName=fields[3]
				score=fields[4]
				exon_start=list(map(int,fields[11].rstrip(',').split(',')))
				exon_start=list(map((lambda x: x + txStart),exon_start))
				exon_end=list(map(int,fields[10].rstrip(',').split(',')))
				exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
			except:
				print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
				continue
		
#if __name__ == "__main__": main()
