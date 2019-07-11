#!/usr/bin/env python

"""
Read a wiggle or bigwig track and print out a series of lines containing
"chrom position score". Ignores track lines, handles bed, variableStep
and fixedStep wiggle lines.
"""
import sys
import bx.wiggle
from bx.bbi.bigwig_file import BigWigFile
import numpy
import collections
from itertools import groupby
from operator import itemgetter
from cmmodule import BED

def wig_to_bgr2(pos2val):
	'''pos2val is dictionary: position: value. position is 0 based '''
	v2p = collections.defaultdict(list) 
	point_num = 1
	count = 0
	coord = min(pos2val.keys())
	range2p={}      #coorindate range to value, bedgraph. #[start]=[len,value]
	for v in list(pos2val.values()):
		coord +=1
		#if v != 0: print >>OUT, "%d\t%.2f" % (coord,v)
		if v != 0: v2p[v].append(coord)
	for v in v2p:
		for k,g in groupby(enumerate(v2p[v]), lambda i_x:i_x[0]-i_x[1]):
			for l in [list(map(itemgetter(1), g))]:
				range2p[l[0]-1] = [len(l),v]
	for i in sorted(range2p):
		yield((i-1, i + range2p[i][0]-1, range2p[i][1]))

def wig_to_bgr1(pos2val):
	'''pos2val is dictionary: position: value. position is 0 based '''
	point_num = 1
	count = 0
	for pos in sorted(pos2val):
		count += 1
		if count ==1:	#initilize value. only run one time
			up_bound = pos+1
			score = pos2val[pos]
			continue
		if pos2val[pos] == score:
			point_num += 1
			up_bound = pos+1
		else:
			yield((up_bound - point_num-1, up_bound-1, score))
			score = pos2val[pos]
			up_bound = pos +1
			point_num = 1
		
		
		
def wig_reader(infile, chrom_sizes=None, informat='wiggle', bin_size = 2000):
	'''infile: either a wiggle or bigwig format file
	   chromsize: chrom_name: size, only needed is format is bigwig
	   format: either 'wiggle' or 'bigwig'
	   return: chrom, position (0-based), value
	'''
	if informat.upper()=='WIGGLE':
		point_num = 1
		count = 0
		for chrom, start, end, strand, score in bx.wiggle.IntervalReader( infile ):
			yield (chrom, start, end, score)
			"""
			count += 1
			if count ==1:
				chrom = fields[0]
				up_bound = fields[1]+1
				score = fields[2]
				continue
			if (fields[0] == chrom) and (fields[1] +1 == up_bound + 1) and (fields[2] == score):
				point_num += 1
				up_bound = fields[1]+1
				continue
			else:
				yield((chrom, up_bound - point_num, up_bound, score))
				chrom = fields[0]
				score = fields[2]
				up_bound = fields[1]+1
				point_num = 1
			"""
			
	
	elif informat.upper() == 'BIGWIG':
		bw_obj = BigWigFile(file=open(infile))
		for chr_name, chr_size in list(chrom_sizes.items()):
			for chrom, st, end in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = bin_size):
				sig_list = bw_obj.get_as_array(chrom,st,end)
				if sig_list is None:
					continue
				sig_list = numpy.nan_to_num(sig_list)
				if numpy.sum(sig_list)==0:
					continue	
				low_bound = st
				point_num = 1
				score = sig_list[0]
				for value in (sig_list[1:]):
					if value == score:
						point_num += 1
					else:
						yield(( chrom, low_bound,  low_bound + point_num, score))
						score = value
						low_bound = low_bound + point_num
						point_num = 1
	else:
		raise Exception("Unknown format. Must be 'wiggle' or 'bigwig'")
	
			
if __name__== '__main__':
	
	# test bigwig input

	#chromSize_hg18={
	#"chr1" : 247249719,
	#}
	
	#for a in wig_reader(infile='GSM686950_reverse.hg18.bw', chrom_sizes = chromSize_hg18, informat='bigwig'):
	#	print '\t'.join(map(str, a))
		
	# test wig input
	for a in wig_reader(infile=sys.argv[1], informat='wiggle'):
		print('\t'.join(map(str, a)))
	
