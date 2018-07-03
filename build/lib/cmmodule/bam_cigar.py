#!/usr/bin/env python

'''manipulate CIGAR string represented as list of tuples (BAM file)
BAM	OP	Description
0	M	alignment match
1	I	insertion to reference
2	D	deletion from reference. region deleted from reference genome
3	N	skipped region from the reference
4	S	soft clipping (clipped sequence present in SEQ)
5	H	hard clipping (clipped sequences NOT present in SEQ) 
6	P	padding (silent deletion from padded reference)
7	=	sequence match
8	X	sequence mismatch

Example:
The tuple [ (0,3), (1,5), (0,2) ] refers to an alignment with 3 matches, 
5 insertions and another 2 matches.
 
NOTE: only deal with Match, Gap, Soft Clip, Insertion, Deletion 
'''

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012 Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Development" #Prototype or Production


def fetch_exon(chrom, st, cigar):
	''' fetch exon regions defined by cigar. st must be zero based
	return list of tuple of (chrom,st, end)
	'''
	#match = re.compile(r'(\d+)(\D)')
	chrom_st = st
	exon_bound =[]
	for c,s in cigar:	#code and size
		if c==0:		#match
			exon_bound.append((chrom, chrom_st,chrom_st + s))
			chrom_st += s
		elif c==1:		#insertion to ref
			continue
		elif c==2:		#deletion to ref
			chrom_st += s
		elif c==3:		#gap or intron
			chrom_st += s
		elif c==4:		#soft clipping. We do NOT include soft clip as part of exon
			chrom_st += s
		else:
			continue
	return exon_bound
	
def fetch_intron(chrom, st, cigar):
	''' fetch intron regions defined by cigar. st must be zero based
	return list of tuple of (chrom,st, end)
	'''
	#match = re.compile(r'(\d+)(\D)')
	chrom_st = st
	intron_bound =[]
	for c,s in cigar:	#code and size
		if c==0:		#match
			chrom_st += s
		elif c==1:		#insertion to ref
			continue
		elif c==2:		#deletion to ref
			chrom_st += s
		elif c==3:		#gap or intron
			intron_bound.append((chrom, chrom_st,chrom_st+s))
			chrom_st += s
		elif c==4:		#soft clipping. We do NOT include soft clip as part of exon
			chrom_st += s
		else:
			continue
	return intron_bound	

def fetch_clip(chrom, st, cigar):
	''' fetch head soft clip regions defined by cigar. st must be zero based
	return list of tuple of (chrom,st, end)
	'''
	#match = re.compile(r'(\d+)(\D)')
	chrom_st = st
	clip_bound =[]
	for c,s in cigar:	#code and size
		if c==0:		#match
			chrom_st += s
		elif c==1:		#insertion to ref
			continue
		elif c==2:		#deletion to ref
			chrom_st += s
		elif c==3:		#gap or intron
			chrom_st += s
		elif c==4:		#soft clipping. We do NOT include soft clip as part of exon
			clip_bound.append((chrom, chrom_st, chrom_st + s))
			chrom_st += s
		else:
			continue
	return clip_bound		

def fetch_deletion(chrom, st, cigar):
	''' fetch deletion regions defined by cigar. st must be zero based
	return list of tuple of (chrom,st, end)
	'''
	#match = re.compile(r'(\d+)(\D)')
	chrom_st = st
	del_bound =[]
	for c,s in cigar:	#code and size
		if c==0:		#match
			chrom_st += s
		elif c==1:		#insertion to ref
			continue
		elif c==2:		#deletion to ref
			del_bound.append((chrom, chrom_st, chrom_st + s))
			chrom_st += s
		elif c==3:		#gap or intron
			chrom_st += s
		elif c==4:		#soft clipping. We do NOT include soft clip as part of exon
			chrom_st += s
		else:
			continue
	return del_bound			
	
def fetch_insertion(chrom, st, cigar):
	''' fetch insertion regions defined by cigar. st must be zero based
	return list of tuple of (chrom,st, end)
	
	NOTE: insertion region does NOT present in reference genome and there
	fore cannot represented using reference coordinates[start, end]. So we
	use [start, SIZE).
	[(100,2)] means 2nt insert between 100 and 101
	'''
	#match = re.compile(r'(\d+)(\D)')
	chrom_st = st
	ins_bound =[]
	for c,s in cigar:	#code and size
		if c==0:		#match
			chrom_st += s
		elif c==1:		#insertion to ref
			ins_bound.append((chrom, chrom_st, s))
			continue
		elif c==2:		#deletion to ref
			chrom_st += s
		elif c==3:		#gap or intron
			chrom_st += s
		elif c==4:		#soft clipping. We do NOT include soft clip as part of exon
			chrom_st += s
		else:
			continue
	return ins_bound			
def list2str (lst):
	'''translate samtools returned cigar_list into cigar_string'''
	code2Char={'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
        
	cigar_str=''
	for i in lst:
		cigar_str += str(i[1]) + code2Char[str(i[0])]
	return cigar_str
	