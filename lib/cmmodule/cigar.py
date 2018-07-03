'''manipulate CIGAR string'''

#import built-in modules
import os,sys
import re
import string
import warnings
import string
import collections
import math
import sets
import random

head_clip = re.compile(r'^(\d+)S')
tail_clip = re.compile(r'(\d+)S$')
insertion = re.compile(r'(\d+)I')
deletion  = re.compile(r'(\d+)D')
matching  = re.compile(r'(\d+)M')
skipping  = re.compile(r'(\d+)N')
read_part = re.compile(r'(\d+)[MIS=X]')
ref_part = re.compile(r'(\d+)[MISND=X]')

prior_insertion = re.compile(r'(.+?)(\d+)I')
prior_deletion = re.compile(r'(.+?)(\d+)D')
prior_intron = re.compile(r'(.+?)(\d+)N')
prior_exon = re.compile(r'(.*?)(\d+)M')

def fetch_head_clip(chr, st, cigar):
	'''return genome coordinate of the head clip part encoded in cigar string
	NOTE: returned coordinates are 0-based.NOTE: st is 0-based'''
	
	block=[]
	chrom_end = int(st)		
	tmp = head_clip.findall(cigar)
	if len(tmp)==0:
		return block
	else:
		chrom_st = chrom_end - int(tmp[0])
		block.append([chr,chrom_st,chrom_end])
	return block

def fetch_tail_clip(chr, st, cigar):
	'''return genome coordinates of the tail clip part encoded in cigar string
	NOTE: returned coordinates are 0-based .  NOTE: st is 0-based'''
	
	block=[]
	h = head_clip.findall(cigar)
	t = tail_clip.findall(cigar)
	if len(t)==0:return block
	else:t_len = int(t[0])
	
	if len(h)==0:h_len=0
	else:h_len=int(h[0])
	ref_length = sum([int(i) for i in ref_part.findall(cigar)])
	#print read_length
	chrom_end = int(st) + (ref_length - h_len)		#because SAM is 1-based
	chrom_st = chrom_end - t_len
	block.append([chr,chrom_st,chrom_end])
	return block
	
def fetch_insertion(chr, st, cigar):
	'''return genome coordinates of the insertion (to reference) encoded in cigar string
	NOTE: returned coordinates are 0-based. Insertion to the reference.  NOTE: st is 0-based'''
	
	block=[]
	h = head_clip.findall(cigar)
	if len(h)==0:h_len=0
	else:h_len=int(h[0])
	
	ref_length=0
	m = prior_insertion.findall(cigar)
	if len(m)==0:
		return block
	else:
		for j in m:
			ref_length += sum([int(i) for i in ref_part.findall(j[0])])
			chrom_st = int(st) + (ref_length - h_len)
			chrom_end = chrom_st + int(j[1])
			block.append([chr,chrom_st,chrom_end])
	return block

def fetch_deletion(chr, st, cigar):
	'''return genome coordinates of the insertion (to reference) encoded in cigar string
	NOTE: returned coordinates are 0-based. Deletion to the reference.  NOTE: st is 0-based'''
	
	block=[]
	h = head_clip.findall(cigar)
	if len(h)==0:h_len=0
	else:h_len=int(h[0])
	
	ref_length=0
	m = prior_deletion.findall(cigar)
	if len(m)==0:
		return block
	else:
		for j in m:
			ref_length += sum([int(i) for i in ref_part.findall(j[0])])
			chrom_st = int(st) + (ref_length - h_len)
			chrom_end = chrom_st + int(j[1])
			block.append([chr,chrom_st,chrom_end])
			ref_length += int(j[1])
	return block

	
def fetch_intron(chr, st, cigar):
	'''return genome coordinates of the introns encoded in cigar string
	NOTE: returned coordinates are 0-based. Deletion to the reference NOTE: 
	st is 0-based '''
	
	block=[]
	h = head_clip.findall(cigar)
	if len(h)==0:h_len=0
	else:h_len=int(h[0])
	
	ref_length=0
	m = prior_intron.findall(cigar)
	if len(m)==0:
		return block
	else:
		for j in m:
			ref_length += sum([int(i) for i in ref_part.findall(j[0])])
			chrom_st = int(st)  + (ref_length - h_len)
			chrom_end = chrom_st + int(j[1])
			block.append([chr,chrom_st,chrom_end])
			ref_length += int(j[1])
	return block
	
def fetch_exon(chr, st, cigar):
	'''return genome coordinates of the exon encoded in cigar string
	NOTE: returned coordinates are 0-based. NOTE: st is 0-based'''
	
	block=[]
	h = head_clip.findall(cigar)
	if len(h)==0:h_len=0
	else:h_len=int(h[0])
	
	ref_length=0
	m = prior_exon.findall(cigar)
	if len(m)==0:
		return block
	else:
		for j in m:
			ref_length += sum([int(i) for i in ref_part.findall(j[0])])
			chrom_st = int(st) + (ref_length - h_len)
			chrom_end = chrom_st + int(j[1])
			block.append([chr,chrom_st,chrom_end])
			ref_length += int(j[1])
	return block
	
def list2str (lst):
	'''translate samtools returned cigar_list into cigar_string'''
	code2Char={'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
	
	cigar_str=''
	for i in lst:
		cigar_str += str(i[1]) + code2Char[str(i[0])]
	return cigar_str
