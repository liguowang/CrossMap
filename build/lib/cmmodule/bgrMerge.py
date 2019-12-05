import os,sys
import subprocess
import random
import string
from time import strftime

from cmmodule  import wig_reader
from cmmodule  import myutils

def randomword(length):
	return ''.join(random.choice(string.ascii_uppercase) for _ in range(length))

def printlog (mesg_lst):
	'''print progress into stderr'''
	if len(mesg_lst)==1:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " +  mesg_lst[0]
	else:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + ' '.join(mesg_lst)
	print(msg, file=sys.stderr)

def read_bed_by_chr(f):
	'''input bed file'''
	
	# sort input file
	tmp_file_name = randomword(10)
	printlog (["Sorting bedGraph file:" + f])
	TMP = open(tmp_file_name,'w')
	sort_cmd = myutils.which('sort')
	try:
		subprocess.call([sort_cmd, "-k1,1", "-k2,2n", f], stdout= TMP)
	except:
		raise Exception("Cannot find GNU \"sort\" command")
	TMP.close()
	
	# generate list
	ret_list=[]
	line_num = 0
	for line in open(tmp_file_name):
		if not line.strip():
			continue
		line_num += 1
		line = line.strip()
		
		if line_num == 1:
			chrom = line.split()[0]
			ret_list.append(line)
			continue
		
		if line.split()[0] != chrom:
			yield ret_list
			ret_list=[]
			chrom = line.split()[0]
		ret_list.append(line)
	yield ret_list
	os.remove(tmp_file_name)

def merge(infile):
	'''merge bedGraph format file. Signal value will be accumulated'''

	line_iter =  read_bed_by_chr(infile)
	for lines in line_iter:
		top_marker = 0
		overlap_pos2val = {}
		for i in range(0, len(lines)-1):
			(chr, start, end, score) = lines[i].split()
			start = int(start)
			end = int(end)
			score = float(score)
			if start < 0 or end < 0 :
				continue
					
			if start  >= top_marker and end <= int(lines[i+1].split()[1]):
				if len(overlap_pos2val) !=0:
					for m,n,p in wig_reader.wig_to_bgr2(overlap_pos2val):
						yield((chr, m, n, p))			
				yield((chr, start, end, score))
				overlap_pos2val = {}
			else:
				for ind in range(start+1, end +1):
					if ind in overlap_pos2val:
						overlap_pos2val[ind] += score
					else:
						overlap_pos2val[ind] = score
			top_marker = max(top_marker, end)
		
		# deal with last line
		else:
			(last_chr, last_start, last_end, last_score) = lines[-1].split()
			try:
				last_start = int(last_start)
				last_end = int(last_end)
				last_score = float(last_score)
			except:
				print(last_chr, last_start,last_end)
				pass
	
			if last_start >= top_marker:
				if len(overlap_pos2val) !=0:
					for m,n,p in wig_reader.wig_to_bgr2(overlap_pos2val):
						yield((last_chr, m, n, p))
				yield( (last_chr, last_start, last_end, last_score ))
				overlap_pos2val = {}
			else:
				for ind in range(last_start+1, last_end +1):
					if ind in overlap_pos2val:
						overlap_pos2val[ind] += last_score
					else:
						overlap_pos2val[ind] = last_score
		
			if len(overlap_pos2val) !=0:
				for m,n,p in wig_reader.wig_to_bgr2(overlap_pos2val):
					yield((last_chr, m, n, p))
			
