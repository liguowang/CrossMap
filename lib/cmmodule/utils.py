import sys
from time import strftime
import pyBigWig
from bx.intervals.intersection import Interval, Intersecter

from cmmodule  import ireader


def printlog (mesg_lst):
	'''
	Print messages into stderr.
	'''
	msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + ' '.join(mesg_lst)
	print(msg, file=sys.stderr)

def parse_header( line ):
	'''
	Parse the header line of wiggle file.
	'''
	return dict( [ field.split( '=' ) for field in line.split()[1:] ] )

def update_chromID(c_temp, c_target):
	'''
	Update chromsome ID to to the template

	Parameters
	----------
	c_temp : str
		Template of chromsome ID

	c_target : str
		Chromosome ID that need to be updated

	Returns
	--------
	Updated chromosome ID

	Examples
	--------
	>>> update_chromID('chrX',1)
	'chr1'
	>>> update_chromID('1','chrY')
	'Y'
	'''
	c_temp = str(c_temp)
	c_target = str(c_target)
	if c_temp.startswith('chr'):
		if c_target.startswith('chr'):
			return c_target
		else:
			return ('chr' + c_target)
	else:
		if c_target.startswith('chr'):
			return c_target.replace('chr','')
		else:
			return c_target

def revcomp_DNA(dna, extended=True):
	'''
	Reverse complement of input DNA sequence.

	Parameters
	----------
	dna : str
		DNA sequences made of 'A', 'C', 'G', 'T', 'N' or 'X'

	extended : bool
		Support full IUPAC nucleotides.

	Examples
	--------
	>>> revcomp_DNA('AACGTG')
	'CACGTT'
	'''
	if extended:
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y': 'R', 'R': 'Y', 'S': 'W', 'W': 'S', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', '.':'.', '*':'*'}
	else:
		complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}

	seq = dna.replace(' ','').upper()
	if ',' not in seq:
		return ''.join([complement[base] for base in reversed(seq)])
	else:
		seqs = seq.split(',')
		comp_seqs = []
		for s in seqs:
			comp_seqs.append(''.join([complement[base] for base in reversed(s)]))
		return ','.join(comp_seqs)

def wiggleReader( f ):
	'''
	Read wiggle (http://genome.ucsc.edu/goldenPath/help/wiggle) file of different styles.

	Parameters
	----------
	f : file
		file in wiggle format. Can be fixedStep, variableStep, or bed4

	Yields
	------
	chrom, start, end, strand, score
	'''
	current_chrom = None
	current_pos = None
	current_step = None

	# always for wiggle data
	strand = '+'

	mode = "bed"
	for line in ireader.reader(f):
		if line.isspace() or line.startswith(("track", "#", "browser")):
			continue
		elif line.startswith( "variableStep" ):
			header = parse_header( line )
			current_chrom = header['chrom']
			current_pos = None
			current_step = None
			if 'span' in header: current_span = int( header['span'] )
			else: current_span = 1
			mode = "variableStep"
		elif line.startswith( "fixedStep" ):
			header = parse_header( line )
			current_chrom = header['chrom']
			current_pos = int( header['start'] ) - 1
			current_step = int( header['step'] )
			if 'span' in header: current_span = int( header['span'] )
			else: current_span = 1
			mode = "fixedStep"
		elif mode == "bed":
			fields = line.split()
			if len( fields ) > 3:
				if len( fields ) > 5:
					yield fields[0], int( fields[1] ), int( fields[2] ), fields[5], float( fields[3] )
				else:
					yield fields[0], int( fields[1] ), int( fields[2] ), strand, float( fields[3] )
		elif mode == "variableStep":
			fields = line.split()
			pos = int( fields[0] ) - 1
			yield current_chrom, pos, pos + current_span, strand, float( fields[1] )
		elif mode == "fixedStep":
			yield current_chrom, current_pos, current_pos + current_span, strand, float( line.split()[0] )
			current_pos += current_step
		else:
			raise "Unexpected input line: %s" % line.strip()

def bigwigReader(infile):
	'''
	Read bigwig (https://genome.ucsc.edu/goldenPath/help/bigWig.html) files.

	Parameters
	----------
	infile: file
		Bigwig format file

	Yields
	------
	chrom, start, end, score

	'''
	bw = pyBigWig.open(infile)
	chrom_sizes = bw.chroms()
	for chr_name, chr_size in list(chrom_sizes.items()):
		for i_st, i_end, i_score in bw.intervals(chr_name, 0, chr_size):
			yield(( chr_name, i_st, i_end, i_score ))

def check_bed12(bedline):
	'''
	Check if bed12 format is correct or not.

	Parameters
	----------
	bedline : str
		line in BED format.

	'''
	fields = bedline.strip().split()
	if len(fields) !=12:
		return False
	if fields[5] not in ['+','-','.']:
		return False
	try:
		chromStart = int(fields[1])
		chromEnd = int(fields[2])
		thickStart = int(fields[6])
		thickEnd = int(fields[7])
		blockCount = int(fields[9])
		blockSizes =  [int(i) for i in fields[10].rstrip(',').split(',')]
		blockStarts = [int(i) for i in fields[11].rstrip(',').split(',')]
	except:
		return False
	if chromStart > chromEnd or thickStart > thickEnd:
		return False
	if thickStart < chromStart or thickEnd > chromEnd:
		return False
	if len(blockSizes) != blockCount:
		return False
	if len(blockStarts) != blockCount:
		return False
	if blockCount <1:
		return False
	for i in blockSizes:
		if i < 0: return False
	for i in blockStarts:
		if i < 0: return False
	return True

def intersectBed(lst1, lst2):
	'''
	Return intersection of two bed regions.

	Parameters
	----------
	lst1 : list
		The 1st genomic region. List of chrom, start, end.
		Example: ['chr1',10, 100]

	lst2 : list
		 The 2nd genomic region. List of chrom, start, end.
		 Example: ['chr1',50, 120]

	Examples
	--------
	>>> intersectBed(['chr1',10, 100],['chr1',50, 120])
	('chr1', 50, 100)
	>>> intersectBed(['chr1',10, 100],['chr1',20, 30])
	('chr1', 20, 30)

	'''
	(chr1, st1, end1) = lst1
	(chr2, st2, end2) = lst2
	if int(st1) > int(end1) or int(st2) > int(end2):
		raise Exception ("Start cannot be larger than end")
	if chr1 != chr2:
		return None
	if int(st1) > int(end2) or int(end1) < int(st2):
		return None
	return (chr1, max(st1, st2), min(end1,end2))

def read_chain_file (chain_file, print_table = False):
	'''
	Read chain file.

	Parameters
	----------
	chain_file : file
		Chain format file. Input chain_file could be either plain text, compressed file
		(".gz",".Z", ".z", ".bz", ".bz2", ".bzip2"), or a URL pointing to the chain file
		("http://","https://", "ftp://"). If url was used, chain file must be plain text.

	print_table : bool, optional
		Print mappings in human readable table.

	Returns
	-------
	maps : dict
		Dictionary with source chrom name as key, IntervalTree object as value. An
		IntervalTree contains many intervals. An interval is a start and end position
		and a value. eg. Interval(11, 12, strand="-", value = "abc")

	target_chromSize : dict
		Chromosome sizes of target genome

	source_chromSize : dict
		Chromosome sizes of source genome
	'''

	printlog(["Read the chain file: ", chain_file])
	maps={}
	target_chromSize={}
	source_chromSize={}
	if print_table:
		blocks=[]

	for line in ireader.reader(chain_file):
		# Example: chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
		if not line.strip():
			continue
		line=line.strip()
		if line.startswith(('#',' ')):continue
		fields = line.split()

		if fields[0] == 'chain' and len(fields) in [12, 13]:
			#score = int(fields[1])		  # Alignment score
			source_name = fields[2]		  # E.g. chrY
			source_size = int(fields[3])  # Full length of the chromosome
			source_strand = fields[4]	  # Must be +
			if source_strand != '+':
				raise Exception("Source strand in a chain file must be +. (%s)" % line)
			source_start = int(fields[5]) # Start of source region
			#source_end = int(fields[6])	  # End of source region

			target_name = fields[7]		  # E.g. chr5
			target_size = int(fields[8])  # Full length of the chromosome
			target_strand = fields[9]	  # + or -
			target_start = int(fields[10])
			#target_end = int(fields[11])
			target_chromSize[target_name]= target_size
			source_chromSize[source_name] = source_size

			if target_strand not in ['+', '-']:
				raise Exception("Target strand must be - or +. (%s)" % line)
			#chain_id = None if len(fields) == 12 else fields[12]
			if source_name not in maps:
				maps[source_name] = Intersecter()

			sfrom, tfrom = source_start, target_start

		# Now read the alignment chain from the file and store it as a list (source_from, source_to) -> (target_from, target_to)
		elif fields[0] != 'chain' and len(fields) == 3:
			size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
			if print_table:
				if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
				elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))

			if target_strand == '+':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
			elif  target_strand == '-':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))

			sfrom += size + sgap
			tfrom += size + tgap

		elif fields[0] != 'chain' and len(fields) == 1:
			size = int(fields[0])
			if print_table:
				if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
				elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))

			if target_strand == '+':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
			elif target_strand == '-':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
		else:
			raise Exception("Invalid chain format. (%s)" % line)
	#if (sfrom + size) != source_end  or (tfrom + size) != target_end:
	#	 raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)

	if print_table:
		for i in blocks:
			print('\t'.join([str(n) for n in i]))

	return (maps,target_chromSize, source_chromSize)


def map_coordinates(mapping, q_chr, q_start, q_end, q_strand = '+', print_match = False):
	'''
	Map coordinates from source (i.e. original) assembly to target (i.e. new) assembly.

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	q_chr : str
		Chromosome ID of query interval

	q_start : int
		Start position of query interval.

	q_end : int
		End position of query interval.

	q_strand : str
		Strand of query interval.

	print_match : bool
		Print match table.
	'''

	matches = []
	complement = {'+':'-','-':'+'}

	if q_chr in mapping:
		targets = mapping[q_chr].find(q_start, q_end)
	elif q_chr.replace('chr','') in mapping:
		targets = mapping[q_chr.replace('chr','')].find(q_start, q_end)
	elif ('chr' + q_chr) in mapping:
		targets = mapping['chr' + q_chr].find(q_start, q_end)
	else:
		return None
	if len(targets)==0:
		return None
	elif len(targets)==1:
		s_start = targets[0].start
		s_end = targets[0].end
		t_chrom = targets[0].value[0]
		t_chrom = update_chromID(q_chr, t_chrom)
		t_start = targets[0].value[1]
		t_end = targets[0].value[2]
		t_strand = targets[0].value[3]

		(chr, real_start, real_end)	 = intersectBed((q_chr,q_start,q_end),(q_chr,s_start,s_end))
		l_offset = abs(real_start - s_start)
		#r_offset = real_end - s_end
		size = abs(real_end - real_start)

		matches.append( (chr, real_start, real_end,q_strand))
		if t_strand == '+':
			i_start = t_start + l_offset
			if q_strand == '+':
				matches.append( (t_chrom, i_start, i_start + size, t_strand))
			else:
				matches.append( (t_chrom, i_start, i_start + size, complement[t_strand]))
		elif t_strand == '-':
			i_start = t_end - l_offset - size
			if q_strand == '+':
				matches.append( (t_chrom, i_start,	i_start + size, t_strand))
			else:
				matches.append( (t_chrom, i_start,	i_start + size, complement[t_strand]))
		else:
			raise Exception("Unknown strand: %s. Can only be '+' or '-'." % q_strand)

	elif len(targets) > 1:
		for t in targets:
			s_start = t.start
			s_end = t.end
			t_chrom = t.value[0]
			t_chrom = update_chromID(q_chr, t_chrom)
			t_start = t.value[1]
			t_end = t.value[2]
			t_strand = t.value[3]

			(chr, real_start, real_end)	 = intersectBed((q_chr,q_start,q_end),(q_chr,s_start,s_end))

			l_offset = abs(real_start - s_start)
			#r_offset = abs(real_end - s_end)
			size = abs(real_end - real_start)
			matches.append( (chr, real_start, real_end,q_strand) )
			if t_strand == '+':
				i_start = t_start + l_offset
				if q_strand == '+':
					matches.append( (t_chrom, i_start, i_start + size, t_strand))
				else:
					matches.append( (t_chrom, i_start, i_start + size, complement[t_strand]))
			elif t_strand == '-':
				i_start = t_end - l_offset - size
				if q_strand == '+':
					matches.append( (t_chrom, i_start,	i_start + size, t_strand))
				else:
					matches.append( (t_chrom, i_start,	i_start + size, complement[t_strand]))
			else:
				raise Exception("Unknown strand: %s. Can only be '+' or '-'." % q_strand)

	if print_match:
		print(matches)
		# input: 'chr1',246974830,247024835
		# output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
		# [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]

	return matches
