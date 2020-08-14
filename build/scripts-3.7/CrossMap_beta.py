#!python
'''
---------------------------------------------------------------------------------------
CrossMap: lift over genomic coordinates between genome assemblies.
Supports BED/BedGraph, GFF/GTF, BAM/SAM/CRAM, BigWig/Wig, VCF, and MAF format files.
---------------------------------------------------------------------------------------
'''

import os,sys
import re
import optparse
from textwrap import wrap
from time import strftime

import pyBigWig
import pysam
from bx.intervals.intersection import Interval, Intersecter
import datetime
from cmmodule  import ireader
from cmmodule  import annoGene
from cmmodule  import sam_header
from cmmodule  import bgrMerge

__author__ = "Liguo Wang, Hao Zhao"
__contributor__="Liguo Wang, Hao Zhao"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPLv2"
__version__="0.4.3"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"

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

def cal_mapped_size(l):
	"""
	Calculate the mapped length of query
	"""
	mapped_size = 0.0
	t_chroms = set()
	t_starts = []
	t_ends = []
	for j in range(1,len(l),2):

		t_chroms.add(l[j][0])
		t_starts.append(l[j][1])
		t_ends.append(l[j][2])

		#q_chrom = l[j-1][0]
		q_start = l[j-1][1]
		q_end = l[j-1][2]
		mapped_size += (q_end - q_start)
	return (mapped_size, t_chroms, min(t_starts), max(t_ends))

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

def crossmap_vcf_file(mapping, infile, outfile, liftoverfile, refgenome):
	'''
	Convert genome coordinates in VCF format.

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	infile : file
		Input file in VCF format. Can be a regular or compressed (*.gz, *.Z,*.z, *.bz,
		*.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to
		remote file.

	outfile : str
		prefix of output files.

	liftoverfile : file
		Chain (https://genome.ucsc.edu/goldenPath/help/chain.html) format file. Can be a
		regular or compressed (*.gz, *.Z,*.z, *.bz, *.bz2, *.bzip2) file, local file or
		URL (http://, https://, ftp://) pointing to remote file.

	refgenome : file
		The genome sequence file of 'target' assembly in FASTA format.
	'''

	#index refegenome file if it hasn't been done
	if not os.path.exists(refgenome + '.fai'):
		printlog(["Creating index for", refgenome])
		pysam.faidx(refgenome)

	refFasta = pysam.Fastafile(refgenome)

	FILE_OUT = open(outfile ,'w')
	UNMAP = open(outfile + '.unmap','w')

	total = 0
	fail = 0
	withChr = False # check if the VCF data lines use 'chr1' or '1'

	for line in ireader.reader(infile):
		if not line.strip():
			continue
		line=line.strip()

		#deal with meta-information lines.
		#meta-information lines needed in both mapped and unmapped files
		if line.startswith('##fileformat'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##INFO'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##FILTER'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##FORMAT'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##ALT'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##SAMPLE'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##PEDIGREE'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)

		#meta-information lines needed in unmapped files
		elif line.startswith('##assembly'):
			print(line, file=UNMAP)
		elif line.startswith('##contig'):
			print(line, file=UNMAP)
			if 'ID=chr' in line:
				withChr = True

		#update contig information
		elif line.startswith('#CHROM'):
			printlog(["Updating contig field ... "])
			target_gsize = dict(list(zip(refFasta.references, refFasta.lengths)))
			for chr_id in sorted(target_gsize):
				if chr_id.startswith('chr'):
					if withChr is True:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % (chr_id, target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)
					else:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % (chr_id.replace('chr',''), target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)
				else:
					if withChr is True:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % ('chr' + chr_id, target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)
					else:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % (chr_id, target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)

			print("##liftOverProgram=CrossMap(https://sourceforge.net/projects/crossmap/)", file=FILE_OUT)
			print("##liftOverFile=" + liftoverfile, file=FILE_OUT)
			print("##new_reference_genome=" + refgenome, file=FILE_OUT)
			print("##liftOverTime=" + datetime.date.today().strftime("%B%d,%Y"), file=FILE_OUT)
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
			printlog(["Lifting over ... "])

		else:
			if line.startswith('#'):continue
			fields = str.split(line,maxsplit=7)
			total += 1

			chrom = fields[0]
			start = int(fields[1])-1	 # 0 based
			end = start + len(fields[3])

			a = map_coordinates(mapping, chrom, start, end,'+')
			if a is None:
				print (line + "\tFail(Unmap)", file=UNMAP)
				fail += 1
				continue

			if len(a) == 2:
				# update chrom
				target_chr = str(a[1][0])	#target_chr is from chain file, could be 'chr1' or '1'
				target_start = a[1][1]
				target_end = a[1][2]
				fields[0] = target_chr

				# update start coordinate
				fields[1] = target_start + 1

				# update ref allele
				target_chr = update_chromID(refFasta.references[0], target_chr)
				fields[3] = refFasta.fetch(target_chr,target_start,target_end).upper()

				if a[1][3] == '-':
					fields[4] = revcomp_DNA(fields[4], True)

				if fields[3] != fields[4]:
					print('\t'.join(map(str, fields)), file=FILE_OUT)
				else:
				   print (line + "\tFail(REF==ALT)", file=UNMAP)
				   fail += 1
			else:
				print (line + "\tFail(Multiple_hits)", file=UNMAP)
				fail += 1
				continue
	FILE_OUT.close()
	UNMAP.close()
	printlog (["Total entries:", str(total)])
	printlog (["Failed to map:", str(fail)])



def crossmap_gvcf_file(mapping, infile, outfile, liftoverfile, refgenome):
	'''
	Convert genome coordinates in GVCF format.

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	infile : file
		Input file in GVCF format. Can be a regular or compressed (*.gz, *.Z,*.z, *.bz,
		*.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to
		remote file.

	outfile : str
		prefix of output files.

	liftoverfile : file
		Chain (https://genome.ucsc.edu/goldenPath/help/chain.html) format file. Can be a
		regular or compressed (*.gz, *.Z,*.z, *.bz, *.bz2, *.bzip2) file, local file or
		URL (http://, https://, ftp://) pointing to remote file.

	refgenome : file
		The genome sequence file of 'target' assembly in FASTA format.
	'''

	#index refegenome file if it hasn't been done
	if not os.path.exists(refgenome + '.fai'):
		printlog(["Creating index for", refgenome])
		pysam.faidx(refgenome)

	refFasta = pysam.Fastafile(refgenome)

	FILE_OUT = open(outfile ,'w')
	UNMAP = open(outfile + '.unmap','w')

	total_var = 0
	failed_var = 0
	total_region = 0
	failed_region = 0
	withChr = False # check if the VCF data lines use 'chr1' or '1'

	for line in ireader.reader(infile):
		if not line.strip():
			continue
		line=line.strip()

		#deal with meta-information lines.
		#meta-information lines needed in both mapped and unmapped files
		if line.startswith('##fileformat'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##INFO'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##FILTER'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##FORMAT'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##ALT'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##SAMPLE'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##PEDIGREE'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##GVCFBlock'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##GATKCommandLine'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
		elif line.startswith('##source'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)

		#meta-information lines needed in unmapped files
		elif line.startswith('##assembly'):
			print(line, file=UNMAP)
		elif line.startswith('##contig'):
			print(line, file=UNMAP)
			if 'ID=chr' in line:
				withChr = True

		#update contig information
		elif line.startswith('#CHROM'):
			printlog(["Updating contig field ... "])
			target_gsize = dict(list(zip(refFasta.references, refFasta.lengths)))
			for chr_id in sorted(target_gsize):
				if chr_id.startswith('chr'):
					if withChr is True:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % (chr_id, target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)
					else:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % (chr_id.replace('chr',''), target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)
				else:
					if withChr is True:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % ('chr' + chr_id, target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)
					else:
						print("##contig=<ID=%s,length=%d,assembly=%s>" % (chr_id, target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)

			print("##liftOverProgram=CrossMap(https://sourceforge.net/projects/crossmap/)", file=FILE_OUT)
			print("##liftOverFile=" + liftoverfile, file=FILE_OUT)
			print("##new_reference_genome=" + refgenome, file=FILE_OUT)
			print("##liftOverTime=" + datetime.date.today().strftime("%B%d,%Y"), file=FILE_OUT)
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
			printlog(["Lifting over ... "])

		else:
			if line.startswith('#'):continue

			# process non-variant region
			if 'END=' in line:
				fields = str.split(line,maxsplit=8)
				total_region += 1
				chrom = fields[0]
				start = int(fields[1])-1	 # 0 based
				try:
					m = re.search(r"END\=(\d+)", line)
					end = int(m[1])
				except:
					print (line + "\tFail(Unmap)", file=UNMAP)
					failed_region += 1
					continue

				a = map_coordinates(mapping, chrom, start, end,'+')
				if a is None:
					print (line + "\tFail(Unmap)", file=UNMAP)
					failed_region += 1
					continue
				if len(a) == 2:
					# update chrom
					target_chr = str(a[1][0])	#target_chr is from chain file, could be 'chr1' or '1'
					target_start = a[1][1]
					target_end = a[1][2]
					fields[0] = target_chr

					# update start coordinate
					fields[1] = target_start + 1

					# update END
					fields[7] = fields[7].replace(('END=' + str(end)), ('END=' + str(target_end)))
					print('\t'.join(map(str, fields)), file=FILE_OUT)

			# process variant line
			else:

				fields = str.split(line,maxsplit=7)
				total_var += 1
				chrom = fields[0]
				start = int(fields[1])-1	 	# 0 based, ref_allele start
				end = start + len(fields[3])	# ref_allele end
				alt_allele = fields[4].replace(' ','').split(',')[0]	# 20  10000598    .   T   A,<NON_REF> 1754.77 .   DP=54;

				a = map_coordinates(mapping, chrom, start, end,'+')
				if a is None:
					print (line + "\tFail(Unmap)", file=UNMAP)
					failed_var += 1
					continue

				if len(a) == 2:
					# update chrom
					target_chr = str(a[1][0])	#target_chr is from chain file, could be 'chr1' or '1'
					target_start = a[1][1]
					target_end = a[1][2]
					fields[0] = target_chr

					# update start coordinate
					fields[1] = target_start + 1

					# update ref allele
					target_chr = update_chromID(refFasta.references[0], target_chr)
					fields[3] = refFasta.fetch(target_chr,target_start,target_end).upper()

					if a[1][3] == '-':
						fields[4] = revcomp_DNA(alt_allele, True) + ',<NON_REF>'

					#ref_allele and alt_alele are different
					if fields[3] != alt_allele:
						print('\t'.join(map(str, fields)), file=FILE_OUT)
					else:
						print (line + "\tFail(REF==ALT)", file=UNMAP)
						failed_var += 1
				else:
					print (line + "\tFail(Multiple_hits)", file=UNMAP)
					failed_var += 1
					continue
	FILE_OUT.close()
	UNMAP.close()
	printlog (["Total variants:", str(total_var)])
	printlog (["Variants failed to map:", str(failed_var)])
	printlog (["Total non-variant regions:", str(total_region)])
	printlog (["Non-variant regions failed to map:", str(failed_region)])

def crossmap_maf_file(mapping, infile, outfile, liftoverfile, refgenome, ref_name):
	'''
	Convert genome coordinates in MAF (mutation annotation foramt) format.

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	infile : file
		Input file in VCF format. Can be a regular or compressed (*.gz, *.Z,*.z, *.bz,
		*.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to
		remote file.

	outfile : str
		prefix of output files.

	liftoverfile : file
		Chain (https://genome.ucsc.edu/goldenPath/help/chain.html) format file. Can be a
		regular or compressed (*.gz, *.Z,*.z, *.bz, *.bz2, *.bzip2) file, local file or
		URL (http://, https://, ftp://) pointing to remote file.

	refgenome : file
		The genome sequence file of 'target' assembly in FASTA format.
	ref_name : str
		The NCBI build name of the target assembly, for example, "GRCh37", "GRCh38".
	'''

	#index refegenome file if it hasn't been done
	if not os.path.exists(refgenome + '.fai'):
		printlog(["Creating index for", refgenome])
		pysam.faidx(refgenome)

	refFasta = pysam.Fastafile(refgenome)

	FILE_OUT = open(outfile ,'w')
	UNMAP = open(outfile + '.unmap','w')

	total = 0
	fail = 0

	for line in ireader.reader(infile):
		if not line.strip():
			continue
		line = line.strip()

		#meta-information lines needed in both mapped and unmapped files
		if line.startswith('#'):
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
			continue
		elif line.startswith('Hugo_Symbol'):
			print("#liftOver: Program=%sv%s, Time=%s, ChainFile=%s, NewRefGenome=%s" % ("CrossMap", __version__, datetime.date.today().strftime("%B%d,%Y"),liftoverfile,refgenome ), file=FILE_OUT)
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
			printlog(["Lifting over ... "])
		else:

			fields = str.split(line,sep = '\t')
			total += 1

			fields[3] = ref_name
			chrom = fields[4]
			start = int(fields[5])-1	 # 0 based
			end = int(fields[6])
			#strand = fields[7]

			a = map_coordinates(mapping, chrom, start, end,'+')

			if a is None:
				print(line, file=UNMAP)
				fail += 1
				continue

			if len(a) == 2:
				target_chr = str(a[1][0])	#target_chr is from chain file, could be 'chr1' or '1'
				target_start = a[1][1]
				target_end = a[1][2]

				# update chrom
				fields[4] = target_chr

				# update start coordinate
				fields[5] = target_start + 1

				# update end
				fields[6] = target_end

				# update ref allele
				target_chr = update_chromID(refFasta.references[0], target_chr)
				fields[10] = refFasta.fetch(target_chr,target_start,target_end).upper()

				if a[1][3] == '-':
					fields[10] = revcomp_DNA(fields[10], True)
				print('\t'.join(map(str, fields)), file=FILE_OUT)

			else:
				print(line, file=UNMAP)
				fail += 1
				continue
	FILE_OUT.close()
	UNMAP.close()
	printlog (["Total entries:", str(total)])
	printlog (["Failed to map:", str(fail)])


def crossmap_bed_file(mapping, inbed, outfile=None, min_ratio = 0.95):
	'''
	Convert genome coordinates (in bed format) between assemblies.
	BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	inbed : file
		Input BED file.

	outfile : str, optional
		Prefix of output files.

	min_ratio : float, optional
		Minimum ratio of query bases that must remap

	'''

	# check if 'outfile' was set. If not set, print to screen, if set, print to file
	if outfile is not None:
		FILE_OUT = open(outfile,'w')
		UNMAP = open(outfile + '.unmap','w')
	else:
		pass

	for line in ireader.reader(inbed):
		if line.startswith(('#', 'track','browser')):continue
		if not line.strip():continue
		line=line.strip()
		fields=line.split()
		strand = '+'

		# filter out line less than 3 columns
		if len(fields)<3:
			print("Less than 3 fields. skip " + line, file=sys.stderr)
			if outfile:
				print(line + '\tInvalidBedFormat', file=UNMAP)
			continue
		try:
			int(fields[1])
		except:
			print("Start corrdinate is not an integer. skip " + line, file=sys.stderr)
			if outfile:
				print(line + '\tInvalidStartPosition', file=UNMAP)
			continue
		try:
			int(fields[2])
		except:
			print("End corrdinate is not an integer. skip " + line, file=sys.stderr)
			if outfile:
				print(line + '\tInvalidEndPosition', file=UNMAP)
			continue
		if int(fields[1]) > int(fields[2]):
			print("\"Start\" is larger than \"End\" corrdinate is not an integer. skip " + line, file=sys.stderr)
			if outfile:
				print(line + '\tStart>End', file=UNMAP)
			continue

		# deal with bed less than 12 columns
		if len(fields)<12:

			# try to reset strand
			try:
				for f in fields:
					if f in ['+','-']:
						strand = f
			except:
				pass

			chrom = fields[0]
			start = int(fields[1])
			end = int(fields[2])
			total_query_length = end - start	#used to calculate q_map_ratio

			a = map_coordinates(mapping, chrom, start, end, strand)

			try:
				if (a is None) or (len(a) % 2 != 0):
					if outfile is None:
						print(line + '\tUnmap')
					else:
						print(line + '\tUnmap', file=UNMAP)
					continue

				#when a == 2, there is one-to-one match (i.e. 100% match)
				if len(a) == 2:
					#reset fields to target assembly
					fields[0] = a[1][0]
					fields[1] = a[1][1]
					fields[2] = a[1][2]
					for i in range(0,len(fields)):	#update the strand information
						if fields[i] in ['+','-']:
							fields[i] = a[1][3]

					if outfile is None:
						print(line + '\t->\t' + '\t'.join([str(i) for i in fields]))
					else:
						print('\t'.join([str(i) for i in fields]), file=FILE_OUT)

				#when a is an even number but bigger than 2, each segment is 100% match,
				# but the whole region is not. In this case, check *min_ratio* of the query
				if len(a) > 2 :
					(q_mapped_size, t_chroms, t_start, t_end) = cal_mapped_size(a)
					map_ratio = q_mapped_size/total_query_length

					if map_ratio < 0 or map_ratio >1:
						if outfile is None:
							print(line + '\tFail')
						else:
							print(line + '\tFail', file=UNMAP)
						continue

					if map_ratio >= min_ratio and len(t_chroms) == 1:
						t_chrom = t_chroms.pop()
						if outfile is None:
							print(line + '\t' + '\t'.join([t_chrom, str(t_start), str(t_end)]) + "\tmap_ratio=%.2f" % map_ratio)
						else:
							print('\t'.join([t_chrom, str(t_start), str(t_end)]) + "\tmap_ratio=%.2f" % map_ratio, file=FILE_OUT)
					else:
						count=0
						for j in range(1,len(a),2):
							count += 1
							fields[0] = a[j][0]
							fields[1] = a[j][1]
							fields[2] = a[j][2]
							for i in range(0,len(fields)):	#update the strand information
								if fields[i] in ['+','-']:
									fields[i] = a[j][3]

							if outfile is None:
								print(line + '\t'+ '(split.' + str(count) + ':' + ':'.join([str(i) for i in a[j-1]]) + ')\t' + '\t'.join([str(i) for i in fields]))
							else:
								print('\t'.join([str(i) for i in fields]), file=FILE_OUT)

			except:
				if outfile is None:
					print(line + '\tFail')
				else:
					print(line + '\tFail', file=UNMAP)
				continue

		# deal with bed12 and bed12+8 (genePred format)
		if len(fields)==12 or len(fields)==20:
			strand = fields[5]
			if strand not in ['+','-']:
				raise Exception("Unknown strand: %s. Can only be '+' or '-'." % strand)
			fail_flag = False
			exons_old_pos = annoGene.getExonFromLine(line)	#[[chr,st,end],[chr,st,end],...]
			#print exons_old_pos
			exons_new_pos = []
			for e_chr, e_start, e_end in exons_old_pos:
				# a has two elements, first is query, 2nd is target. # [('chr1', 246974830, 246974833,'+'), ('chr1', 248908207, 248908210,'+')]
				a = map_coordinates(mapping, e_chr, e_start, e_end, strand)
				if a is None:
					fail_flag =True
					break

				if len(a) == 2:
					exons_new_pos.append(a[1])
				else:
					fail_flag =True
					break

			if not fail_flag:
				# check if all exons were mapped to the same chromosome and the same strand
				chr_id = set()
				exon_strand = set()

				for e_chr, e_start, e_end, e_strand in exons_new_pos:
					chr_id.add(e_chr)
					exon_strand.add(e_strand)
				if len(chr_id) != 1 or len(exon_strand) != 1:
					fail_flag = True

				if not fail_flag:
					# build new bed
					cds_start_offset = int(fields[6]) - int(fields[1])
					cds_end_offset = int(fields[2]) - int(fields[7])
					new_chrom = exons_new_pos[0][0]
					new_chrom_st = exons_new_pos[0][1]
					new_chrom_end = exons_new_pos[-1][2]
					new_name = fields[3]
					new_score = fields[4]
					new_strand = exons_new_pos[0][3]
					new_thickStart = new_chrom_st + cds_start_offset
					new_thickEnd = new_chrom_end - cds_end_offset
					new_ittemRgb = fields[8]
					new_blockCount = len(exons_new_pos)
					new_blockSizes = ','.join([str(o - n) for m,n,o,p in exons_new_pos])
					new_blockStarts = ','.join([str(n - new_chrom_st) for m,n,o,p in exons_new_pos])

					new_bedline = '\t'.join(str(i) for i in (new_chrom,new_chrom_st,new_chrom_end,new_name,new_score,new_strand,new_thickStart,new_thickEnd,new_ittemRgb,new_blockCount,new_blockSizes,new_blockStarts))
					if check_bed12(new_bedline) is False:
						fail_flag = True
					else:
						if outfile is None:
							print(line + '\t->\t' + new_bedline)
						else:
							print(new_bedline, file=FILE_OUT)

			if fail_flag:
				if outfile is None:
					print(line + '\tFail')
				else:
					print(line, file=UNMAP)

def crossmap_gff_file(mapping, ingff,outfile = None):
	'''
	Description
	-----------
	Convert genome coordinates (in GFF/GTF format) between assemblies.
	GFF (General Feature Format) lines have nine required fields that must be Tab-separated:

	1. seqname - The name of the sequence. Must be a chromosome or scaffold.
	2. source - The program that generated this feature.
	3. feature - The name of this type of feature. Some examples of standard feature types
	   are "CDS", "start_codon", "stop_codon", and "exon".
	4. start - The starting position of the feature in the sequence. The first base is numbered 1.
	5. end - The ending position of the feature (inclusive).
	6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1
	   for this annotation data set, the score value will determine the level of gray in
	   which this feature is displayed (higher numbers = darker gray). If there is no score
	   value, enter ".".
	7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
	8. frame - If the feature is a coding exon, frame should be a number between 0-2 that
	   represents the reading frame of the first base. If the feature is not a coding exon,
	   the value should be '.'.
	9. group - All lines with the same group are linked together into a single item.

	GFF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format3

	GTF (Gene Transfer Format) is a refinement to GFF that tightens the specification. The
	first eight GTF fields are the same as GFF. The group field has been expanded into a
	list of attributes. Each attribute consists of a type/value pair. Attributes must end
	in a semi-colon, and be separated from any following attribute by exactly one space.

	GTF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format4

	We do NOT check if features (exon, CDS, etc) originally belonging to the same gene  were
	converted into the same chromosome/strand.

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	ingff : file
		Input GFF/GTF file.

	outfile : str, optional
		Prefix of output files.
	'''

	if outfile is not None:
		FILE_OUT = open(outfile,'w')
		UNMAP = open(outfile + '.unmap', 'w')

	for line in ireader.reader(ingff):
		if line.startswith(('#','track','browser','visibility')):continue
		if not line.strip():continue

		line=line.strip()
		fields=line.split('\t')
		try:
			start = int(fields[3]) - 1	#0-based
			end =  int(fields[4])/1
			feature_size = end - start
		except:
			print('Cannot recognize \"start\" and \"end\" coordinates. Skip ' + line, file=sys.stderr)
			if outfile:
				print(line, file=UNMAP)
			continue
		if fields[6] not in ['+','-','.']:
			print('Cannot recognize \"strand\". Skip ' + line, file=sys.stderr)
			if outfile:
				print(line, file=UNMAP)
			continue

		strand = '-' if fields[6] == '-' else '+'

		chrom = fields[0]
		a = map_coordinates(mapping, chrom,start,end,strand)

		if a is None:
			if outfile is None:
				print(line + '\tfail (no match to target assembly)')
			else:
				print(line, file=UNMAP)
			continue
		if len(a) !=2:
			if outfile is None:
				print(line + '\tfail (multpile match to target assembly)')
			else:
				print(line, file=UNMAP)
		else:
			if (int(a[1][2]) - int(a[1][1])) != feature_size:	# check if it is exact match
				if outfile is None:
					print(line + '\tfail (not exact match)')
				else:
					print(line, file=UNMAP)
			fields[0] = a[1][0]			# chrom
			fields[3] = int(a[1][1]) + 1	 # start, 1-based
			fields[4] = int(a[1][2])
			fields[6] = a[1][3]

			if outfile is None:
				print(line + '\t->\t' + '\t'.join([str(i) for i in fields]))
			else:
				print('\t'.join([str(i) for i in fields]), file=FILE_OUT)

def crossmap_bam_file(mapping, chainfile, infile,  outfile_prefix, chrom_size, IS_size=200, IS_std=30.0, fold=3, addtag = True):
	'''

	Description
	-----------
	Convert genome coordinates (in BAM/SAM format) between assemblies.
	BAM/SAM format: http://samtools.sourceforge.net/
	chrom_size is target chromosome size

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	chainfile : file
		Input chain format file.

	infile : file
		Input BAM, SAM or CRAM foramt file.

	outfile_prefix : str
		Output prefix.

	chrom_size : dict
		Chromosome size of the *target* assembly, used to build bam header.

	IS_size : int
		Average insert size of pair-end sequencing.

	IS_std : float
		Stanadard deviation of insert size.

	fold : float
		A mapped pair is considered as \"proper pair\" if both ends mapped to
		different strand and the distance between them is less then fold * stdev
		from the mean.

	addtag : bool
		if addtag is set to True, will add tags to each alignmnet:
			Q = QC (QC failed)
			N = unmapped (originally unmapped or originally mapped but failed
			    to liftover to new assembly)
			M = multiple mapped (alignment can be liftover to multiple places)
			U = unique mapped (alignment can be liftover to only 1 place)

		tags for pair-end sequencing include:
			QF: QC failed
			NN: both read1 and read2 unmapped
			NU: read1 unmapped, read2 unique mapped
			NM: read1 unmapped, multiple mapped
			UN: read1 uniquely mapped, read2 unmap
			UU: both read1 and read2 uniquely mapped
			UM: read1 uniquely mapped, read2 multiple mapped
			MN: read1 multiple mapped, read2 unmapped
			MU: read1 multiple mapped, read2 unique mapped
			MM: both read1 and read2 multiple mapped

		tags for single-end sequencing include:
			QF: QC failed
			SN: unmaped
			SM: multiple mapped
			SU: uniquely mapped
	'''

	# determine the input file format (BAM, CRAM or SAM)
	file_type = ''
	if infile.lower().endswith('.bam'):
		file_type = 'BAM'
		comments=['ORIGINAL_BAM_FILE=' + infile]
		samfile = pysam.Samfile(infile,'rb')
		if len(samfile.header) == 0:
			print("BAM file has no header section. Exit!", file=sys.stderr)
			sys.exit(1)
	elif infile.lower().endswith('.cram'):
		file_type = 'CRAM'
		comments=['ORIGINAL_CRAM_FILE=' + infile]
		samfile = pysam.Samfile(infile,'rc')
		if len(samfile.header) == 0:
			print("CRAM file has no header section. Exit!", file=sys.stderr)
			sys.exit(1)
	elif infile.lower().endswith('.sam'):
		file_type = 'SAM'
		comments=['ORIGINAL_SAM_FILE=' + infile]
		samfile = pysam.Samfile(infile,'r')
		if len(samfile.header) == 0:
			print("SAM file has no header section. Exit!", file=sys.stderr)
			sys.exit(1)
	else:
		print("Unknown file type! Input file must have suffix '.bam','.cram', or '.sam'.", file=sys.stderr)
		sys.exit(1)
	comments.append('CHAIN_FILE=' + chainfile)

	sam_ori_header = samfile.header.to_dict()

	# chromosome ID style of the original BAM file
	chrom_style = sam_ori_header['SQ'][0]['SN']	# either 'chr1' or '1'

	# update chrom_size of target genome
	target_chrom_sizes = {}
	for n,l in chrom_size.items():
		target_chrom_sizes[update_chromID(chrom_style, n)] = l

	(new_header, name_to_id) = sam_header.bam_header_generator(orig_header = sam_ori_header, chrom_size = target_chrom_sizes, prog_name="CrossMap",prog_ver = __version__, format_ver=1.0,sort_type = 'coordinate',co=comments)

	# write to file
	if outfile_prefix is not None:
		if file_type == 'BAM':
			OUT_FILE = pysam.Samfile( outfile_prefix + '.bam', "wb", header = new_header )
			printlog (["Liftover BAM file:", infile, '==>', outfile_prefix + '.bam'])
		elif file_type == 'CRAM':
			OUT_FILE = pysam.Samfile( outfile_prefix + '.bam', "wb", header = new_header )
			printlog (["Liftover CRAM file:", infile, '==>', outfile_prefix + '.bam'])
		elif file_type == 'SAM':
			OUT_FILE = pysam.Samfile( outfile_prefix + '.sam', "wh", header = new_header )
			printlog (["Liftover SAM file:", infile, '==>',	 outfile_prefix + '.sam'])
		else:
			print("Unknown file type! Input file must have suffix '.bam','.cram', or '.sam'.", file=sys.stderr)
			sys.exit(1)
	# write to screen
	else:
		if file_type == 'BAM':
			OUT_FILE = pysam.Samfile( '-', "wb", header = new_header )
			printlog (["Liftover BAM file:", infile])
		elif file_type == 'CRAM':
			OUT_FILE = pysam.Samfile( '-', "wb", header = new_header )
			printlog (["Liftover CRAM file:", infile])
		elif file_type == 'SAM':
			OUT_FILE = pysam.Samfile( '-', "w", header = new_header )
			printlog (["Liftover SAM file:", infile])
		else:
			print("Unknown file type! Input file must have suffix '.bam','.cram', or '.sam'.", file=sys.stderr)
			sys.exit(1)
	QF = 0
	NN = 0
	NU = 0
	NM = 0
	UN = 0
	UU = 0
	UM = 0
	MN = 0
	MU = 0
	MM = 0
	SN = 0
	SM = 0
	SU = 0
	total_item = 0
	try:
		while(1):
			total_item += 1
			old_alignment = next(samfile)
			new_alignment = pysam.AlignedRead() # create AlignedRead object

			new_alignment.query_name = old_alignment.query_name               # 1st column. read name.
			new_alignment.query_sequence = old_alignment.query_sequence       # 10th column. read sequence. all bases.
			new_alignment.query_qualities = old_alignment.query_qualities	  # 11th column. read sequence quality. all bases.
			new_alignment.set_tags(old_alignment.get_tags() )                 # 12 - columns


			# by default pysam will change RG:Z to RG:A, which can cause downstream failures with GATK and freebayes
			# Thanks Wolfgang Resch <wresch@helix.nih.gov> identified this bug and provided solution.
			try:
				rg, rgt = old_alignment.get_tag("RG", with_value_type=True)
			except KeyError:
				pass
			else:
				new_alignment.set_tag("RG", str(rg), rgt)


			## Pair-end sequencing
			if old_alignment.is_paired:
				new_alignment.flag = 0x1 #pair-end in sequencing
				if old_alignment.is_read1:
					new_alignment.flag = new_alignment.flag | 0x40
				elif old_alignment.is_read2:
					new_alignment.flag = new_alignment.flag | 0x80

				if old_alignment.is_qcfail:
					new_alignment.flag = new_alignment.flag | 0x200
					new_alignment.reference_id = -1                       #3
					new_alignment.reference_start = 0                     #4
					new_alignment.mapping_quality = 255				      #5
					new_alignment.cigartuples = old_alignment.cigartuples #6
					new_alignment.next_reference_id = -1   #7
					new_alignment.next_reference_start = 0 #8
					new_alignment.template_length = 0      #9

					QF += 1
					if addtag: new_alignment.set_tag(tag="QF", value=0)
					OUT_FILE.write(new_alignment)
					continue
				#==================================
				# R1 originally unmapped
				#==================================
				elif old_alignment.is_unmapped:
					new_alignment.flag = new_alignment.flag | 0x4         #2
					new_alignment.reference_id = -1                       #3
					new_alignment.reference_start = 0                     #4
					new_alignment.mapping_quality = 255			          #5
					new_alignment.cigartuples = old_alignment.cigartuples #6

					# R1 & R2 originally unmapped
					if old_alignment.mate_is_unmapped:
						new_alignment.next_reference_id = -1   #7
						new_alignment.next_reference_start = 0 #8
						new_alignment.template_length = 0      #9

						NN += 1
						if addtag: new_alignment.set_tag(tag="NN", value=0)
						OUT_FILE.write(new_alignment)
						continue
					# R1 unmap, R2 is mapped
					else:
						try:
							read2_chr = samfile.get_reference_name(old_alignment.next_reference_id)
							read2_strand = '-' if old_alignment.mate_is_reverse else '+'
							read2_start = old_alignment.next_reference_start
							read2_end = read2_start + 1
							read2_maps = map_coordinates(mapping, read2_chr, read2_start, read2_end, read2_strand)
						except:
							read2_maps = None

						#------------------------------------
						# R1 unmapped, R2 failed to liftover
						#------------------------------------
						if read2_maps is None:
							new_alignment.next_reference_id = -1   #7
							new_alignment.next_reference_start = 0 #8
							new_alignment.template_length = 0      #9

							NN += 1
							if addtag: new_alignment.set_tag(tag="NN", value=0)
							OUT_FILE.write(new_alignment)
							continue

						#------------------------------------
						# R1 unmapped, R2 unique
						#------------------------------------
						elif len(read2_maps) == 2:
							# 2-9
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							new_alignment.reference_id = name_to_id[read2_maps[1][0]]	#recommend to set the RNAME of unmapped read to its mate's
							new_alignment.reference_start = read2_maps[1][1]				#recommend to set the POS of unmapped read to its mate's
							new_alignment.mapping_quality = old_alignment.mapping_quality
							new_alignment.cigartuples = old_alignment.cigartuples
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.next_reference_start = read2_maps[1][1]
							new_alignment.template_length = 0

							NU += 1
							if addtag: new_alignment.set_tag(tag="NU", value=0)
							OUT_FILE.write(new_alignment)
							continue

						#------------------------------------
						# R1 unmapped, R2 multiple
						#------------------------------------
						else:
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							# 2-9
							new_alignment.flag = new_alignment.flag | 0x100
							new_alignment.reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.reference_start = read2_maps[1][1]
							new_alignment.mapping_quality = old_alignment.mapping_quality
							new_alignment.cigartuples = old_alignment.cigartuples
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.next_reference_start =  read2_maps[1][1]
							new_alignment.template_length = 0

							NM += 1
							if addtag: new_alignment.set_tag(tag="NM", value=0)
							OUT_FILE.write(new_alignment)
							continue
				#==================================
				# R1 is originally mapped
				#==================================
				else:
					try:
						read1_chr = samfile.get_reference_name(old_alignment.reference_id)
						read1_strand = '-' if old_alignment.is_reverse else '+'
						read1_start = old_alignment.reference_start
						read1_end = old_alignment.reference_end
						read1_maps = map_coordinates(mapping, read1_chr, read1_start, read1_end, read1_strand)
					except:
						read1_maps = None

					if not old_alignment.mate_is_unmapped:
						try:
							read2_chr = samfile.get_reference_name(old_alignment.next_reference_id)
							read2_strand = '-' if old_alignment.mate_is_reverse else '+'
							read2_start = old_alignment.next_reference_start
							read2_end = read2_start + 1
							read2_maps = map_coordinates(mapping, read2_chr, read2_start, read2_end, read2_strand)
						except:
							read2_maps = None
					#------------------------------------
					# R1 failed to liftover
					#------------------------------------
					if read1_maps is None:
						# read2 is unmapped or failed to convertion
						if old_alignment.mate_is_unmapped or (read2_maps is None):
							# col2 - col9
							new_alignment.flag = new_alignment.flag | 0x4  #2
							new_alignment.reference_id = -1                #3
							new_alignment.reference_start = 0              #4
							new_alignment.mapping_quality = 255            #5
							new_alignment.cigartuples = old_alignment.cigartuples #6
							new_alignment.next_reference_id = -1	           #7
							new_alignment.next_reference_start = 0         #8
							new_alignment.template_length = 0              #9

							if addtag: new_alignment.set_tag(tag="NN", value=0)
							NN += 1
							OUT_FILE.write(new_alignment)
							continue

						# read2 is unique mapped
						elif len(read2_maps) == 2:
							# col2 - col9
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							new_alignment.reference_id = name_to_id[read2_maps[1][0]]	 #recommend to set the RNAME of unmapped read to its mate's
							new_alignment.reference_start = read2_maps[1][1]				#recommend to set the POS of unmapped read to its mate's
							new_alignment.mapping_quality = old_alignment.mapping_quality
							new_alignment.cigartuples = old_alignment.cigartuples
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.next_reference_start = read2_maps[1][1]	#start
							new_alignment.template_length = 0

							NU += 1
							if addtag: new_alignment.set_tag(tag="NU", value=0)
							OUT_FILE.write(new_alignment)
							continue

						# read2 is multiple mapped
						else:
							# col2 - col9
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							new_alignment.flag = new_alignment.flag | 0x100
							new_alignment.reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.reference_start = read2_maps[1][1]
							new_alignment.mapping_quality = 255			# mapq not available
							new_alignment.cigartuples = old_alignment.cigartuples
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.next_reference_start =  read2_maps[1][1] #start
							new_alignment.template_length = 0

							NM += 1
							if addtag:new_alignment.set_tag(tag="NM", value=0)
							OUT_FILE.write(new_alignment)
							continue

					#------------------------------------
					# R1 uniquely mapped
					#------------------------------------
					elif len(read1_maps) == 2:
						# col2 - col5
						if read1_maps[1][3] == '-':new_alignment.flag = new_alignment.flag | 0x10
						new_alignment.reference_id = name_to_id[read1_maps[1][0]]
						new_alignment.reference_start = read1_maps[1][1]
						new_alignment.mapping_quality = old_alignment.mapping_quality

						if read1_maps[0][3] != read1_maps[1][3]:	# opposite strand
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples[::-1]			#reverse cigar tuple
							# 10
							new_alignment.query_sequence = revcomp_DNA(old_alignment.query_sequence)		#reverse complement read sequence
							# 11
							new_alignment.query_qualities = old_alignment.query_qualities[::-1]			#reverse quality string
						elif read1_maps[0][3] == read1_maps[1][3]:	#  same strand
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples



						# R2 unmapped before or after conversion
						if (old_alignment.mate_is_unmapped) or (read2_maps is None):
							#2,7-9
							new_alignment.flag = new_alignment.flag | 0x8
							new_alignment.next_reference_id = name_to_id[read1_maps[1][0]]
							new_alignment.next_reference_start =  read1_maps[1][1]
							new_alignment.template_length = 0

							UN += 1
							if addtag: new_alignment.set_tag(tag="UN", value=0)
							OUT_FILE.write(new_alignment)
							continue

						# R2 is unique mapped
						elif len(read2_maps)==2:
							# 2,7-9
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]	#chrom
							new_alignment.next_reference_start =  read2_maps[1][1]
							new_alignment.template_length = abs(new_alignment.reference_start - new_alignment.next_reference_start) + old_alignment.reference_length
							# 2
							if (read2_maps[1][3] != read1_maps[1][3]) and (new_alignment.template_length <= IS_size + fold * IS_std) and (new_alignment.template_length >= IS_size - fold * IS_std):
								new_alignment.flag = new_alignment.flag | 0x2

							UU += 1
							if addtag: new_alignment.set_tag(tag="UU", value=0)
							OUT_FILE.write(new_alignment)
							continue

						# R2 is multiple mapped
						else:
							# 2 (strand)
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							# 2 (secondary alignment)
							new_alignment.flag = new_alignment.flag | 0x100

							#7-9
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]
							new_alignment.next_reference_start =  read2_maps[1][1]
							new_alignment.template_length = 0

							UM += 1
							if addtag: new_alignment.set_tag(tag="UM", value=0)
							OUT_FILE.write(new_alignment)
							continue
					#------------------------------------
					# R1 multiple mapped
					#-----------------------------------
					elif len(read1_maps) > 2 and len(read1_maps) % 2 ==0:
						# 2
						new_alignment.flag = new_alignment.flag | 0x100
						if read1_maps[1][3] == '-':
							new_alignment.flag = new_alignment.flag | 0x10
						# 3-5
						new_alignment.tid = name_to_id[read1_maps[1][0]]	#chrom
						new_alignment.pos = read1_maps[1][1]	#start
						new_alignment.mapq = 255

						if read1_maps[0][3] != read1_maps[1][3]:	# opposite strand
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples[::-1]			        #reverse cigar tuple
							# 10
							new_alignment.query_sequence = revcomp_DNA(old_alignment.query_sequence)		#reverse complement read sequence
							# 11
							new_alignment.query_qualities = old_alignment.query_qualities[::-1]			#reverse quality string
						elif read1_maps[0][3] == read1_maps[1][3]:	#  same strand
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples

						# (1) R2 is unmapped
						if (old_alignment.mate_is_unmapped) or (read2_maps is None):
							#2,7-9
							new_alignment.flag = new_alignment.flag | 0x8
							new_alignment.next_reference_id = name_to_id[read1_maps[1][0]]
							new_alignment.next_reference_start =  read1_maps[1][1]
							new_alignment.template_length = 0

							MN += 1
							if addtag: new_alignment.set_tag(tag="MN", value=0)
							OUT_FILE.write(new_alignment)
							continue

						# (2) read2 is unique mapped
						elif len(read2_maps)==2:
							# 2,7-9
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]	#chrom
							new_alignment.next_reference_start =  read2_maps[1][1]
							new_alignment.template_length = 0

							MU += 1
							if addtag: new_alignment.set_tag(tag="MU", value=0)
							OUT_FILE.write(new_alignment)
							continue

						# (3) R2 is multiple mapped
						else:
							# 2,7-9
							if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
							# 2 (secondary alignment)
							new_alignment.flag = new_alignment.flag | 0x100
							new_alignment.next_reference_id = name_to_id[read2_maps[1][0]]	#chrom
							new_alignment.next_reference_start =  read2_maps[1][1]
							new_alignment.template_length = 0

							MM += 1
							if addtag: new_alignment.set_tag(tag="MM", value=0)
							OUT_FILE.write(new_alignment)
							continue

			# Singel end sequencing
			else:
				# 7-9
				new_alignment.next_reference_id = -1
				new_alignment.next_reference_start = 0
				new_alignment.template_length = 0

				# (1) originally unmapped
				if old_alignment.is_unmapped:
					# 2-6
					new_alignment.flag = new_alignment.flag | 0x4
					new_alignment.reference_id = -1
					new_alignment.reference_start = 0
					new_alignment.mapping_quality = 255
					new_alignment.cigartuples = old_alignment.cigartuples

					SN += 1
					if addtag: new_alignment.set_tag(tag="SN",value=0)
					OUT_FILE.write(new_alignment)
					continue
				else:
					new_alignment.flag = 0x0
					read_chr = samfile.get_reference_name(old_alignment.reference_id)
					read_strand = '-' if old_alignment.is_reverse else '+'
					read_start = old_alignment.reference_start
					read_end = old_alignment.reference_end
					read_maps = map_coordinates(mapping, read_chr, read_start, read_end, read_strand)

					# (2) unmapped afte liftover
					if read_maps is None:
						new_alignment.flag = new_alignment.flag | 0x4
						new_alignment.reference_id = -1
						new_alignment.reference_start = 0
						new_alignment.mapping_quality = 255

						SN += 1
						if addtag: new_alignment.set_tag(tag="SN",value=0)
						OUT_FILE.write(new_alignment)
						continue

					# (3) unique mapped
					if len(read_maps)==2:
						if read_maps[1][3] == '-':
							new_alignment.flag = new_alignment.flag | 0x10
						if read_maps[0][3] != read_maps[1][3]:
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples[::-1]			#reverse cigar tuple
							# 10
							new_alignment.query_sequence = revcomp_DNA(old_alignment.query_sequence)		#reverse complement read sequence
							# 11
							try:
								new_alignment.query_qualities = old_alignment.query_qualities[::-1]			#reverse quality string
							except:
								new_alignment.query_qualities = []
						else:
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples

						# 3-5
						new_alignment.reference_id = name_to_id[read_maps[1][0]]
						new_alignment.reference_start = read_maps[1][1]
						new_alignment.mapping_quality = old_alignment.mapping_quality

						SU += 1
						if addtag: new_alignment.set_tag(tag="SU",value=0)
						OUT_FILE.write(new_alignment)
						continue

					# (4) multiple mapped
					if len(read_maps) > 2 and len(read_maps) % 2 ==0:
						new_alignment.flag = new_alignment.flag | 0x100
						if read_maps[1][3] == '-':
							new_alignment.flag = new_alignment.flag | 0x10
						if read_maps[0][3] != read_maps[1][3]:
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples[::-1]			#reverse cigar tuple
							# 10
							new_alignment.query_sequence = revcomp_DNA(old_alignment.query_sequence)		#reverse complement read sequence
							# 11
							new_alignment.query_qualities = old_alignment.query_qualities[::-1]			#reverse quality string
						else:
							# 6
							new_alignment.cigartuples = old_alignment.cigartuples

						# 3-5
						new_alignment.tid = name_to_id[read_maps[1][0]]
						new_alignment.pos = read_maps[1][1]
						new_alignment.mapq = old_alignment.mapq

						SM += 1
						if addtag:	new_alignment.set_tag(tag="SM",value=0)
						OUT_FILE.write(new_alignment)
						continue
	except StopIteration:
		printlog(["Done!"])
	OUT_FILE.close()

	if outfile_prefix is not None:
		if file_type == "BAM" or file_type == "CRAM":
			try:
				printlog (['Sort "%s" and save as "%s"' % (outfile_prefix + '.bam', outfile_prefix + '.sorted.bam')])
				pysam.sort("-o",  outfile_prefix + '.sorted.bam', outfile_prefix + '.bam')
			except:
				printlog(["Warning: ","output BAM file was NOT sorted"])
			try:
				printlog (['Index "%s" ...' % (outfile_prefix + '.sorted.bam')])
				pysam.index(outfile_prefix + '.sorted.bam',outfile_prefix + '.sorted.bam.bai')
			except:
				printlog(["Warning: ","output BAM file was NOT indexed."])

	print("Total alignments:" + str(total_item-1))
	print("	 QC failed: " + str(QF))
	if max(NN,NU, NM, UN, UU, UM, MN, MU, MM) > 0:
		print("	 Paired-end reads:")
		print("\tR1 unique, R2 unique (UU): " + str(UU))
		print("\tR1 unique, R2 unmapp (UN): " + str(UN))
		print("\tR1 unique, R2 multiple (UM): " + str(UM))

		print("\tR1 multiple, R2 multiple (MM): " + str(MM))
		print("\tR1 multiple, R2 unique (MU): " + str(MU))
		print("\tR1 multiple, R2 unmapped (MN): " + str(MN))

		print("\tR1 unmap, R2 unmap (NN): " + str(NN))
		print("\tR1 unmap, R2 unique (NU): " + str(NU))
		print("\tR1 unmap, R2 multiple (NM): " + str(NM))
	if max(SN,SU,SM) > 0:
		print("	 Single-end reads:")
		print("\tUniquley mapped (SU): " +	str(SU))
		print("\tMultiple mapped (SM): " +	str(SM))
		print("\tUnmapped (SN): " + str(SN))


def crossmap_wig_file(mapping, in_file, out_prefix, taget_chrom_size, in_format, binSize=100000):
	'''
	Description
	-----------
	Convert genome coordinates (in wiggle/bigwig format) between assemblies.
	wiggle format: http://genome.ucsc.edu/goldenPath/help/wiggle.html
	bigwig format: http://genome.ucsc.edu/goldenPath/help/bigWig.html

	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.

	in_file : file
		Input file in wig or bigwig format. Both "variableStep" and "fixedStep" wiggle
		lines are supported.

	out_prefix : str
		Prefix of output files.

	taget_chrom_size : dict
		Chromosome size of target genome. Key is chromsome ID, value is the length of the
		chromosome.

	in_format : str
		Either "wiggle" or "bigwig"

	binSize : int
		The chunk size when reading bigwig file in each iteration.
	'''

	OUT_FILE1 = open(out_prefix + '.bgr','w')	# original bgr file
	OUT_FILE2 = open(out_prefix + '.sorted.bgr','w')	# sorted bgr file
	OUT_FILE3 = pyBigWig.open(out_prefix + '.bw', "w")	# bigwig file

	# bigwig header
	target_chroms_sorted = []
	for k in sorted(taget_chrom_size.keys()):
		i_chrom = update_chromID('chr1',k)
		i_value = taget_chrom_size[k]
		target_chroms_sorted.append((i_chrom, i_value))

	if in_format.upper() == "WIGGLE":
		printlog (["Liftover wiggle file:", in_file, '==>', out_prefix + '.bgr'])

		for chrom, start, end, strand, score in wiggleReader (in_file):
			maps = map_coordinates(mapping, chrom, start, end, '+')
			if maps is None:
				continue
			if len(maps) == 2:
				print('\t'.join([str(i) for i in [maps[1][0],maps[1][1],maps[1][2], score]]), file=OUT_FILE1)
			else:
				continue
			maps[:]=[]
		OUT_FILE1.close()

		printlog (["Merging overlapped entries in bedGraph file ..."])
		for (chrom, start, end, score) in bgrMerge.merge(out_prefix + '.bgr'):
			print('\t'.join([str(i) for i in (chrom, start, end, score )]), file=OUT_FILE2)
		OUT_FILE2.close()

		os.remove(out_prefix + '.bgr')	#remove .bgr, keep .sorted.bgr

		# add header to bigwig file
		printlog (["Writing header to \"%s\" ..." % (out_prefix + '.bw') ])
		OUT_FILE3.addHeader(target_chroms_sorted)

		# add entries to bigwig file
		printlog (["Writing entries to \"%s\" ..." % (out_prefix + '.bw') ])
		for line in ireader.reader(out_prefix + '.sorted.bgr'):
			r_chr,r_st,r_end,r_value = line.split()
			OUT_FILE3.addEntries([r_chr], [int(r_st)], ends=[int(r_end)], values=[float(r_value)])

		OUT_FILE3.close()

	elif in_format.upper() == "BIGWIG":

		printlog (["Liftover bigwig file:", in_file, '==>', out_prefix + '.bgr'])
		for chrom, start, end, score in bigwigReader(in_file):
			maps = map_coordinates(mapping, chrom, start, end, '+')
			try:
				if maps is None: continue
				if len(maps) == 2:
					print('\t'.join([str(i) for i in [maps[1][0],maps[1][1],maps[1][2], score]]), file=OUT_FILE1)
				else:
					continue
			except:
				continue
			maps[:]=[]
		OUT_FILE1.close()

		printlog (["Merging overlapped entries in bedGraph file ..."])
		for (chrom, start, end, score) in bgrMerge.merge(out_prefix + '.bgr'):
			print('\t'.join([str(i) for i in (chrom, start, end, score )]), file=OUT_FILE2)
		OUT_FILE2.close()
		os.remove(out_prefix + '.bgr')	#remove .bgr, keep .sorted.bgr

		# add header to bigwig file
		printlog (["Writing header to \"%s\" ..." % (out_prefix + '.bw') ])
		OUT_FILE3.addHeader(target_chroms_sorted)

		# add entries to bigwig file
		printlog (["Writing entries to \"%s\" ..." % (out_prefix + '.bw') ])
		for line in ireader.reader(out_prefix + '.sorted.bgr'):
			r_chr,r_st,r_end,r_value = line.split()
			OUT_FILE3.addEntries([r_chr], [int(r_st)], [int(r_end)], [float(r_value)])
		OUT_FILE3.close()
	else:
		raise Exception("Unknown foramt. Must be 'wiggle' or 'bigwig'")


def general_help(cmds):
	desc="""CrossMap is a program to convert genome coordinates between different reference assemblies (e.g. from human hg18 to hg19 or vice versa). \
The supported file formats include BAM, BED, BigWig, CRAM, GFF, GTF, GVCF, MAF (mutation annotation format), SAM, Wiggle and VCF."""

	print("Program: %s (v%s)" % ("CrossMap", __version__), file=sys.stderr)
	print("\nDescription: \n%s" % '\n'.join('  '+i for i in wrap(desc,width=80)), file=sys.stderr)
	print("\nUsage: CrossMap.py <command> [options]\n", file=sys.stderr)
	for k in sorted(cmds):
		print('	 ' + k + '\t' + cmds[k], file=sys.stderr)
	print(file=sys.stderr)

def gff_help():
	msg =[
	('Usage', "CrossMap.py gff chain_file input.gff output_file"),
	('Description', "Convert GFF or GTF format file. The\"chain_file\" can be regular or compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file. Input file must be in GFF or GTF format. GFF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format3 GTF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format4"),
	('Example1 (write output to file)', "CrossMap.py gff  hg19ToHg18.over.chain.gz test.hg19.gtf test.hg18.gtf"),
	('Example2 (write output to screen)', "CrossMap.py gff	hg19ToHg18.over.chain.gz test.hg19.gtf"),
	]
	for i,j in msg:
		print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

def wig_help():
	msg =[
	('Usage', "CrossMap.py wig chain_file input.wig output_prefix"),
	('Description', "Convert wiggle format file. The \"chain_file\" can be regular or compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file.  Both \"variableStep\" and \"fixedStep\" wiggle lines are supported. Wiggle format: http://genome.ucsc.edu/goldenPath/help/wiggle.html"),
	('Example', "CrossMap.py wig hg18ToHg19.over.chain.gz test.hg18.wig test.hg19"),
	]
	for i,j in msg:
		print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

def bigwig_help():
	msg =[
	('Usage', "CrossMap.py bigwig chain_file input.bigwig output_prefix"),
	('Description', "Convert BigWig format file. The \"chain_file\" can be regular or compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file. Bigwig format: http://genome.ucsc.edu/goldenPath/help/bigWig.html"),
	('Example', "CrossMap.py bigwig hg18ToHg19.over.chain.gz test.hg18.bw test.hg19"),
	]
	for i,j in msg:
		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)


def vcf_help():
	msg =[
	("usage","CrossMap.py vcf chain_file input.vcf refGenome.fa output_file"),
	("Description", "Convert VCF format file. The \"chain_file\" and \"input.vcf\" can be regular or compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file. \"refGenome.fa\" is genome sequence file of the *target assembly*."),
	("Example", " CrossMap.py vcf hg19ToHg18.over.chain.gz test.hg19.vcf hg18.fa test.hg18.vcf"),
	]
	for i,j in msg:
		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

def gvcf_help():
	msg =[
	("usage","CrossMap.py gvcf chain_file input.gvcf refGenome.fa output_file"),
	("Description", "Convert GVCF format file. The \"chain_file\" and \"input.gvcf\" can be regular or compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file. \"refGenome.fa\" is genome sequence file of the *target assembly*."),
	("Example", " CrossMap.py gvcf hg19ToHg18.over.chain.gz test.hg19.gvcf hg18.fa test.hg18.gvcf"),
	]
	for i,j in msg:
		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)


def maf_help():
	msg =[
	("usage","CrossMap.py maf chain_file input.maf refGenome.fa build_name output_file"),
	("Description", "Convert MAF format file. The \"chain_file\" and \"input.maf\" can be regular or compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file. \"refGenome.fa\" is genome sequence file of *target assembly*. \"build_name\" is the name of the *target_assembly* (eg \"GRCh38\")"),
	("Example", " CrossMap.py  maf	hg19ToHg38.over.chain.gz  test.hg19.maf	 hg38.fa  GRCh38 test.hg38.maf"),
	]
	for i,j in msg:
		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)


if __name__=='__main__':
	commands = {
	'bed':'convert BED, bedGraph or other BED-like file.',
	'bam':'convert BAM, CRAM or SAM format file.',
	'gff':'convert GFF or GTF format file.',
	'wig':'convert Wiggle or bedGraph format file.',
	'bigwig':'convert BigWig file.',
	'vcf':'convert VCF file.',
	'gvcf':'convert GVCF file.',
	'maf':'convert MAF (mutation annotation format) file.'
	}
	kwds = list(commands.keys())

	if len(sys.argv) == 1:
		general_help(commands)
		sys.exit(0)
	elif len(sys.argv) >=2:
		# deal with bed input
		if sys.argv[1].lower() == 'bed':
			usage="\nCrossMap.py bed chain_file bed_file [output_file] [options]\n\nExamples:\nCrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed test.hg19.bed    # write to file\nCrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed    # write to screen"
			parser = optparse.OptionParser(usage, add_help_option=False)
			parser.add_option("-r", "--ratio", action="store",type="float",dest="min_map_ratio", default=1.0, help="Minimum ratio of bases that must remap. [default=%default]" )
			(options,args)=parser.parse_args()

			if len(args) == 3:
				chain_file = args[1]
				in_file = args[2]
				out_file = None
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file, print_table = False)
				crossmap_bed_file(mapTree, in_file, out_file, min_ratio = options.min_map_ratio)
			elif len(args) == 4:
				chain_file = args[1]
				in_file = args[2]
				out_file = args[3]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_bed_file(mapTree, in_file, out_file, min_ratio = options.min_map_ratio)
			else:
				parser.print_help()
				sys.exit(0)

		elif sys.argv[1].lower() == 'gff':
			if len(sys.argv) == 4:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_gff_file(mapTree, in_file, None)
			elif len(sys.argv) == 5:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				out_file = sys.argv[4]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_gff_file(mapTree, in_file, out_file)
			else:
				gff_help()
				sys.exit(0)
		elif sys.argv[1].lower() == 'wig':
			if len(sys.argv) == 5:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				out_file = sys.argv[4]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_wig_file(mapTree, in_file, out_file, targetChromSizes, in_format = 'wiggle')
			else:
				wig_help()
				sys.exit(0)
		elif sys.argv[1].lower() == 'bigwig':
			if len(sys.argv) == 5:
				chain_file = sys.argv[2]

				in_file = sys.argv[3]
				try:
					bw = pyBigWig.open(in_file)
				except:
					print ("\nPlease check if \"%s\" is in bigWig format!\n" % in_file, file=sys.stderr)
					sys.exit(0)

				out_file = sys.argv[4]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_wig_file(mapTree, in_file, out_file, targetChromSizes, in_format = 'bigwig')
			else:
				bigwig_help()
				sys.exit(0)
		elif sys.argv[1].lower() == 'bam':
			usage="\nCrossMap.py bam chain_file input_file [output_file] [options]\n\nNote:\nIf output_file is \"STDOUT\", \"-\" or missing, CrossMap will write BAM file to the screen"
			parser = optparse.OptionParser(usage, add_help_option=False)
			parser.add_option("-m", "--mean", action="store",type="float",dest="insert_size", default=200.0, help="Average insert size of pair-end sequencing (bp). [default=%default]")
			parser.add_option("-s", "--stdev", action="store",type="float",dest="insert_size_stdev", default=30.0, help="Stanadard deviation of insert size. [default=%default]" )
			parser.add_option("-t", "--times", action="store",type="float",dest="insert_size_fold", default=3.0, help="A mapped pair is considered as \"proper pair\" if both ends mapped to different strand and the distance between them is less then '-t' * stdev from the mean. [default=%default]")
			parser.add_option("-a","--append-tags",action="store_true",dest="add_tags",help="Add tag to each alignment.")
			(options,args)=parser.parse_args()

			if len(args) >= 3:
				print("Insert size = %f" % (options.insert_size), file=sys.stderr)
				print("Insert size stdev = %f" % (options.insert_size_stdev), file=sys.stderr)
				print("Number of stdev from the mean = %f" % (options.insert_size_fold), file=sys.stderr)
				if options.add_tags:
					print("Add tags to each alignment = %s" % ( options.add_tags), file=sys.stderr)
				else:
					print("Add tags to each alignment = %s" % ( False), file=sys.stderr)
				chain_file = args[1]
				in_file = args[2]
				out_file = args[3] if len(args) >= 4 else None
				(mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)

				if out_file in ["STDOUT","-"]:
					out_file = None
				if options.add_tags:
					crossmap_bam_file(mapping = mapTree, chainfile = chain_file, infile = in_file, outfile_prefix = out_file, chrom_size = targetChromSizes, IS_size=options.insert_size, IS_std=options.insert_size_stdev, fold=options.insert_size_fold,addtag=True)
				else:
					crossmap_bam_file(mapping = mapTree, chainfile = chain_file, infile = in_file, outfile_prefix = out_file, chrom_size = targetChromSizes, IS_size=options.insert_size, IS_std=options.insert_size_stdev, fold=options.insert_size_fold,addtag=False)
			else:
				parser.print_help()

		elif sys.argv[1].lower() == 'vcf':
			if len(sys.argv) == 6:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				genome_file = sys.argv[4]
				out_file = sys.argv[5]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_vcf_file(mapping = mapTree, infile= in_file, outfile = out_file, liftoverfile = sys.argv[2], refgenome = genome_file)
			else:
				vcf_help()
				sys.exit(0)

		elif sys.argv[1].lower() == 'gvcf':
			if len(sys.argv) == 6:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				genome_file = sys.argv[4]
				out_file = sys.argv[5]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_gvcf_file(mapping = mapTree, infile= in_file, outfile = out_file, liftoverfile = sys.argv[2], refgenome = genome_file)
			else:
				gvcf_help()
				sys.exit(0)

		elif sys.argv[1].lower() == 'maf':	#mapping, infile, outfile, liftoverfile, refgenome, ref_name
			if len(sys.argv) == 7:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				genome_file = sys.argv[4]
				build_name = sys.argv[5]
				out_file = sys.argv[6]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_maf_file(mapping = mapTree, infile= in_file, outfile = out_file, liftoverfile = chain_file, refgenome = genome_file, ref_name = build_name)
			else:
				maf_help()
				sys.exit(0)

		else:
			general_help()
			sys.exit(0)

