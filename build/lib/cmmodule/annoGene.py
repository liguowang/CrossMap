import collections
from bx.intervals import *
from cmmodule import BED

'''Compare given bed entry to reference gene model'''

def getCDSExonFromFile(bedfile):
	
	'''Only Extract CDS exon regions from input bed file (must be 12-column).'''		
	ret_lst=[]
	for f in open(bedfile,'r'):
		f = f.strip().split()
		chrom = f[0]
		chrom_start = int(f[1])
		name = f[4]
		strand = f[5]
		cdsStart = int(f[6])
		cdsEnd = int(f[7])
		blockCount = int(f[9])
		blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
		blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
          # grab cdsStart - cdsEnd
		cds_exons = []
		cds_seq = ''
		genome_seq_index = []
		
		chrom = chrom + ':' + strand
		for base,offset in zip( blockStarts, blockSizes ):
			if (base + offset) < cdsStart: continue
			if base > cdsEnd: continue
			exon_start = max( base, cdsStart )
			exon_end = min( base+offset, cdsEnd ) 
			#cds_exons.append( (exon_start, exon_end) )
			ret_lst.append([chrom,exon_start,exon_end])
	return ret_lst

def getUTRExonFromFile(bedfile,utr=35):
	'''Only Extract UTR regions from input bed file (must be 12-column). output is 6-column bed format.
	When utr=35 [default], extract both 5' and 3' UTR. When utr=3, only extract 3' UTR. When utr=5,
	only extract 5' UTR'''
	
	ret_lst=[]
	for line in open(bedfile,'r'):
		if line.startswith('#'):continue
		if line.startswith('track'):continue
		if line.startswith('browser'):continue
		fields=line.rstrip('\r\n').split()
		chrom=fields[0]
		strand=fields[5]
		txStart=int(fields[1])
		txEnd=int(fields[2])
		cdsStart=int(fields[6])
		cdsEnd=int(fields[7])		
		exon_start=list(map(int,fields[11].rstrip(',').split(',')))
		exon_start=list(map((lambda x: x + txStart),exon_start))
			
		exon_end=list(map(int,fields[10].rstrip(',').split(',')))
		exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
			   
		chrom = chrom + ':' + strand
		if (utr==35 or utr==5):
			for st,end in zip(exon_start,exon_end):
				if st < cdsStart:
					utr_st = st
					utr_end = min(end,cdsStart)
					ret_lst.append([chrom,utr_st,utr_end])					
		if (utr==35 or utr==3):
			for st,end in zip(exon_start,exon_end):
				if end > cdsEnd:
					utr_st = max(st, cdsEnd)
					utr_end = end
					ret_lst.append([chrom,utr_st,utr_end]) 
	return ret_lst



def getExonFromFile(bedfile):
	'''Extract ALL exon regions from input bed file (must be 12-column). return list of [chrom:+ st end]'''
	
	ret_lst=[]
	for line in open(bedfile,'r'):
		try:
			if line.startswith('#'):continue
			if line.startswith('track'):continue
			if line.startswith('browser'):continue
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
		chrom = chrom + ':' + strand
		for st,end in zip(exon_start,exon_end):
			ret_lst.append([chrom,st,end])
	return ret_lst

def getExonFromFile2(bedfile):
	'''Extract ALL exon regions from input bed file (must be 12-column). return dict'''
	
	ret_dict_full = collections.defaultdict(set)	
	#ret_dict_inner = collections.defaultdict(set)	#trim off start_of_1st_exon and end_of_last_exon
	for line in open(bedfile,'r'):
		tmp=[]
		try:
			if line.startswith('#'):continue
			if line.startswith('track'):continue
			if line.startswith('browser'):continue
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
			key = chrom + ":" + txstart + "-" + txEnd + ":" + strand + ':' + geneName
		except:
			print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
			continue
		for st,end in zip(exon_start,exon_end):
			tmp.append(exon_start,exon_end)
		ret_dict_full[key] = set(tmp)
		#ret_dict_inner[key] = set(tmp[1:-1])
	return ret_dict_full



def getUTRExonFromLine(bedline,utr=35):
	'''Extract UTR regions from input bed line. When utr=35 [default], extract both
	5' and 3' UTR. When utr=3, only extract 3' UTR. When utr=5,only extract 5' UTR'''
	
	ret_lst=[]
	line = bedline
	if line.startswith('#'):return None
	if line.startswith('track'):return None
	if line.startswith('browser'):return None
	if not line.strip():return None
	fields=line.rstrip('\r\n').split()
	chrom=fields[0]
	strand=fields[5]
	txStart=int(fields[1])
	txEnd=int(fields[2])
	cdsStart=int(fields[6])
	cdsEnd=int(fields[7])		
	exon_start=list(map(int,fields[11].rstrip(',').split(',')))
	exon_start=list(map((lambda x: x + txStart),exon_start))
		
	exon_end=list(map(int,fields[10].rstrip(',').split(',')))
	exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
	
	chrom = chromm + ':' + strand
	if (utr==35 or utr==5):
		for st,end in zip(exon_start,exon_end):
			if st < cdsStart:
				utr_st = st
				utr_end = min(end,cdsStart)
				ret_lst.append([chrom,utr_st,utr_end])					
	if (utr==35 or utr==3):
		for st,end in zip(exon_start,exon_end):
			if end > cdsEnd:
				utr_st = max(st, cdsEnd)
				utr_end = end
				ret_lst.append([chrom,utr_st,utr_end])
	return ret_lst
		

def getCDSExonFromLine(bedline):
	
	'''Extract CDS exon regions from input bed line (must be 12-column).'''		
	ret_lst=[]
	line = bedline
	if line.startswith('#'):return None
	if line.startswith('track'):return None
	if line.startswith('browser'):return None
	if not line.strip():return None
	f = line.strip().split()
	chrom = f[0]
	chrom_start = int(f[1])
	name = f[4]
	strand = f[5]
	cdsStart = int(f[6])
	cdsEnd = int(f[7])
	blockCount = int(f[9])
	blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
	blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
        # grab cdsStart - cdsEnd
	cds_exons = []
	cds_seq = ''
	genome_seq_index = []
	chrom = chromm + ':' + strand
	for base,offset in zip( blockStarts, blockSizes ):
		if (base + offset) < cdsStart: continue
		if base > cdsEnd: continue
		exon_start = max( base, cdsStart )
		exon_end = min( base+offset, cdsEnd ) 
		#cds_exons.append( (exon_start, exon_end) )
		ret_lst.append([chrom,exon_start,exon_end])	
	return ret_lst

def getExonFromLine(bedline):
	'''Extract ALL exon regions from input bed line (must be 12-column). return list of [chrom st end]'''
	
	ret_lst=[]
	line = bedline
	#if line.startswith('#'):continue
	#if line.startswith('track'):continue
	#if line.startswith('browser'):continue
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
	#chrom = chrom + ':' + strand
	for st,end in zip(exon_start,exon_end):
		ret_lst.append([chrom, st, end])
	return ret_lst

	
def annotateBed(inputbed,refbed,outfile):
	'''compare inputbed to refbed'''
	ref_exon_ranges = {}
	ref_exon_starts = collections.defaultdict(set)		#key='chrom:+', value=set()
	ref_exon_ends = collections.defaultdict(set)
	OF = open(outfile,'w')

	#read reference bed file
	#Extract CDS exons from reference bed
	tmp = getCDSExonFromFile(refbed)
	for i in tmp:	#chr:+, st, end (0-base)
		ref_exon_starts[i[0]].add(int(i[1]))
		ref_exon_ends[i[0]].add(int(i[2]))
		if i[0] not in ref_exon_ranges:
			ref_exon_ranges[i[0]] = Intersecter()
		ref_exon_ranges[i[0]].add_interval( Interval( int(i[1]), int(i[2]) ))
		
	#Extract UTR exons from reference bed
	tmp = getUTRExonFromFile(refbed)
	for i in tmp:	#chr:+, st, end (0-base)
		ref_exon_starts[i[0]].add(int(i[1]))
		ref_exon_ends[i[0]].add(int(i[2]))
		if i[0] not in ref_exon_ranges:
			ref_exon_ranges[i[0]] = Intersecter()
		ref_exon_ranges[i[0]].add_interval( Interval( int(i[1]), int(i[2]) ))
	
	#prepare data structure 
	ref_exon_chain = getExonFromFile2(refbed)

	#read input bed
	for line in open(inputbed,'r'):
		if line.startswith('#'):continue
		if line.startswith('track'):continue
		if line.startswith('browser'):continue
		if not line.strip(): continue
		line = line.strip()
		fields=line.split()
		
		chrom = fields[0]
		strand = fields[5]
		tx_start = int(fields[1])
		tx_end = int(fields[2])
		key = chrom + ":" +strand
		if key in ref_exon_ranges:
			if  len(ref_exon_ranges[key].find(tx_start,tx_end))==0:	#input gene does NOT overlap with any known exons
				print(line + '\t' + 'novel(intergenic)')
			else:
				input_exon_chain=getExonFromLine(line)
				#print line + '\t' + 'overlap'
				#utr_3_exons = getUTRExon(line,utr=3)
				#utr_5_exons = getUTRExon(line,utr=5)
				#cds_exons = getCDSExon(line)
		else:
			print(line + '\t' + 'unknownChrom')
		
		
		#for utr3 in utr_3_exons:
		#	(chrom, st, end) = (utr3[0], int(utr3[1]),int(utr3[2]))
		#	if chrom in ref_exon_ranges:		
		#		if len(ref_exon_ranges[chrom].find(st,end))>0 :		#input exon overlap with known exon		
		#	else:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
