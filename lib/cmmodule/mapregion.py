import sys
from cmmodule  import ireader
from cmmodule.utils import map_coordinates

def crossmap_region_file(mapping, inbed, outfile=None, min_ratio = 0.85):
	'''
	Convert large genomic regions (in bed format) between assemblies.
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
		# input: 'chr1',246974830,247024835
		# output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
		# [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]



		if (a is None) or (len(a) % 2 != 0):
			if outfile is None:
				print(line + '\tFail\tUnmap')
			else:
				print(line + '\tFail\tUnmap', file=UNMAP)
			continue

		#when a == 2, there is one-to-one match (i.e. 100% match)
		if len(a) == 2:
			#reset fields to target assembly
			fields[0] =  a[1][0]
			fields[1] =  a[1][1]
			fields[2] =  a[1][2]
			for i in range(0,len(fields)):	#update the strand information
				if fields[i] in ['+','-']:
					fields[i] = a[1][3]

			if outfile is None:
				print(line + '\t->\t' + '\t'.join([str(i) for i in fields]) + "\tmap_ratio=1.0000")
			else:
				print('\t'.join([str(i) for i in fields]) + "\tmap_ratio=1.0000", file=FILE_OUT)

		#when a is an even number but bigger than 2, each segment is 100% match,
		# but the whole region is not. In this case, check *min_ratio* of the query
		if len(a) > 2 :
			a_query =  a[::2] #EVEN: [('chr1', 246974830, 246974833, '+'), ('chr1', 247024833, 247024835, '+')]
			a_query_mapped_nt = sum([i[2]-i[1] for i in a_query]) #sum([3,2])
			a_target = a[1::2] #ODDS: [('chr1', 248908207, 248908210, '+'), ('chr1', 249058210, 249058212, '+')]
			a_target_chroms = set([i[0] for i in a_target])
			a_target_chroms = set([i[0] for i in a_target])
			a_target_starts = [i[1] for i in a_target]
			a_target_ends = [i[2] for i in a_target]
			#print (a_target_ends)
			map_ratio = a_query_mapped_nt/total_query_length

			#map_ratio > cutoff
			if map_ratio >= min_ratio:
				if len(a_target_chroms) == 1:
					t_chrom = a_target_chroms.pop()
					fields[0] = t_chrom
					fields[1] = min(a_target_starts)
					fields[2] = max(a_target_ends)
					if outfile is None:
						print(line + '\t->\t' + '\t'.join([str(i) for i in fields]) + ("\tmap_ratio=%.4f" % map_ratio))
					else:
						print('\t'.join([str(i) for i in fields]) + ("\tmap_ratio=%.4f" % map_ratio), file=FILE_OUT)
				else:
					if outfile is None: print(line + '\tFail\tCrossChroms')
					else: print(line + '\tFail\tCrossChroms', file=UNMAP)
			# map_ratio > 0 but < cutoff
			elif map_ratio >0 and map_ratio < min_ratio:
				if outfile is None: print(line + '\tFail' + ("\tmap_ratio=%.4f" % map_ratio))
				else: print(line + '\tFail' + ("\tmap_ratio=%.4f" % map_ratio), file=UNMAP)
