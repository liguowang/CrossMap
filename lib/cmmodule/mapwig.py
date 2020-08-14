import os
from cmmodule  import ireader
from cmmodule  import bgrMerge
import pyBigWig
from cmmodule.utils import printlog,update_chromID
from cmmodule.utils import map_coordinates,wiggleReader,bigwigReader

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
