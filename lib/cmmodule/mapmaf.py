import os
import pysam
import datetime
from cmmodule  import ireader
from cmmodule.utils import printlog,update_chromID,revcomp_DNA
from cmmodule.utils import map_coordinates
from cmmodule.meta_data import __version__


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