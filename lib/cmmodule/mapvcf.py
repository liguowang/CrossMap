import os
import pysam
import re
import datetime
import subprocess
import logging
from cmmodule  import ireader
from cmmodule.utils import update_chromID,revcomp_DNA
from cmmodule.utils import map_coordinates
from cmmodule.meta_data import __version__

def crossmap_vcf_file(mapping, infile, outfile, liftoverfile, refgenome, noCompAllele = False, compress = False):
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

	noCompAllele : bool
		A logical value indicates whether to compare ref_allele to alt_allele after
		liftover. If True, the variant will be marked as "unmap" if
		ref_allele == alt_allele.
	'''

	if noCompAllele:
		logging.info("Keep variants [reference_allele == alternative_allele] ...")
	else:
		logging.info("Filter out variants [reference_allele == alternative_allele] ...")

	#index refegenome file if it hasn't been done
	if not os.path.exists(refgenome + '.fai'):
		logging.info("Creating index for: %s" % refgenome)
		pysam.faidx(refgenome)
	if os.path.getctime(refgenome + '.fai') < os.path.getctime(refgenome):
		logging.info("Index file is older than reference genome. Re-creating index for: %s" % refgenome)
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
			logging.info("Updating contig field ... ")
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

			print("##liftOverProgram=<CrossMap,version=%s,website=https://sourceforge.net/projects/crossmap>" % __version__, file=FILE_OUT)
			print("##liftOverChainFile=<%s>" % liftoverfile, file=FILE_OUT)
			print("##originalFile=<%s>" % infile, file=FILE_OUT)
			print("##targetRefGenome=<%s>" % refgenome, file=FILE_OUT)
			print("##liftOverDate=<%s>" % datetime.date.today().strftime("%B%d,%Y"), file=FILE_OUT)
			print(line, file=FILE_OUT)
			print(line, file=UNMAP)
			logging.info("Lifting over ... ")

		else:
			if line.startswith('#'):continue
			fields = str.split(line,maxsplit=7)
			total += 1

			chrom = fields[0]
			start = int(fields[1])-1	 # 0 based
			end = start + len(fields[3])

			a = map_coordinates(mapping, chrom, start, end, '+')
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
				try:
					fields[3] = refFasta.fetch(target_chr,target_start,target_end).upper()
				except:
					print (line + "\tFail(KeyError)", file=UNMAP)
					fail += 1
					continue

				# update END if any
				fields[7] = re.sub('END\=\d+','END='+str(target_end),fields[7])


				if a[1][3] == '-':
					fields[4] = revcomp_DNA(fields[4], True)

				# check if ref_allele is the same as alt_allele
				if noCompAllele:
					print('\t'.join(map(str, fields)), file=FILE_OUT)
				else:
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

	logging.info ("Total entries: %d" % total)
	logging.info ("Failed to map: %d" % fail)

	if compress:
		try:
			logging.info("Compressing \"%s\" ..." % outfile)
			subprocess.call("gzip " + outfile, shell=True)
		except:
			pass
