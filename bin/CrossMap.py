#!/usr/bin/env python
'''
-------------------------------------------------------------------------------------
CrossMap: lift over genomic coordinates between genome assemblies.
Supports BED/BedGraph, GFF/GTF, BAM/SAM/CRAM, BigWig/Wig, VCF, and MAF format files.
-------------------------------------------------------------------------------------
'''

import sys
import optparse
import pyBigWig
#import logging
from textwrap import wrap
from cmmodule.utils     import read_chain_file
from cmmodule.mapvcf    import crossmap_vcf_file
from cmmodule.mapgvcf   import crossmap_gvcf_file
from cmmodule.mapmaf    import crossmap_maf_file
from cmmodule.mapbed    import crossmap_bed_file
from cmmodule.mapregion import crossmap_region_file
from cmmodule.mapbam    import crossmap_bam_file
from cmmodule.mapgff    import crossmap_gff_file
from cmmodule.mapwig    import crossmap_wig_file

__author__ = "Liguo Wang, Hao Zhao"
__contributor__="Liguo Wang, Hao Zhao"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPLv2"
__version__="0.5.2"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"


def general_help(cmds):
	desc=("CrossMap is a program to convert genome coordinates between different reference assemblies"
		"(e.g. from human hg19 to hg38 or vice versa). The supported file formats include BAM, BED, "
		"BigWig, CRAM, GFF, GTF, GVCF, MAF (mutation annotation format), SAM, Wiggle, and VCF.")

	print("Program: %s (v%s)" % ("CrossMap", __version__), file=sys.stderr)
	print("\nDescription: \n%s" % '\n'.join('  '+i for i in wrap(desc,width=80)), file=sys.stderr)
	print("\nUsage: CrossMap.py <command> [options]\n", file=sys.stderr)
	for k in sorted(cmds):
		print('	 ' + k + '\t' + cmds[k], file=sys.stderr)
	print(file=sys.stderr)

def bed_help():
	msg =[
	('Usage', "CrossMap.py  bed  <chain_file>  <input.bed>  [output_file]"),
	('Description', ("Convert BED format file. The \"chain_file\" and \"input.bed\" can be regular or compressed"
					"(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
					"pointing to remote file. BED format file must have at least 3 columns (chrom, start, end). "
					"If  no \"output_file\" is specified, output will be directed to the screen (console).")),
	('Example1 (write output to file)', "CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed test.hg19.bed"),
	('Example2 (write output to screen)', "CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed"),
	]
	for i,j in msg:
		print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

def gff_help():
	msg =[
	('Usage', "CrossMap.py  gff  <chain_file>  <input.gff>  <output_file>"),
	('Description', ("Convert GFF or GTF format file. The\"chain_file\" can be regular or compressed "
					"(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://)"
					" pointing to remote file. Input file must be in GFF or GTF format. GFF format: "
					"http://genome.ucsc.edu/FAQ/FAQformat.html#format3 GTF format: http://genome.ucsc.edu/FAQ/"
					"FAQformat.html#format4")),
	('Example1 (write output to file)', "CrossMap.py gff  hg19ToHg18.over.chain.gz test.hg19.gtf test.hg18.gtf"),
	('Example2 (write output to screen)', "CrossMap.py gff	hg19ToHg18.over.chain.gz test.hg19.gtf"),
	]
	for i,j in msg:
		print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

def wig_help():
	msg =[
	('Usage', "CrossMap.py  wig  <chain_file>  <input.wig>  <output_prefix>"),
	('Description', ("Convert wiggle format file. The \"chain_file\" can be regular or compressed (*.gz, *.Z, *.z, "
					"*.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file.  "
					"Both \"variableStep\" and \"fixedStep\" wiggle lines are supported. Wiggle format: "
					"http://genome.ucsc.edu/goldenPath/help/wiggle.html")),
	('Example', "CrossMap.py wig hg18ToHg19.over.chain.gz test.hg18.wig test.hg19"),
	]
	for i,j in msg:
		print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

def bigwig_help():
	msg =[
	('Usage', "CrossMap.py  bigwig  <chain_file>  <input.bigwig>  <output_prefix>"),
	('Description', ("Convert BigWig format file. The \"chain_file\" can be regular or compressed "
					"(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
					"pointing to remote file. Bigwig format: http://genome.ucsc.edu/goldenPath/help/bigWig.html")),
	('Example', "CrossMap.py bigwig hg18ToHg19.over.chain.gz test.hg18.bw test.hg19"),
	]
	for i,j in msg:
		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

#def vcf_help():
#	msg =[
#	("usage","CrossMap.py  vcf  <chain_file>  <input.vcf>  <refGenome.fa>  <output_file>"),
#	("Description", ("Convert VCF format file. The \"chain_file\" and \"input.vcf\" can be regular or compressed "
#				  "(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
#				  "pointing to remote file. \"refGenome.fa\" is genome sequence file of the *target assembly*.")),
#	("Example", " CrossMap.py vcf hg19ToHg18.over.chain.gz test.hg19.vcf hg18.fa test.hg18.vcf"),
#	]
#	for i,j in msg:
#		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)

#def gvcf_help():
#	msg =[
#	("usage","CrossMap.py gvcf <chain_file>  <input.gvcf>  <refGenome.fa>  <output_file>"),
#	("Description", ("Convert GVCF format file. The \"chain_file\" and \"input.gvcf\" can be regular or compressed "
#				  "(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
#				  "pointing to remote file. \"refGenome.fa\" is genome sequence file of the *target assembly*.")),
#	("Example", " CrossMap.py gvcf hg19ToHg18.over.chain.gz test.hg19.gvcf hg18.fa test.hg18.gvcf"),
#	]
#	for i,j in msg:
#		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)


def maf_help():
	msg =[
	("usage","CrossMap.py maf  <chain_file>  <input.maf>  <refGenome.fa>  <build_name>  <output_file>"),
	("Description", ("Convert MAF format file. The \"chain_file\" and \"input.maf\" can be regular or compressed "
				  "(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
				  "pointing to remote file. \"refGenome.fa\" is genome sequence file of *target assembly*. "
				  "\"build_name\" is the name of the *target_assembly* (eg \"GRCh38\")")),
	("Example", " CrossMap.py  maf	hg19ToHg38.over.chain.gz  test.hg19.maf	 hg38.fa  GRCh38 test.hg38.maf"),
	]
	for i,j in msg:
		 print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j,width=80)]), file=sys.stderr)


if __name__=='__main__':

	commands = {
	'bed':'convert BED, bedGraph or other BED-like files.',
	'bam':'convert BAM, CRAM or SAM format file.',
	'gff':'convert GFF or GTF format file.',
	'wig':'convert Wiggle or bedGraph format file.',
	'bigwig':'convert BigWig file.',
	'vcf':'convert VCF file.',
	'gvcf':'convert GVCF file.',
	'maf':'convert MAF (mutation annotation format) file.',
	'region':'convert big genomic regions (in BED format) such as CNV blocks.'
	}

	kwds = list(commands.keys())

	if len(sys.argv) == 1:
		general_help(commands)
		sys.exit(0)
	elif len(sys.argv) >=2:
		# deal with bed input
		if sys.argv[1].lower() == 'bed':
			if len(sys.argv) == 4:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				out_file = None
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file, print_table = False)
				crossmap_bed_file(mapTree, in_file, out_file)

			elif len(sys.argv) == 5:
				chain_file = sys.argv[2]
				in_file = sys.argv[3]
				out_file = sys.argv[4]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_bed_file(mapTree, in_file, out_file)
			else:
				bed_help()
				sys.exit(0)
		elif sys.argv[1].lower() == 'region':
			usage=("\nCrossMap.py region <chain_file>  <regions.bed> [output_file] [options]\n\nExamples:\nCrossMap.py "
				  "region hg18ToHg19.over.chain.gz CNV.hg18.bed CNV.hg19.bed   # write to file\nCrossMap.py region "
				  "hg18ToHg19.over.chain.gz CNV.hg18.bed   # write to screen")
			parser = optparse.OptionParser(usage, add_help_option=False)
			parser.add_option("-r", "--ratio", action="store",type="float",dest="min_map_ratio", default=0.85, help=
							"Minimum ratio of bases that must remap. {default=%default}" )
			(options,args)=parser.parse_args()

			if len(args) == 3:
				chain_file = args[1]
				in_file = args[2]
				out_file = None
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file, print_table = False)
				crossmap_region_file(mapTree, in_file, out_file, min_ratio = options.min_map_ratio)
			elif len(args) == 4:
				chain_file = args[1]
				in_file = args[2]
				out_file = args[3]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_region_file(mapTree, in_file, out_file, min_ratio = options.min_map_ratio)
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
			usage=("\nCrossMap.py bam  <chain_file>  <input.bam> [output_file] [options]\n\nNote:\nIf output_file is "
					"'STDOUT','-' or missing, CrossMap will write BAM file to the screen")
			parser = optparse.OptionParser(usage, add_help_option=False)
			parser.add_option("-m", "--mean", action="store",type="float",dest="insert_size", default=200.0, help=
							"Average insert size of pair-end sequencing (bp). {default=%default}")
			parser.add_option("-s", "--stdev", action="store",type="float",dest="insert_size_stdev", default=30.0, help=
							"Stanadard deviation of insert size. {default=%default}" )
			parser.add_option("-t", "--times", action="store",type="float",dest="insert_size_fold", default=3.0, help=
							"A mapped pair is considered as \"proper pair\" if both ends mapped to different strand and\
							the distance between them is less then '-t' * stdev from the mean. {default=%default}")
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
					crossmap_bam_file(mapping = mapTree, chainfile = chain_file, infile = in_file, outfile_prefix = out_file,
									chrom_size = targetChromSizes, IS_size=options.insert_size, IS_std=options.insert_size_stdev,
									fold=options.insert_size_fold,addtag=True)
				else:
					crossmap_bam_file(mapping = mapTree, chainfile = chain_file, infile = in_file, outfile_prefix = out_file,
									chrom_size = targetChromSizes, IS_size=options.insert_size, IS_std=options.insert_size_stdev,
									fold=options.insert_size_fold,addtag=False)
			else:
				parser.print_help()

		elif sys.argv[1].lower() == 'vcf':
			usage=("\nCrossMap.py vcf <chain_file>  <input.vcf>  <refGenome.fa>  <output_file> [options]\n\nExamples:\n"
				  "CrossMap.py vcf hg19ToHg18.over.chain.gz test.hg19.vcf hg18.fa test.hg18.vcf                     #comparing ref_allele to alt_allele to make sure they are different.\n"
				  "CrossMap.py vcf hg19ToHg18.over.chain.gz test.hg19.vcf hg18.fa test.hg18.vcf  --no-comp-alleles  #do NOT compare ref_allele to alt_allele.")
			parser = optparse.OptionParser(usage, add_help_option=False)
			parser.add_option("--no-comp-alleles", action="store_true",dest="no_comp_alleles", help=
							"If set, CrossMap does NOT check if the reference allele is different from the alternate allele.")
			parser.add_option("--compress", action="store_true",dest="compression", help=
							"If set, compress the output VCF file by calling the system \"gzip\".")
			(options,args)=parser.parse_args()

			if options.no_comp_alleles is None:
				options.no_comp_alleles = False

			if len(args) == 5:
				chain_file = args[1]
				in_file = args[2]
				genome_file = args[3]
				out_file = args[4]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_vcf_file(mapping = mapTree, infile= in_file, outfile = out_file, liftoverfile = chain_file, refgenome = genome_file, noCompAllele = options.no_comp_alleles, compress = options.compression)
			else:
				parser.print_help()
				sys.exit(0)

		elif sys.argv[1].lower() == 'gvcf':
			usage=("\nCrossMap.py gvcf <chain_file>  <input.gvcf>  <refGenome.fa>  <output_file> [options]\n\nExamples:\n"
				  "CrossMap.py gvcf hg19ToHg18.over.chain.gz test.hg19.gvcf hg18.fa test.hg18.gvcf                     #comparing ref_allele to alt_allele to make sure they are different.\n"
				  "CrossMap.py gvcf hg19ToHg18.over.chain.gz test.hg19.gvcf hg18.fa test.hg18.gvcf  --no-comp-alleles  #do NOT compare ref_allele to alt_allele.")
			parser = optparse.OptionParser(usage, add_help_option=False)
			parser.add_option("--no-comp-alleles", action="store_true",dest="no_comp_alleles", help=
							"If set, CrossMap does NOT check if the reference allele is different from the alternate allele." )
			parser.add_option("--compress", action="store_true",dest="compression", help=
							"If set, compress the output gVCF file by calling the system \"gzip\".")
			(options,args)=parser.parse_args()

			if options.no_comp_alleles is None:
				options.no_comp_alleles = False

			if len(args) == 5:
				chain_file = args[1]
				in_file = args[2]
				genome_file = args[3]
				out_file = args[4]
				(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file)
				crossmap_gvcf_file(mapping = mapTree, infile= in_file, outfile = out_file, liftoverfile = chain_file, refgenome = genome_file, noCompAllele = options.no_comp_alleles, compress = options.compression)
			else:
				parser.print_help()
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
			general_help(commands)
			sys.exit(0)
