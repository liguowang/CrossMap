#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from textwrap import wrap
from cmmodule.meta_data import __version__

"""
Created on Fri Apr 16 13:51:10 2021

@author: m102324
"""
def print_help(lst):
	for i,j in lst:
		print('\n' + i + '\n' + '-'*len(i) + '\n' + '\n'.join(['  ' + k for k in wrap(j, width=80)]), file=sys.stderr)

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
	print_help(msg)

def viewchain_help():
	msg =[
	('Usage', "CrossMap.py  viewchain  <chain_file>"),
	('Description', ("print chain file into a human readable, tab-separated, 8-column file. The first 4 columns represent"
					" 'chrom','start','end','strand' of the source genome assembly, and the last 4 columns represent "
					" 'chrom','start','end','strand' of the target genome assembly.")),
	('Example', "CrossMap.py viewchain hg18ToHg19.over.chain.gz"),
	]
	print_help(msg)

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
	print_help(msg)

def wig_help():
	msg =[
	('Usage', "CrossMap.py  wig  <chain_file>  <input.wig>  <output_prefix>"),
	('Description', ("Convert wiggle format file. The \"chain_file\" can be regular or compressed (*.gz, *.Z, *.z, "
					"*.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file.  "
					"Both \"variableStep\" and \"fixedStep\" wiggle lines are supported. Wiggle format: "
					"http://genome.ucsc.edu/goldenPath/help/wiggle.html")),
	('Example', "CrossMap.py wig hg18ToHg19.over.chain.gz test.hg18.wig test.hg19"),
	]
	print_help(msg)

def bigwig_help():
	msg =[
	('Usage', "CrossMap.py  bigwig  <chain_file>  <input.bigwig>  <output_prefix>"),
	('Description', ("Convert BigWig format file. The \"chain_file\" can be regular or compressed "
					"(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
					"pointing to remote file. Bigwig format: http://genome.ucsc.edu/goldenPath/help/bigWig.html")),
	('Example', "CrossMap.py bigwig hg18ToHg19.over.chain.gz test.hg18.bw test.hg19"),
	]
	print_help(msg)

def maf_help():
	msg =[
	("usage","CrossMap.py maf  <chain_file>  <input.maf>  <refGenome.fa>  <build_name>  <output_file>"),
	("Description", ("Convert MAF format file. The \"chain_file\" and \"input.maf\" can be regular or compressed "
				  "(*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) "
				  "pointing to remote file. \"refGenome.fa\" is genome sequence file of *target assembly*. "
				  "\"build_name\" is the name of the *target_assembly* (eg \"GRCh38\")")),
	("Example", " CrossMap.py  maf	hg19ToHg38.over.chain.gz  test.hg19.maf	 hg38.fa  GRCh38 test.hg38.maf"),
	]
	print_help(msg)