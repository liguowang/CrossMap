#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
CrossMap: lift over genomic coordinates between genome assemblies.
Supports BED/BedGraph, GFF/GTF, BAM/SAM/CRAM, BigWig/Wig, VCF, and MAF format
files.
-------------------------------------------------------------------------------

Created on Sun Aug 29 18:25:14 2021
"""

import sys
import pyBigWig
import logging
import argparse
from cmmodule.utils import read_chain_file
from cmmodule.mapvcf import crossmap_vcf_file
from cmmodule.mapgvcf import crossmap_gvcf_file
from cmmodule.mapmaf import crossmap_maf_file
from cmmodule.mapbed import crossmap_bed_file
from cmmodule.mapregion import crossmap_region_file
from cmmodule.mapbam import crossmap_bam_file
from cmmodule.mapgff import crossmap_gff_file
from cmmodule.mapwig import crossmap_wig_file


__author__ = "Liguo Wang, Hao Zhao"
__contributor__ = "Liguo Wang, Hao Zhao"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPLv2"
__version__ = "0.6.5"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"


debug = False
if debug:
    logging.basicConfig(format="%(asctime)s [%(levelname)s]  %(message)s",
                        datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
else:
    logging.basicConfig(format="%(asctime)s [%(levelname)s]  %(message)s",
                        datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)


if __name__ == '__main__':
    # sub commands and help.
    commands = {
        'bed':
            'converts BED, bedGraph or other BED-like files. Only genome \
            coordinates (i.e., the first 3 columns) will be updated. Regions \
            mapped to multiple locations to the new assembly will be split. \
            Use the \"region\" command to liftover large genomic regions. \
            Use the \"wig\" command if you need bedGraph/bigWig output.',
        'bam':
            'converts BAM, CRAM, or SAM format file. Genome coordinates, \
            header section, all SAM flags, insert size will be updated.',
        'gff':
            'converts GFF or GTF format file. Genome coordinates will be \
            updated.',
        'wig':
            'converts Wiggle or bedGraph format file. Genome coordinates \
            will be updated.',
        'bigwig':
            'converts BigWig file. Genome coordinates will be updated.',
        'vcf':
            'converts VCF file. Genome coordinates, header section, reference \
            alleles will be updated.',
        'gvcf':
            'converts GVCF file. Genome coordinates, header section, \
            reference alleles will be updated.',
        'maf':
            'converts MAF (mutation annotation format) file. Genome \
            coordinates and reference alleles will be updated.',
        'region':
            'converts big genomic regions (in BED format) such as CNV blocks. \
            Genome coordinates will be updated.',
        'viewchain':
            'prints out the content of a chain file into a human readable, \
            block-to-block format.'
    }

    general_help = 'CrossMap (v%s) is a program to convert (liftover) genome \
        coordinates between different reference assemblies (e.g., from human \
        GRCh37/hg19 to GRCh38/hg38 or vice versa). Supported file formats: \
        BAM, BED, BigWig, CRAM, GFF, GTF, GVCF, MAF (mutation annotation \
        format), SAM, Wiggle, and VCF.' % (__version__)

    chromid_help = 'The style of chromosome IDs. "a" = "as-is"; "l" = \
        "long style" (eg. \"chr1\", \"chrX\"); "s" = "short style" \
        (eg. \"1\", \"X\").'
    chain_help = 'Chain file \
        (https://genome.ucsc.edu/goldenPath/help/chain.html) describes \
        pairwise alignments between two genomes. The input chain file can be \
        a plain text file or compressed (.gz, .Z, .z, .bz, .bz2, .bzip2) file.'

    # creat a parser
    parser = argparse.ArgumentParser(
        description=general_help,
        epilog='https://crossmap.readthedocs.io/en/latest/')

    # create sub-parser
    sub_parsers = parser.add_subparsers(
        help='sub-command help')

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%s %s' % ('CrossMap', __version__))

    # create the parser for the "bed" sub-command
    parser_bed = sub_parsers.add_parser(
        'bed', help=commands['bed'])
    parser_bed.add_argument(
        "chain", type=str,
        metavar="input.chain",
        help=chain_help)
    parser_bed.add_argument(
        "in_bed",
        type=str,
        metavar='input.bed',
        help="The input BED file. The first 3 columns must be “chrom”, \
            “start”, and “end”. The input BED file can be plain text file, \
            compressed file with extension of .gz, .Z, .z, .bz, .bz2 and \
            .bzip2, or even a URL pointing to accessible remote files \
            (http://, https:// and ftp://). Compressed remote files are not \
            supported.")
    parser_bed.add_argument(
        "out_bed",
        type=str,
        nargs='?',
        help="Output BED file. if argument is missing, CrossMap will write \
            BED file to the STDOUT.")
    parser_bed.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)
    parser_bed.add_argument(
        '--unmap-file',
        type=str,
        action="store",
        dest="unmap_file",
        default=None,
        help="file to save unmapped entries. This will be ignored if \
            [out_bed] was not provided.")

    # create the parser for the "bam" sub-command
    parser_bam = sub_parsers.add_parser(
        'bam',
        help=commands['bam'])
    parser_bam.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_bam.add_argument(
        "in_bam",
        type=str,
        metavar='input.bam',
        help="Input BAM file \
            (https://genome.ucsc.edu/FAQ/FAQformat.html#format5.1).")
    parser_bam.add_argument(
        "out_bam",
        type=str,
        nargs='?',
        help="Output BAM file. if argument is missing, CrossMap will write \
            BAM file to the STDOUT.")
    parser_bam.add_argument(
        "-m", "--mean",
        action="store",
        type=float,
        dest="insert_size",
        default=200,
        help="Average insert size of pair-end sequencing (bp).",)
    parser_bam.add_argument(
        "-s", "--stdev",
        action="store",
        type=float,
        dest="insert_size_stdev",
        default=30,
        help="Stanadard deviation of insert size.")
    parser_bam.add_argument(
        "-t", "--times",
        action="store",
        type=float,
        dest="insert_size_fold",
        default=3.0,
        help="A mapped pair is considered as \"proper pair\" if both ends \
            mapped to different strand and the distance between them is less \
            then '-t' * stdev from the mean.")
    parser_bam.add_argument(
        "-a", "--append-tags",
        action="store_true",
        dest="add_tags",
        help="Add tag to each alignment in BAM file. Tags for pair-end \
            alignments include: QF = QC failed, NN = both read1 and read2 \
            unmapped, NU = read1 unmapped, read2 unique mapped, NM = read1 \
            unmapped, multiple mapped, UN = read1 uniquely mapped, read2 \
            unmap, UU = both read1 and read2 uniquely mapped, UM = read1 \
            uniquely mapped, read2 multiple mapped, MN = read1 multiple \
            mapped, read2 unmapped, MU = read1 multiple mapped, read2 unique \
            mapped, MM = both read1 and read2 multiple mapped. Tags for \
            single-end alignments include: QF = QC failed, SN = unmaped, \
            SM = multiple mapped, SU = uniquely mapped.")
    parser_bam.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)

    # create the parser for the "gff" sub-command
    parser_gff = sub_parsers.add_parser(
        'gff',
        help=commands['gff'])
    parser_gff.add_argument(
        "chain", type=str,
        metavar="input.chain",
        help=chain_help)
    parser_gff.add_argument(
        "in_gff",
        type=str,
        metavar='input.gff',
        help="The input GFF \
            (General Feature Format, \
             http://genome.ucsc.edu/FAQ/FAQformat.html#format3) or GTF \
            (Gene Transfer Format, \
             http://genome.ucsc.edu/FAQ/FAQformat.html#format4) file. The \
            input GFF/GTF file can be plain text file, compressed file with \
            extension of .gz, .Z, .z, .bz, .bz2 and .bzip2, or even a URL \
            pointing to accessible remote files (http://, https:// and \
            ftp://). Compressed remote files are not supported.")
    parser_gff.add_argument(
        "out_gff",
        type=str,
        nargs='?',
        help="Output GFF/GTF file. if argument is missing, CrossMap will \
            write GFF/GTF file to the STDOUT.")
    parser_gff.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)

    # create the parser for the "wig" sub-command
    parser_wig = sub_parsers.add_parser(
        'wig',
        help=commands['wig'])
    parser_wig.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_wig.add_argument(
        "in_wig",
        type=str,
        metavar='input.wig',
        help='The input wiggle/bedGraph format file \
            (http://genome.ucsc.edu/goldenPath/help/wiggle.html). Both \
            "variableStep" and "fixedStep" wiggle lines are supported. The \
            input wiggle/bedGraph file can be plain text file, compressed \
            file with extension of .gz, .Z, .z, .bz, .bz2 and .bzip2, or \
            even a URL pointing to accessible remote files (http://, \
            https:// and ftp://). Compressed remote files are not supported.')
    parser_wig.add_argument(
        "out_wig",
        type=str,
        help='Output bedGraph file. Regardless of the input is wiggle or \
            bedGraph, the output file is always in bedGraph format.')
    parser_wig.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)

    # create the parser for the "bigwig" sub-command
    parser_bigwig = sub_parsers.add_parser(
        'bigwig',
        help=commands['bigwig'])
    parser_bigwig.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_bigwig.add_argument(
        "in_bw",
        type=str,
        metavar='input.bw',
        help='The input bigWig format file \
            (https://genome.ucsc.edu/goldenPath/help/bigWig.html).')
    parser_bigwig.add_argument(
        "out_bw",
        type=str,
        metavar='output.bw',
        help='Output bigWig file.')
    parser_bigwig.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)

    # create the parser for the "vcf" sub-command
    parser_vcf = sub_parsers.add_parser(
        'vcf',
        help=commands['vcf'])
    parser_vcf.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_vcf.add_argument(
        "in_vcf",
        type=str,
        metavar='input.vcf',
        help='Input VCF (variant call format, \
            https://samtools.github.io/hts-specs/VCFv4.2.pdf). The VCF file \
            can be plain text file, compressed file with extension of .gz, \
            .Z, .z, .bz, .bz2 and .bzip2, or even a URL pointing to \
            accessible remote files (http://, https:// and ftp://). \
            Compressed remote files are not supported.')
    parser_vcf.add_argument(
        "in_refgenome",
        type=str,
        metavar='refgenome.fa',
        help='Chromosome sequences of target assembly in FASTA \
            (https://en.wikipedia.org/wiki/FASTA_format) format.')
    parser_vcf.add_argument(
        "out_vcf",
        type=str,
        help='Output VCF file.')
    parser_vcf.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)
    parser_vcf.add_argument(
        "--no-comp-alleles",
        action="store_true",
        dest="no_comp_alleles",
        help="If set, CrossMap does NOT check if the reference allele is \
            different from the alternate allele.")
    parser_vcf.add_argument(
        "--compress",
        action="store_true",
        dest="compression",
        help='If set, compress the output VCF file by calling the \
            system \"gzip\".')

    # create the parser for the "gvcf" sub-command
    parser_gvcf = sub_parsers.add_parser(
        'gvcf',
        help=commands['gvcf'])
    parser_gvcf.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_gvcf.add_argument(
        "in_gvcf",
        type=str,
        metavar='input.gvcf',
        help='Input gVCF (genomic variant call format, \
            https://samtools.github.io/hts-specs/VCFv4.2.pdf). The gVCF file \
            can be plain text file, compressed file with extension of .gz, \
            .Z, .z, .bz, .bz2 and .bzip2, or even a URL pointing to \
            accessible remote files (http://, https:// and ftp://). \
            Compressed remote files are not supported.')
    parser_gvcf.add_argument(
        "in_refgenome",
        type=str,
        metavar='refgenome.fa',
        help='Chromosome sequences of target assembly in FASTA \
            (https://en.wikipedia.org/wiki/FASTA_format) format.')
    parser_gvcf.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)
    parser_gvcf.add_argument(
        "out_gvcf",
        type=str,
        help='Output gVCF file.')
    parser_gvcf.add_argument(
        "--no-comp-alleles",
        action="store_true",
        dest="no_comp_alleles",
        help="If set, CrossMap does NOT check if the reference allele is \
            different from the alternate allele.")
    parser_gvcf.add_argument(
        "--compress",
        action="store_true",
        dest="compression",
        help='If set, compress the output VCF file by calling the \
            system \"gzip\".')

    # create the parser for the "maf" sub-command
    parser_maf = sub_parsers.add_parser(
        'maf',
        help=commands['maf'])
    parser_maf.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_maf.add_argument(
        "in_maf",
        type=str,
        metavar='input.maf',
        help='Input MAF \
            (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) \
            format file. The MAF file can be plain text file, compressed \
            file with extension of .gz, .Z, .z, .bz, .bz2 and .bzip2, or \
            even a URL pointing to accessible remote files \
            (http://, https:// and ftp://). Compressed remote files are \
            not supported.')
    parser_maf.add_argument(
        "in_refgenome",
        type=str,
        metavar='refgenome.fa',
        help='Chromosome sequences of target assembly in FASTA \
            (https://en.wikipedia.org/wiki/FASTA_format) format.')
    parser_maf.add_argument(
        "build_name",
        type=str,
        metavar='build_name',
        help='the name of the *target_assembly* (eg "GRCh38").')
    parser_maf.add_argument(
        "out_maf",
        type=str,
        help='Output MAF file.')
    parser_maf.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)

    # create the parser for the "region" sub-command
    parser_region = sub_parsers.add_parser(
        'region',
        help=commands['region'])
    parser_region.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)
    parser_region.add_argument(
        "in_bed",
        type=str,
        metavar='input.bed',
        help="The input BED file. The first 3 columns must be “chrom”, \
            “start”, and “end”. The input BED file can be plain text file, \
            compressed file with extension of .gz, .Z, .z, .bz, .bz2 and \
            .bzip2, or even a URL pointing to accessible remote files \
            (http://, https:// and ftp://). Compressed remote files are not \
            supported.")
    parser_region.add_argument(
        "out_bed",
        type=str,
        nargs='?',
        help="Output BED file. if argument is missing, CrossMap will write \
            BED file to the STDOUT.")
    parser_region.add_argument(
        '--chromid',
        type=str,
        choices=['a', 's', 'l'],
        dest="cstyle",
        default='a',
        help=chromid_help)
    parser_region.add_argument(
        "-r", "--ratio",
        action="store",
        type=float,
        dest="min_map_ratio",
        default=0.85,
        help="Minimum ratio of bases that must remap.")

    # create the parser for the "viewchain" sub-command
    parser_viewchain = sub_parsers.add_parser(
        'viewchain',
        help=commands['viewchain'])
    parser_viewchain.add_argument(
        "chain",
        type=str,
        metavar="input.chain",
        help=chain_help)

    # print help message if no argument is provided

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    if len(sys.argv) >= 2:
        command = sys.argv[1].lower()
        args = parser.parse_args()

        if command == 'bed':
            out_file = args.out_bed
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_bed_file(
                mapping=mapTree,
                inbed=args.in_bed,
                outfile=out_file,
                unmapfile=args.unmap_file,
                cstyle=args.cstyle)

        elif command == 'bam':
            out_file = args.out_bam
            if out_file in ["STDOUT", "-"]:
                out_file = None
            print("Add tags: %s" % args.add_tags)
            print("Insert size = %f"
                  % (args.insert_size), file=sys.stderr)
            print("Insert size stdev = %f"
                  % (args.insert_size_stdev), file=sys.stderr)
            print("Number of stdev from the mean = %f"
                  % (args.insert_size_fold), file=sys.stderr)
            if args.add_tags:
                print("Add tags to each alignment = %s"
                      % (args.add_tags), file=sys.stderr)
            else:
                print("Add tags to each alignment = False", file=sys.stderr)
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_bam_file(
                mapping=mapTree,
                chainfile=args.chain,
                infile=args.in_bam,
                outfile_prefix=out_file,
                chrom_size=targetChromSizes,
                IS_size=args.insert_size,
                IS_std=args.insert_size_stdev,
                fold=args.insert_size_fold,
                addtag=args.add_tags,
                cstyle=args.cstyle)

        elif command == 'gff':
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_gff_file(
                mapping=mapTree,
                ingff=args.in_gff,
                outfile=args.out_gff,
                cstyle=args.cstyle)

        elif command == 'wig':
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_wig_file(
                mapping=mapTree,
                in_file=args.in_wig,
                out_prefix=args.out_wig,
                taget_chrom_size=targetChromSizes,
                in_format='wiggle',
                cstyle=args.cstyle)

        elif command == 'bigwig':
            # check bigwig file
            try:
                bw = pyBigWig.open(args.in_bw)
            except:
                print("\nPlease check if \"%s\" is in bigWig format!\n"
                      % args.in_bw, file=sys.stderr)
                sys.exit(0)
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_wig_file(
                mapping=mapTree,
                in_file=args.in_bw,
                out_prefix=args.out_bw,
                taget_chrom_size=targetChromSizes,
                in_format='bigwig',
                cstyle=args.cstyle)

        elif command == 'vcf':
            if args.no_comp_alleles is None:
                args.no_comp_alleles = False
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_vcf_file(
                mapping=mapTree,
                infile=args.in_vcf,
                outfile=args.out_vcf,
                liftoverfile=args.chain,
                refgenome=args.in_refgenome,
                noCompAllele=args.no_comp_alleles,
                compress=args.compression,
                cstyle=args.cstyle)

        elif command == 'gvcf':
            if args.no_comp_alleles is None:
                args.no_comp_alleles = False
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_gvcf_file(
                mapping=mapTree,
                infile=args.in_gvcf,
                outfile=args.out_gvcf,
                liftoverfile=args.chain,
                refgenome=args.in_refgenome,
                noCompAllele=args.no_comp_alleles,
                compress=args.compression,
                cstyle=args.cstyle)

        elif command == 'maf':
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_maf_file(
                mapping=mapTree,
                infile=args.in_maf,
                outfile=args.out_maf,
                liftoverfile=args.chain,
                refgenome=args.in_refgenome,
                ref_name=args.build_name,
                cstyle=args.cstyle)

        elif command == 'region':
            (mapTree, targetChromSizes, sourceChromSizes) = read_chain_file(
                args.chain)
            crossmap_region_file(
                mapping=mapTree,
                inbed=args.in_bed,
                outfile=args.out_bed,
                min_ratio=args.min_map_ratio,
                cstyle=args.cstyle)

        elif command == 'viewchain':
            chain_file = args.chain
            read_chain_file(chain_file, print_table=True)

        else:
            print("\n**Error**\nCannot recognize command: %s!\n" % command)
            parser.print_help(sys.stderr)
            sys.exit(0)
