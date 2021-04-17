.. toctree::
   :maxdepth: 2
   

.. CrossMap documentation master file, created by
   sphinx-quickstart on Thu Nov 06,  2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/logo.png
   :height: 150px
   :width: 750 px
   :scale: 50 %
   :alt: alternate text

What is CrossMap ?
====================

CrossMap is a program for genome coordinates conversion between *different assemblies*
(such as `hg18 (NCBI36) <http://www.ncbi.nlm.nih.gov/assembly/2928/>`_   <=> `hg19 (GRCh37) <http://www.ncbi.nlm.nih.gov/assembly/2758/>`_). 
It supports commonly used file formats including `BAM <https://samtools.github.io/hts-specs/SAMv1.pdf>`_, `CRAM <https://en.wikipedia.org/wiki/CRAM_(file_format)>`_, `SAM <https://en.wikipedia.org/wiki/SAM_(file_format)>`_, `Wiggle <https://genome.ucsc.edu/goldenPath/help/wiggle.html>`_, `BigWig <https://genome.ucsc.edu/goldenPath/help/bigWig.html>`_, `BED <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_, `GFF <https://genome.ucsc.edu/FAQ/FAQformat.html#format3>`_, `GTF <https://genome.ucsc.edu/FAQ/FAQformat.html#format4>`_, `MAF <https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_ `VCF <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_, and `gVCF <https://sites.google.com/site/gvcftools/home/about-gvcf>`_.


How CrossMap works?
===================

.. image:: _static/howitworks.png
   :height: 250px
   :width: 600 px
   :scale: 85 %
   :alt: alternate text


Release history
===================

**4/16/2021: Release version 0.5.3**

Add :code:`CrossMap.py viewchain` to convert chain file into block-to-block, more readable format. 

**12/08/2020: Release version 0.5.2**

Add '--no-comp-alleles' flag to :code:`CrossMap.py vcf` and :code:`CrossMap.py gvcf`. If set, CrossMap does not check if the "reference allele" is different from the "alternative allele".  

**08/19/2020: Release version 0.5.1**

In :code:`CrossMap.py region`: keep additional columns (columns after the 3rd column) of the original BED file after conversion. 

**08/14/2020: Release version 0.5.0**

Add :code:`CrossMap.py region` function to convert large genomic regions. Unlike the :code:`CrossMap.py bed` function which splits big genomic regions, :code:`CrossMap.py region` tries to convert the big genomic region as a whole. 

**07/09/2020: Release version 0.4.3**
 
Structural Variants VCF files often use INFO/END field to indicate the end of a deletion. v0.4.3 updates "END" coordinate in the INFO field. 
 
**05/04/2020: Release version 0.4.2**

Support `GVCF <https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format>`__ file conversion.
 
**03/24/2020: Release version 0.4.1**

Deal with consecutive TABs in the input MAF file.

**10/09/2019: Release version 0.3.8**

The University of California holds the copyrights in the UCSC chain files. As requested by UCSC, all UCSC generated chain files will be permanently removed from this website and the CrossMap distributions.

**07/22/2019: Release version 0.3.6**
  
1. Support MAF (mutation annotation format). 
2. Fix error "TypeError: AlignmentHeader does not support item assignment (use header.to_dict()" when lifting over BAM files. User does not need to downgrade pysam to 0.13.0 to lift over BAM files. 

**04/01/2019: Release version 0.3.4**

Fix bugs when chromosome IDs (of the source genome) in chain file do not have 'chr' prefix (such as "GRCh37ToHg19.over.chain.gz"). This version also allows CrossMap to detect if a VCF mapping was inverted, and if so reverse complements the alternative allele (Thanks to Andrew Yates). Improve wording. 

**01/07/2019: Release version 0.3.3**
 
Version 0.3.3 is exactly the same as Version 0.3.2. The reason to release this version is that CrossMap-0.3.2.tar.gz was broken when uploading to pypi.

**12/14/18: Release version 0.3.2**
 
Fix the key error problem (e.g  *KeyError: "sequence 'b'7_KI270803v1_alt'' not present"*). This error happens when a locus from the original assembly is mapped to an "alternative", "unplaced" or "unlocalized" contig in the target assembly, and this "target contig" does not exist in your target_ref.fa. In version 0.3.2, such loci will be silently skipped and saved to the ".unmap" file. 
 
**11/05/18: Release version 0.3.0**

1. v0.3.0 or newer will Support Python3. Previous versions support Python2.7.\*
2. add `pyBigWig <https://github.com/deeptools/pyBigWig>`_ as a dependency.  

Installation
==================

::

 pip3 install CrossMap	#Install CrossMap supporting Python3
 pip3 install CrossMap --upgrade	#upgrade CrossMap supporting Python3
 
 pip2 install CrossMap	#Install CrossMap supporting Python2.7.*
 pip2 install CrossMap --upgrade	#upgrade CrossMap supporting Python2.7.*


Input and Output
=================

Chain file
-----------

A `chain file <https://genome.ucsc.edu/goldenPath/help/chain.html>`_ describes a pairwise alignment between two reference assemblies. `UCSC <https://genome.ucsc.edu/>`_ and `Ensembl <https://uswest.ensembl.org/index.html>`_ chain files are available:


**UCSC chain files** 

 * Chain files from hg38 (GRCh38) to hg19 and all other organisms: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
 * Chain File from hg19 (GRCh37) to hg17/hg18/hg38 and all other organisms: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
 * Chain File from mm10 (GRCm38) to mm9 and all other organisms: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/ 

**Ensembl chain files** 

 * Human to Human: ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/
 * Mouse to Mouse: ftp://ftp.ensembl.org/pub/assembly_mapping/mus_musculus/
 * Other organisms: ftp://ftp.ensembl.org/pub/assembly_mapping/


User Input file
----------------

CrossMap supports the following file formats.
 
1. `BAM <http://samtools.sourceforge.net/SAMv1.pdf>`__, `CRAM <https://samtools.github.io/hts-specs/CRAMv3.pdf>`__, or `SAM <http://samtools.sourceforge.net/SAMv1.pdf/>`__
2. `BED <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`__ or BED-like. (BED file must have at least 'chrom', 'start', 'end')
3. `Wiggle <http://genome.ucsc.edu/goldenPath/help/wiggle.html>`__ ("variableStep", "fixedStep" and "bedGraph" formats are supported)
4. `BigWig <http://genome.ucsc.edu/goldenPath/help/bigWig.html>`__
5. `GFF <http://genome.ucsc.edu/FAQ/FAQformat.html#format3>`__ or `GTF <http://genome.ucsc.edu/FAQ/FAQformat.html#format4>`__
6. `VCF <http://vcftools.sourceforge.net/index.html>`__  
7. `GVCF <https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format>`__
8. `MAF <https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`__


Output file
-----------

The format of output files depends on the input

==============  =========================================================================================
Input_format        Output_format         
==============  =========================================================================================
BED             BED (Genome coordinates will be updated)
BAM             BAM (Genome coordinates, header section, all SAM flags, insert size will be updated)
CRAM			BAM (require pysam >= 0.8.2)
SAM             SAM (Genome coordinates, header section, all SAM flags, insert size will be updated)
Wiggle          BigWig
BigWig          BigWig
GFF		        GFF (Genome coordinates will be updated to the target assembly)
GTF             GTF (Genome coordinates will be updated to the target assembly)
VCF             VCF (header section, Genome coordinates, reference alleles will be updated)
GVCF			GVCF (header section, Genome coordinates, reference alleles will be updated)
MAF				MAF (Genome coordinates and reference alleles will be updated)
==============  =========================================================================================


Usage
=====

Run CrossMap.py **without** any arguments will print a help message
::
 
 $CrossMap.py
 
 Program: CrossMap (v0.5.2)

 Description:
   CrossMap is a program to convert genome coordinates between different reference
   assemblies(e.g. from human hg19 to hg38 or vice versa). The supported file
   formats include BAM, BED, BigWig, CRAM, GFF, GTF, GVCF, MAF (mutation annotation
   format), SAM, Wiggle, and VCF.

 Usage: CrossMap.py <command> [options]

 	 bam	convert BAM, CRAM or SAM format file.
 	 bed	convert BED, bedGraph or other BED-like files.
 	 bigwig	convert BigWig file.
 	 gff	convert GFF or GTF format file.
 	 gvcf	convert GVCF file.
 	 maf	convert MAF (mutation annotation format) file.
 	 region	convert big genomic regions (in BED format) such as CNV blocks.
 	 vcf	convert VCF file.
 	 wig	convert Wiggle or bedGraph format file.


Convert BED format files
------------------------
A `BED <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ (Browser Extensible Data) file
is a tab-delimited text file describing genome regions or gene annotations.
It consists of one line per feature, each containing 3-12 columns.
CrossMap converts BED files with less than 12 columns to a different assembly by updating the
chromosome and genome coordinates only; all other columns remain unchanged. Regions from old
assembly mapping to multiple locations to the new assembly will be split.  For 12-columns BED
files, all columns will be updated accordingly except the 4th column (name of bed line), 5th
column (score value) and 9th column (RGB value describing the display color). 12-column BED
files usually define multiple blocks (e.g. exons); if any of the exons fails to map to a new
assembly, the whole BED line is skipped. 

The input BED file can be plain text file, compressed file with extension of .gz, .Z, .z,
.bz, .bz2 and .bzip2, or even a URL pointing to accessible remote files (http://, https://
and ftp://). Compressed remote files are not supported. The output is a BED format file with
exact the same number of columns as the original one.

Standard `BED <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`__ format has 12 columns, but CrossMap also supports BED-like formats:

* BED3: The first 3 columns ("chrom", "start", "end") of BED format file.
* BED6: The first 6 columns ("chrom", "start", "end", "name", "score", "strand") of BED format file.
* Other: Format has at least 3 columns ("chrom", "start", "end") and no more than 12 columns. All other columns are arbitrary.

.. NOTE::

   1. For BED-like formats mentioned above, CrossMap only updates the "chrom", "start", "end", and "strand" columns. All other columns will be kept AS-IS.
   2.  Lines starting with '#', 'browser', 'track' will be skipped.
   3.  Lines less than 3 columns will be skipped.
   4.  The 2nd and 3rd columns must be integers.
   5.  The "+" strand is assumed if no strand information is found.
   6.  For standard BED format (12 columns). If any of the defined exon blocks cannot be uniquely mapped to target assembly, the whole entry will be skipped.
   7. The "input_chain_file" and "input_bed_file" can be regular or compressed (.gz, .Z, .z, .bz, .bz2, .bzip2) file, local file or URL (http://, https://, ftp://) pointing to remote file.
   8. If the output_file is not specified, results will be printed to screen (console). In this case, the original bed entries (include entries failed to convert) were also printed out.
   9. If the input region cannot be consecutively mapped to the target assembly, it will be split. 
   10. The \*.unmap file contains regions that cannot be unambiguously converted. 

**Example 1** 

run :code:`CrossMap bed` with **no** *output_file*::

 $ CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed3

 # Conversion results were printed to screen directly (column1-3 are hg18 based, column5-7 are hg19 based)::
 chr1	142614848	142617697	->	chr1	143903503	143906352
 chr1	142617697	142623312	->	chr1	143906355	143911970
 chr1	142623313	142623350	->	chr1	143911971	143912008

**Example 2**

run :code:`CrossMap bed` with *output_file* (test.hg19.bed3) specified::

 $ CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed3 test.hg19.bed3

 $ cat test.hg19.bed3
 chr1	143903503	143906352
 chr1	143906355	143911970
 chr1	143911971	143912008

**Example 3**

One input region was split because it cannot be mapped consecutively to the target assembly::

 $ CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed3
 
 chr10	81369946	81370453	+	->	chr10	81380000	81380507	+
 chr10	81370483	81371363	+	->	chr10	81380539	81381419	+
 chr10	81371363	81371365	+	->	chr10	62961832	62961834	+
 chr10	81371412	81371432	+	(split.1:chr10:81371412:81371422:+)	chr10	62961775	62961785	+
 chr10	81371412	81371432	+	(split.2:chr10:81371422:81371432:+)	chrX	63278348	63278358	+


**Example 4**

`BedGraph <https://genome.ucsc.edu/goldenPath/help/bedgraph.html>`_ format file can be converted using either :code:`CrossMap bed` or :code:`CrossMap wig`,
however, the output formats are different:

* Use :code:`CrossMap bed` command to convert a bedGraph file, output is a **bedGraph** file. 
* Use :code:`CrossMap wig` command to convert a bedGraph file, output is a **bigWig** file.

::

 $ CrossMap.py bed hg19ToHg38.over.chain.gz 4_hg19.bgr
 
 chrX	5873316	5873391	2.0	->	chrX	5955275	5955350	2.0
 chrX	5873673	5873710	0.8	->	chrX	5955632	5955669	0.8
 chrX	5873710	5873785	1.4	->	chrX	5955669	5955744	1.4
 chrX	5873896	5873929	0.9	->	chrX	5955855	5955888	0.9
 chrX	5873929	5874004	1.5	->	chrX	5955888	5955963	1.5
 chrX	5874230	5874471	0.3	->	chrX	5956189	5956430	0.3
 chrX	5874471	5874518	0.9	->	chrX	5956430	5956477	0.9

 $ python3 CrossMap.py wig hg19ToHg38.over.chain.gz 4_hg19.bgr output_hg38
 @ 2018-11-06 00:09:11: Read chain_file:  hg19ToHg38.over.chain.gz
 @ 2018-11-06 00:09:12: Liftover wiggle file: 4_hg19.bgr ==> output_hg38.bgr
 @ 2018-11-06 00:09:12: Merging overlapped entries in bedGraph file ...
 @ 2018-11-06 00:09:12: Sorting bedGraph file:output_hg38.bgr
 @ 2018-11-06 00:09:12: Writing header to "output_hg38.bw" ...
 @ 2018-11-06 00:09:12: Writing entries to "output_hg38.bw" ...

**Example 5** 

Use :code:`CrossMap region` command to convert large genomic regions (such as `CNV <https://en.wikipedia.org/wiki/Copy_number_variation>`_ blocks) in BED format. ::
 
 # a genomic region of 3.48Mb
 $ cat test.bed
 chr2	239716679	243199373
 
If we use :code:`CrossMap bed` command to convert this 3.48 Mb region. It will be split into 74 small blocks::

 $CrossMap.py bed GRCh37_to_GRCh38.chain.gz  test.bed
 
 chr2	239716679	243199373	(split.1:chr2:239716679:239801978:+)	chr2	238808038	238893337
 chr2	239716679	243199373	(split.2:chr2:239831978:240205681:+)	chr2	238910282	239283985
 chr2	239716679	243199373	(split.3:chr2:240205681:240319336:+)	chr2	239283986	239397641
 ... (split 74 times)
 
If we use :code:`CrossMap region` command to convert this 3.48Mb region. Note: "-r" (the minimum ratio of bases that must remap) is 0.85 by default::

 $CrossMap.py region GRCh37_to_GRCh38.chain.gz  test.bed
 
 chr2	239716679	243199373	->	chr2	238808038	242183529	map_ratio=0.9622

If we increase -r to 0.99, this region will fail::

 $CrossMap.py region GRCh37_to_GRCh38.chain.gz  test.bed -r 0.99
 
 chr2	239716679	243199373	Fail	map_ratio=0.9622
 
 
.. _bam_conversion:

Convert BAM/CRAM/SAM format files
---------------------------------
`SAM <http://samtools.sourceforge.net/samtools.shtml#5>`_ (Sequence Alignment Map) format
is a generic format for storing sequencing alignments, and BAM is the binary and compressed
version of SAM (`Li et al., 2009 <http://bioinformatics.oxfordjournals.org/content/25/16/2078.full>`_).
`CRAM <https://samtools.github.io/hts-specs/CRAMv3.pdf>`_ was designed to be an efficient reference-based
alternative to the `SAM <http://samtools.sourceforge.net/samtools.shtml#5>`_ and BAM file formats
Most high-throughput sequencing  (HTS) alignments were in SAM/BAM format and many HTS analysis
tools work with SAM/BAM format. CrossMap updates chromosomes, genome coordinates, header
sections, and all SAM flags accordingly.  CrossMap's version number is inserted into
the header section, along with the names of the original BAM file and the chain file.  For
pair-end sequencing, insert size is also recalculated. The input BAM file should be sorted
and indexed properly using Samtools (`Li et al., 2009 <http://bioinformatics.oxfordjournals.org/content/25/16/2078.full>`_).
The output format is determined by the input format and BAM output will be sorted and indexed automatically.


Typing :code:`CrossMap.py bam` without any arguments will print help a message::

 Usage:
 CrossMap.py bam  <chain_file>  <input.bam> [output_file] [options]
 
 .. NOTE::
    If output_file is 'STDOUT','-' or missing, CrossMap will write the BAM file to Stdout
 
 Options:
   -m INSERT_SIZE, --mean=INSERT_SIZE
                         Average insert size of pair-end sequencing (bp).
                         {default=200.0}
   -s INSERT_SIZE_STDEV, --stdev=INSERT_SIZE_STDEV
                        Stanadard deviation of insert size. {default=30.0}
   -t INSERT_SIZE_FOLD, --times=INSERT_SIZE_FOLD
                         A mapped pair is considered as "proper pair" if both
                         ends mapped to different strand and
                         the distance between them is less than '-t' * stdev
                         from the mean. {default=3.0}
   -a, --append-tags     Add tag to each alignment.
  
**Example**

Convert BAM from hg19 to hg18::

 # add optional tags using '-a' (recommend always use '-a' option)
 
 $ CrossMap.py bam -a ../data/hg19ToHg18.over.chain.gz test.hg19.bam test.hg18		
 Insert size = 200.000000
 Insert size stdev = 30.000000
 Number of stdev from the mean = 3.000000
 Add tags to each alignment = True
 @ 2016-10-07 15:29:06: Read chain_file:  ../data/hg19ToHg18.over.chain.gz
 @ 2016-10-07 15:29:07: Liftover BAM file: test.hg19.bam ==> test.hg18.bam
 @ 2016-10-07 15:29:14: Done!
 @ 2016-10-07 15:29:14: Sort "test.hg18.bam" ...
 @ 2016-10-07 15:29:15: Index "test.hg18.sorted.bam" ...
 Total alignments:99914
	QC failed: 0
	R1 unique, R2 unique (UU): 96094
	R1 unique, R2 unmapp (UN): 3579
	R1 unique, R2 multiple (UM): 0
	R1 multiple, R2 multiple (MM): 0
	R1 multiple, R2 unique (MU): 233
	R1 multiple, R2 unmapped (MN): 8
	R1 unmap, R2 unmap (NN): 0
	R1 unmap, R2 unique (NU): 0
	R1 unmap, R2 multiple (NM): 0
  
  
  
# BAM/SAM header sections was updated::

 $ samtools view -H  test.hg19.bam 
 @SQ	SN:chr1	LN:249250621
 @SQ	SN:chr2	LN:243199373
 @SQ	SN:chr3	LN:198022430
 ...
 @SQ	SN:chrX	LN:155270560
 @SQ	SN:chrY	LN:59373566
 @SQ	SN:chrM	LN:16571
 @RG	ID:Sample_618545BE	SM:Sample_618545BE	LB:Sample_618545BE	PL:Illumina
 @PG	ID:bwa	PN:bwa	VN:0.6.2-r126

 $ samtools view -H  test.hg18.bam
 @HD	VN:1.0	SO:coordinate
 @SQ	SN:chr1	LN:247249719
 @SQ	SN:chr10	LN:135374737
 @SQ	SN:chr11	LN:134452384
 ...
 @SQ	SN:chrX	LN:154913754
 @SQ	SN:chrX_random	LN:1719168
 @SQ	SN:chrY	LN:57772954
 @RG	ID:Sample_618545BE	SM:Sample_618545BE	LB:Sample_618545BE	PL:Illumina
 @PG	PN:bwa	ID:bwa	VN:0.6.2-r126
 @PG	ID:CrossMap	VN:0.5.0
 @CO	Liftover from original BAM/SAM file: test.hg19.bam
 @CO	Liftover is based on the chain file: ../test/hg19ToHg18.over.chain.gz 


**Optional tags:**

Q
  QC. QC failed.
N
  Unmapped. Originally unmapped or originally mapped but failed to lift over to new assembly.
M
  Multiple mapped. Alignment can be lifted over to multiple places.
U
  Unique mapped. Alignment can be lifted over to only 1 place.
		
**Tags for pair-end sequencing include:**
		
- QF = QC failed
- NN = both read1 and read2 unmapped
- NU = read1 unmapped, read2 unique mapped
- NM = read1 unmapped, multiple mapped
- UN = read1 uniquely mapped, read2 unmap
- UU = both read1 and read2 uniquely mapped
- UM = read1 uniquely mapped, read2 multiple mapped
- MN = read1 multiple mapped, read2 unmapped
- MU = read1 multiple mapped, read2 unique mapped
- MM = both read1 and read2 multiple mapped
		
**Tags for single-end sequencing include:**
		
- QF = QC failed
- SN = unmaped
- SM = multiple mapped
- SU = uniquely mapped

                         
.. note::

   1. All alignments (mapped, partial mapped, unmapped, QC failed) will write to one file. Users can filter them by tags.
   2. The header section will be updated to the target assembly.
   3. Genome coordinates and all SAM flags in the alignment section will be updated to the target assembly.
   4. If the input is a CRAM file, pysam version should >= 0.8.2 
   5. Optional fields in the alignment section will not update.

Convert Wiggle format files
---------------------------
`Wiggle <http://genome.ucsc.edu/goldenPath/help/wiggle.html>`_ (WIG) format is useful for
displaying continuous data such as GC content and the reads intensity of high-throughput sequencing data.
BigWig is a self-indexed binary-format Wiggle file and has the advantage of supporting random access.
Input wiggle data can be in variableStep (for data with irregular intervals) or fixedStep
(for data with regular intervals). Regardless of the input, the output files are always in bedGraph
format. We export files in bedGraph format because it's more compact than wiggle format, and more importantly,
CrossMap internally transforms wiggle into bedGraph to increase running speed.
 

Typing :code:`CrossMap.py wig` without any arguments will print a help message::
 
 $ CrossMap.py  wig
 
 Usage
 -----
   CrossMap.py  wig  <chain_file>  <input.wig>  <output_prefix>
 
 Description
 -----------
   Convert wiggle format file. The "chain_file" can be regular or compressed (*.gz,
   *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://,
   ftp://) pointing to remote file.  Both "variableStep" and "fixedStep" wiggle
   lines are supported. Wiggle format:
   http://genome.ucsc.edu/goldenPath/help/wiggle.html
 
 Example
 -------
   CrossMap.py wig hg18ToHg19.over.chain.gz test.hg18.wig test.hg19     

.. note::

   To improve performance, this script calls `GNU "sort" <http://www.gnu.org/software/coreutils/manual/html_node/sort-invocation.html>`_ command internally. If the "sort" command is not callable, CrossMap will exit.    

Convert BigWig format files
----------------------------
If an input file is in BigWig format, the output is BigWig format if UCSC’s
'`wigToBigWig <http://hgdownload.cse.ucsc.edu/admin/exe/>`_' executable can be found;
otherwise, the output file will be in bedGraph format.

After v0.3.0, UCSC's :code:`wigToBigWig` command is no longer needed.

Typing :code:`CrossMap.py bigwig` without any arguments will print a help message::
 
 $ CrossMap.py bigwig
 
 Usage
 -----
   CrossMap.py  bigwig  <chain_file>  <input.bigwig>  <output_prefix>
 
 Description
 -----------
   Convert BigWig format file. The "chain_file" can be regular or compressed (*.gz,
   *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://, https://,
   ftp://) pointing to remote file. Bigwig format:
   http://genome.ucsc.edu/goldenPath/help/bigWig.html
 
 Example
 -------
   CrossMap.py bigwig hg18ToHg19.over.chain.gz test.hg18.bw test.hg19
    
Example (Convert BigWig file from hg18 to hg19)::

 $ python CrossMap.py bigwig  hg19ToHg18.over.chain.gz  test.hg19.bw test.hg18
 @ 2013-11-17 22:12:42: Read chain_file:  ../data/hg19ToHg18.over.chain.gz
 @ 2013-11-17 22:12:44: Liftover bigwig file: test.hg19.bw ==> test.hg18.bgr
 @ 2013-11-17 22:15:38: Merging overlapped entries in bedGraph file ...
 @ 2013-11-17 22:15:38: Sorting bedGraph file:test.hg18.bgr
 @ 2013-11-17 22:15:39: Convert wiggle to bigwig ...

.. note:: 

   To improve performance, this script calls `GNU "sort" <http://www.gnu.org/software/coreutils/manual/html_node/sort-invocation.html>`_ command internally. If "sort" command does not exist, CrossMap will exit.

  
Convert GFF/GTF format files
----------------------------
`GFF <http://genome.ucsc.edu/FAQ/FAQformat.html#format3>`_ (General Feature Format) is another
plain text file used to describe gene structure. `GTF <http://genome.ucsc.edu/FAQ/FAQformat.html#format4>`_
(Gene Transfer Format) is a refined version of GTF. The first eight fields are the same as
GFF. Plain text, compressed plain text, and URLs pointing to remote files are all supported.
Only chromosome and genome coordinates are updated. The format of the output is determined from
the input.

Typing :code:`CrossMap.py gff` without any arguments will print a help message::
 
 $ CrossMap.py  gff
  
 Usage
 -----
   CrossMap.py  gff  <chain_file>  <input.gff>  <output_file>
 
 Description
 -----------
   Convert GFF or GTF format file. The"chain_file" can be regular or compressed
   (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL (http://,
   https://, ftp://) pointing to remote file. Input file must be in GFF or GTF
   format. GFF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format3 GTF
   format: http://genome.ucsc.edu/FAQ/FAQformat.html#format4
 
 Example1 (write output to file)
 -------------------------------
   CrossMap.py gff  hg19ToHg18.over.chain.gz test.hg19.gtf test.hg18.gtf 
 
 Example2 (write output to the screen)
 ---------------------------------
   CrossMap.py gff hg19ToHg18.over.chain.gz test.hg19.gtf
  
Example (Convert GTF file from hg19 to hg18)::

 $ python CrossMap.py gff  hg19ToHg18.over.chain.gz test.hg19.gtf test.hg18.gtf
 @ 2013-11-17 20:44:47: Read chain_file:  ../data/hg19ToHg18.over.chain.gz
 
 $ head test.hg19.gtf 
 chr1	hg19_refGene	CDS	48267145	48267291	0.000000	-	0	gene_id "NM_001194986"; transcript_id "NM_001194986"; 
 chr1	hg19_refGene	exon	66081691	66081907	0.000000	+	.	gene_id "NM_002303"; transcript_id "NM_002303"; 
 chr1	hg19_refGene	CDS	145334684	145334792	0.000000	+	2	gene_id "NM_001039703"; transcript_id "NM_001039703"; 
 chr1	hg19_refGene	exon	172017752	172017890	0.000000	+	.	gene_id "NM_001136127"; transcript_id "NM_001136127"; 
 chr1	hg19_refGene	CDS	206589249	206589333	0.000000	+	2	gene_id "NM_001170637"; transcript_id "NM_001170637"; 
 chr1	hg19_refGene	exon	210573812	210574006	0.000000	+	.	gene_id "NM_001170580"; transcript_id "NM_001170580"; 
 chr1	hg19_refGene	CDS	235850249	235850347	0.000000	-	0	gene_id "NM_000081"; transcript_id "NM_000081"; 
 chr1	hg19_refGene	CDS	235880012	235880078	0.000000	-	1	gene_id "NM_000081"; transcript_id "NM_000081"; 
 chr1	hg19_refGene	exon	3417741	3417872	0.000000	-	.	gene_id "NM_001409"; transcript_id "NM_001409"; 
 chr1	hg19_refGene	exon	10190773	10190871	0.000000	+	.	gene_id "NM_006048"; transcript_id "NM_006048"; 
 
 $ head test.hg18.gtf
 chr1	hg19_refGene	CDS	48039732	48039878	0.000000	-	0	gene_id "NM_001194986"; transcript_id "NM_001194986";
 chr1	hg19_refGene	exon	65854279	65854495	0.000000	+	.	gene_id "NM_002303"; transcript_id "NM_002303";
 chr1	hg19_refGene	CDS	144046041	144046149	0.000000	+	2	gene_id "NM_001039703"; transcript_id "NM_001039703";
 chr1	hg19_refGene	exon	170284375	170284513	0.000000	+	.	gene_id "NM_001136127"; transcript_id "NM_001136127";
 chr1	hg19_refGene	CDS	204655872	204655956	0.000000	+	2	gene_id "NM_001170637"; transcript_id "NM_001170637";
 chr1	hg19_refGene	exon	208640435	208640629	0.000000	+	.	gene_id "NM_001170580"; transcript_id "NM_001170580";
 chr1	hg19_refGene	CDS	233916872	233916970	0.000000	-	0	gene_id "NM_000081"; transcript_id "NM_000081";
 chr1	hg19_refGene	CDS	233946635	233946701	0.000000	-	1	gene_id "NM_000081"; transcript_id "NM_000081";
 chr1	hg19_refGene	exon	3407601	3407732	0.000000	-	.	gene_id "NM_001409"; transcript_id "NM_001409";
 chr1	hg19_refGene	exon	10113360	10113458	0.000000	+	.	gene_id "NM_006048"; transcript_id "NM_006048"; 


.. note::

   1. Each feature  (exon, intron, UTR, etc) is processed separately and independently, and we do NOT check if features originally belonging to the same gene were converted into the same gene.
   2. If a user wants to lift over gene annotation files, use BED12 format.
   3. If no output file was specified, the output will be printed to screen (console). In this case, items that failed to convert are also printed out.
  
Convert VCF format files
-------------------------
`VCF <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_
(variant call format) is a flexible and extendable line-oriented text format developed by
the `1000 Genome Project <http://www.1000genomes.org/>`_. It is useful for representing single
nucleotide variants, indels, copy number variants, and structural variants. Chromosomes,
coordinates, and reference alleles are updated to a new assembly, and all the other fields
are not changed.

Typing :code:`CrossMap.py vcf` without any arguments will print a help message::

 $ CrossMap.py  vcf

 Usage:
 CrossMap.py vcf <chain_file>  <input.vcf>  <refGenome.fa>  <output_file> [options]
 
 Options:
   --no-comp-alleles  If set, CrossMap does NOT check if the reference allele is
                      different from the alternate allele.
  --compress          If set, compress the output VCF file by calling the
                      system "gzip".

Example: filter out variants [reference_allele == alternative_allele]::

 $ CrossMap.py  vcf  GRCh37_to_GRCh38.chain.gz  test02_hg19.vcf  hg38.fa  out.hg38.vcf
 @ 2020-12-08 22:33:16: Read the chain file:  ../data/human/GRCh37_to_GRCh38.chain.gz
 @ 2020-12-08 22:33:17: Filter out variants [reference_allele == alternative_allele] ...
 @ 2020-12-08 22:33:17: Updating contig field ...
 @ 2020-12-08 22:33:17: Lifting over ...
 @ 2020-12-08 22:33:17: Total entries: 882
 @ 2020-12-08 22:33:17: Failed to map: 2

Example: Keep variants [reference_allele == alternative_allele]. Turn on :code:`--no-comp-allele`::

 $ CrossMap.py  vcf  GRCh37_to_GRCh38.chain.gz  test02_hg19.vcf  hg38.fa  out.hg38.vcf --no-comp-allele
 @ 2020-12-08 22:36:51: Read the chain file:  ../data/human/GRCh37_to_GRCh38.chain.gz
 @ 2020-12-08 22:36:51: Keep variants [reference_allele == alternative_allele] ...
 @ 2020-12-08 22:36:51: Updating contig field ...
 @ 2020-12-08 22:36:51: Lifting over ...
 @ 2020-12-08 22:36:51: Total entries: 882
 @ 2020-12-08 22:36:51: Failed to map: 1


.. note::

   1. Genome coordinates and reference alleles will be updated to target assembly.
   2. Reference genome is the genome sequences of target assembly.
   3. If the reference genome sequence file (../database/genome/hg18.fa) was not indexed, CrossMap will automatically index it (only the first time you run CrossMap). 
   4. Output files: *output_file* and *output_file.unmap*. 
   5. In the output VCF file, whether the chromosome IDs contain "chr" or not depends on the format of the input VCF file. 


**Interpretation of Failed tags**:

* Fail(Multiple_hits) : This genomic location was mapped to two or more locations to the target assembly.
* Fail(REF==ALT) : After liftover, the reference allele is the same as the alternative allele (i.e. this is NOT an SNP/variant after liftover). In version 0.5.2, this checking can be turned off by setting '--no-comp-alleles'.
* Fail(Unmap) : Unable to map this location to the target assembly. 
* Fail(KeyError) : Unable to find the contig ID (or chromosome ID) from the reference genome sequence (of the target assembly).  

  
Convert MAF format files
-------------------------
`MAF <https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format/>`_ 
(mutation annotation format) files are tab-delimited
files that contain somatic and/or germline mutation annotations. Please do not confuse
with the `Multiple Alignment Format <https://genome.ucsc.edu/FAQ/FAQformat.html#format5>`_.


Typing :code:`CrossMap.py maf` without any arguments will print a help message::

 $ CrossMap.py  maf
 
 usage
 -----
   CrossMap.py maf  <chain_file>  <input.maf>  <refGenome.fa>  <build_name>
   <output_file>
 
 Description
 -----------
   Convert MAF format file. The "chain_file" and "input.maf" can be regular or
   compressed (*.gz, *.Z, *.z, *.bz, *.bz2, *.bzip2) file, local file or URL
   (http://, https://, ftp://) pointing to remote file. "refGenome.fa" is genome
   sequence file of *target assembly*. "build_name" is the name of the
   *target_assembly* (eg "GRCh38")
 
 Example
 -------
  CrossMap.py  maf       hg19ToHg38.over.chain.gz  test.hg19.maf  hg38.fa  GRCh38 test.hg38.maf


Convert GVCF format files
-------------------------

GVCF file format is described in `here <https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format>`_.

Typing :code:`CrossMap.py gvcf` without any arguments will print a help message::

 $ CrossMap.py  gvcf

 Usage:
 CrossMap.py gvcf <chain_file>  <input.gvcf>  <refGenome.fa>  <output_file> [options]
 
  Options:
   --no-comp-alleles  If set, CrossMap does NOT check if the reference allele is
                      different from the alternate allele.
  --compress          If set, compress the output gVCF file by calling the
                      system "gzip".

   
Example (Convert GVCF file from hg19 to hg38)::

 $ CrossMap.py  gvcf  GRCh37_to_GRCh38.chain.gz  test10_hg19.gvcf   hg38.fa  out.hg38.gvcf 
 @ 2020-12-08 22:19:44: Read the chain file:  ../data/human/GRCh37_to_GRCh38.chain.gz
 @ 2020-12-08 22:19:44: Filter out variants [reference_allele == alternative_allele] ...
 @ 2020-12-08 22:19:44: Updating contig field ...
 @ 2020-12-08 22:19:44: Lifting over ...
 @ 2020-12-08 22:19:44: Total variants: 10
 @ 2020-12-08 22:19:44: Variants failed to map: 1
 @ 2020-12-08 22:19:44: Total non-variant regions: 22
 @ 2020-12-08 22:19:44: Non-variant regions failed to map: 0

Convert large genomic regions
------------------------------


For **large genomic regions** such as CNV blocks, the :code:`CrossMap.py bed` will split each large region into smaller blocks that are 100% matched to the target assembly. 
:code:`CrossMap.py region` will NOT split large regions, instead, it will calculate the **map ratio** (i.e. {bases mapped to target genome} / {total bases in query region}). If the
**map ratio** is larger than the threshold specified by :code:`-r`, the coordinates will be converted to the target genome, otherwise, fails. 

Typing :code:`CrossMap.py region` without any arguments will print a help message::

 Usage
 -----
 CrossMap.py region <chain_file>  <regions.bed> [output_file] [options]
 
 Examples:
 CrossMap.py region hg18ToHg19.over.chain.gz CNV.hg18.bed CNV.hg19.bed   # write to file
 CrossMap.py region hg18ToHg19.over.chain.gz CNV.hg18.bed   # write to screen
 
 Options:
   -r MIN_MAP_RATIO, --ratio=MIN_MAP_RATIO
                         Minimum ratio of bases that must remap. {default=0.85}

Example::

 $CrossMap.py region  GRCh37_to_GRCh38.chain.gz test11_hg19_region.bed
 
 @ 2020-08-14 16:46:04: Read the chain file:  ../data/human/GRCh37_to_GRCh38.chain.gz
 chr1	0	2500000	->	chr1	10000	2568561	map_ratio=0.9360
 chr1	145394955	145807817	->	chr1	145627235	146040039	map_ratio=0.9994
 chr1	146527987	147394444	->	chr1	147056425	147922330	map_ratio=0.9989
 chr10	82045472	88931651	->	chr10	80285716	87171894	map_ratio=1.0000
 chr11	43940000	46020000	->	chr11	43918450	45998449	map_ratio=1.0000
 chr15	22805313	23094530	Fail	map_ratio=0.3607
 chr15	22805313	28390339	->	chr15	22598414	28145193	map_ratio=0.8967
 chr15	31080645	32462776	->	chr15	30788442	32170575	map_ratio=1.0000
 chr15	72900171	78151253	->	chr15	72607830	77858911	map_ratio=1.0000
 chr15	83219735	85722039	->	chr15	82550985	85178808	map_ratio=0.9800
 chr16	15511655	16293689	->	chr16	15417798	16199832	map_ratio=1.0000
 chr16	21950135	22431889	->	chr16	21938814	22420568	map_ratio=1.0000
 chr16	28823196	29046783	->	chr16	28811875	29035462	map_ratio=1.0000
 chr16	29650840	30200773	->	chr16	29639519	30189452	map_ratio=1.0000
 chr17	1247834	1303556	->	chr17	1344540	1400262	map_ratio=1.0000
 chr17	2496923	2588909	->	chr17	2593629	2685615	map_ratio=1.0000
 chr17	16812771	20211017	->	chr17	16909457	20307704	map_ratio=1.0000
 chr17	29107491	30265075	->	chr17	30780473	31938056	map_ratio=1.0000
 chr17	34815904	36217432	Unmap
 chr17	43705356	44164691	->	chr17	45627990	46087325	map_ratio=1.0000
 chr2	50145643	51259674	->	chr2	49918505	51032536	map_ratio=1.0000
 chr2	96742409	97677516	->	chr2	96076661	97011779	map_ratio=1.0000
 chr2	111394040	112012649	->	chr2	110636463	111255072	map_ratio=1.0000
 chr2	239716679	243199373	->	chr2	238808038	242183529	map_ratio=0.9622
 chr22	19037332	21466726	Fail	map_ratio=0.8490
 chr22	21920127	23653646	->	chr22	21565838	23311459	map_ratio=0.9996
 chr22	51113070	51171640	->	chr22	50674642	50733212	map_ratio=1.0000
 chr3	195720167	197354826	->	chr3	195993296	197627955	map_ratio=1.0000
 chr4	1552030	2091303	->	chr4	1550303	2089576	map_ratio=1.0000
 chr5	175720924	177052594	->	chr5	176293921	177625593	map_ratio=1.0000
 chr7	72744915	74142892	->	chr7	73330912	74728554	map_ratio=0.9997
 chr8	8098990	11872558	->	chr8	8241468	12015049	map_ratio=1.0000
 chr9	140513444	140730578	->	chr9	137618992	137836126	map_ratio=1.0000

.. note::

   1. Input BED file should have at least 3 columns (chrom, start, end). Additional columns will be kept as is. 


Convert chain file 
-------------------

Typing :code:`CrossMap.py viewchain` without any arguments will print a help message::

 Usage
 -----
   CrossMap.py  viewchain  <chain_file>
 
 Description
 -----------
   print chain file into a human readable, tab-separated, 8-column file. The first
   4 columns represent 'chrom','start','end','strand' of the source genome
   assembly, and the last 4 columns represent  'chrom','start','end','strand' of
   the target genome assembly.


Example::

 $CrossMap.py viewchain ../data/human/GRCh37_to_GRCh38.chain.gz >chain.tab
 $head chain.tab
 1  10000 177417   +  1  10000 177417   +
 1  227417   267719   +  1  257666   297968   +
 1  317719   471368   +  1  347968   501617   -
 1  521368   1566075  +  1  585988   1630695  +
 1  1566075  1569784  +  1  1630696  1634405  +
 1  1569784  1570918  +  1  1634408  1635542  +
 1  1570918  1570922  +  1  1635546  1635550  +
 1  1570922  1574299  +  1  1635560  1638937  +
 1  1574299  1583669  +  1  1638938  1648308  +
 1  1583669  1583878  +  1  1648309  1648518  +


Compare to UCSC liftover tool
==============================

To access the accuracy of CrossMap, we randomly generated 10,000 genome intervals (download from `here <https://sourceforge.net/projects/crossmap/files/hg19.rand.bed.gz/download>`_) with the
fixed interval size of 200 bp from hg19. Then we converted them into hg18 using CrossMap
and `UCSC liftover tool <http://genome.ucsc.edu/cgi-bin/hgLiftOver>`_ with default configurations. We compare CrossMap
to `UCSC liftover tool <http://genome.ucsc.edu/cgi-bin/hgLiftOver>`_ because it is the most widely
used tool to convert genome coordinates.

CrossMap failed to convert 613 intervals, and the UCSC liftover tool failed to convert 614
intervals. All failed intervals are exactly the same except for one region (chr2 90542908 90543108).
UCSC failed to convert it because this region needs to be split twice:

==========================   ===========================   ====================================
Original (hg19)              Split (hg19)                  Target (hg18)
==========================   ===========================   ====================================
chr2 90542908  90543108 -    chr2 90542908 90542933 -      chr2    89906445        89906470 -
chr2 90542908  90543108 -    chr2 90542933 90543001 -      chr2    87414583        87414651 -
chr2 90542908  90543108 -    chr2 90543010 90543108 -      chr2    87414276        87414374 -
==========================   ===========================   ====================================

For genome intervals that were successfully converted to hg18, the start and end coordinates are
exactly the same between UCSC conversion and CrossMap conversion.

.. image:: _static/CrossMap_vs_UCSC.png
   :height: 400 px
   :width: 700 px
   :scale: 100 %
   :alt: CrossMap_vs_UCSC_liftover.png
   
   
Citation
=========
`Zhao, H., Sun, Z., Wang, J., Huang, H., Kocher, J.-P., & Wang, L. (2013). CrossMap: a versatile tool for coordinate conversion between genome assemblies. Bioinformatics (Oxford, England), btt730 <https://pubmed.ncbi.nlm.nih.gov/24351709/>`_   

LICENSE
==========
CrossMap is distributed under `GNU General Public License <http://www.gnu.org/copyleft/gpl.html>`_

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version. This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details. You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA


Contact                        
=======

* Wang.Liguo AT mayo.edu

