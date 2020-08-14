import sys
import pysam
from cmmodule.utils import printlog,update_chromID,revcomp_DNA
from cmmodule.utils import map_coordinates
from cmmodule  import sam_header
from cmmodule.meta_data import __version__

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
