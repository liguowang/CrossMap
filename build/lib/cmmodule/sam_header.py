#def sam_header_generator(orig_header, chrom_size,prog_name,prog_ver,format_ver=1.0,sort_type = 'coordinate'):
#	'''generates header section for SAM file
#	chrom_size: dictionary of chromosome sizes
#	prog_name: program name
#	prog_ver: program version
#	assembly: genome assembly identifier
#	species: Species
#	sort_type: 'unknown', 'unsorted','queryname' or 'coordinate'
#	'''
#	samHeaderLine=orig_header
#	header_line = '@HD\tVN:' + str(format_ver) + '\tSO:' + sort_type
#	samHeaderLine.append(header_line)
#	for ref_id in sorted(chrom_size):
#		samHeaderLine.append('@SQ\t' + 'SN:' + ref_id + '\t' + 'LN:' + str(chrom_size[ref_id]))
#	prog_line = '@PG\t' + 'ID:' + prog_name + '\t' + 'VN:' + str(prog_ver)
#	samHeaderLine.append(prog_line)
#	
#	#for i in samHeaderLine:print i
#	return samHeaderLine

def bam_header_generator(orig_header,chrom_size,prog_name,prog_ver,co,format_ver=1.0,sort_type = 'coordinate'):
	'''generates header section for BAM file'''
	bamHeaderLine=orig_header
	name2id={}
	id = 0
	# replace 'HD'
	bamHeaderLine['HD'] = {'VN':format_ver,'SO':sort_type}
	
	# replace SQ
	tmp=[]
	for ref_name in sorted(chrom_size):
		tmp.append({'LN':chrom_size[ref_name],'SN':ref_name})
		name2id[ref_name] = id
		id += 1
	bamHeaderLine['SQ'] =  tmp
	if 'PG' in bamHeaderLine:
		bamHeaderLine['PG'] .append( {'ID':prog_name,'VN':prog_ver})
	else:
		bamHeaderLine['PG'] = [{'ID':prog_name,'VN':prog_ver}]
	
	for comment in co:
		if 'CO' in bamHeaderLine:
			bamHeaderLine['CO'].append(comment)
		else:
			bamHeaderLine['CO'] = [comment]
	return (bamHeaderLine, name2id)