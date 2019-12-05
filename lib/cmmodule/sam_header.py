
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