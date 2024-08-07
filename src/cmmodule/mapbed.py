import sys
from cmmodule import ireader
from cmmodule.utils import check_bed12
from cmmodule.utils import map_coordinates
from cmmodule import annoGene


def crossmap_bed_file(mapping, inbed,
                      outfile=None, unmapfile=None, cstyle='a', naive_bed_parsing=False):
    '''
    Convert genome coordinates (in bed format) between assemblies.
    BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

    Parameters
    ----------
    mapping : dict
        Dictionary with source chrom name as key, IntervalTree object as value.

    inbed : file
        Input BED file.

    outfile : str, optional
        Prefix of output files.

    unmapfile: str, optional
        Name of file to save unmapped entries. This option will be ignored if
        outfile is None.

    cstyle : str, optional
        Chromosome ID style. Must be one of ['a', 's', 'l'], where
        'a' : as-is. The chromosome ID of the output file is in the same style
            of the input file.
        's' : short ID, such as "1", "2", "X.
        'l' : long ID, such as "chr1", "chr2", "chrX.
        'n' : no-change. do not modify the Chromosome ID in any way.
    
    naive_bed_parsing : bool, optional
        Parse inbed with simpler assumptions (as if the file had <12 columns,
        even if it has >=12 columns)
    '''

    # check if 'outfile' was set. If not set, print to screen,
    # if set, print to file
    if outfile is not None:
        FILE_OUT = open(outfile, 'w')
        if unmapfile is not None:
            UNMAP = open(unmapfile, 'w')
        else:
            UNMAP = open(outfile + '.unmap', 'w')
    else:
        pass

    for line in ireader.reader(inbed):
        if line.startswith(('#', 'track', 'browser')):
            continue
        if not line.strip():
            continue
        line = line.strip()
        fields = line.split()
        strand = '+'

        # filter out line less than 3 columns
        if len(fields) < 3:
            print("Less than 3 fields. skip " + line, file=sys.stderr)
            if outfile:
                print(line + '\tInvalidBedFormat', file=UNMAP)
            continue
        try:
            int(fields[1])
        except:
            print(
                "Start coordinate is not an integer. skip "
                + line, file=sys.stderr)
            if outfile:
                print(line + '\tInvalidStartPosition', file=UNMAP)
            continue
        try:
            int(fields[2])
        except:
            print(
                "End coordinate is not an integer. skip "
                + line, file=sys.stderr)
            if outfile:
                print(line + '\tInvalidEndPosition', file=UNMAP)
            continue
        if int(fields[1]) > int(fields[2]):
            print(
                "\"Start\" is larger than \"End\". skip "
                + line, file=sys.stderr)
            if outfile:
                print(line + '\tStart>End', file=UNMAP)
            continue

        # deal with bed less than 12 columns, which uses less assumptions
        # about format and content of various columns. It focuses on the first
        # 3 columns (chrom, start, end). It also searches for strand
        # information by looking for a field (one of the columns) containing
        # exactly a '+' or '-'. The strand info will be kept from the latest
        # (right-most) occurance of such information. This lesser-set of
        # assumptions can also be used for files with 12+ columns if
        # naive_bed_parsing is set to True.
        if len(fields) < 12 or naive_bed_parsing:
            # try to reset strand
            try:
                for f in fields:
                    if f in ['+', '-']:
                        strand = f
            except:
                pass

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            a = map_coordinates(
                mapping, chrom, start, end, strand, chrom_style=cstyle)

            try:
                if (a is None) or (len(a) % 2 != 0):
                    if outfile is None:
                        print(line + '\tUnmap')
                    else:
                        print(line + '\tUnmap', file=UNMAP)
                    continue

                if len(a) == 2:
                    # reset fields
                    fields[0] = a[1][0]
                    fields[1] = a[1][1]
                    fields[2] = a[1][2]
                    # update the strand information
                    for i in range(0, len(fields)):
                        if fields[i] in ['+', '-']:
                            fields[i] = a[1][3]

                    if outfile is None:
                        print(
                            line + '\t->\t'
                            + '\t'.join([str(i) for i in fields]))
                    else:
                        print(
                            '\t'.join([str(i) for i in fields]),
                            file=FILE_OUT)
                if len(a) > 2:
                    count = 0
                    for j in range(1, len(a), 2):
                        count += 1
                        fields[0] = a[j][0]
                        fields[1] = a[j][1]
                        fields[2] = a[j][2]
                        # update the strand information
                        for i in range(0, len(fields)):
                            if fields[i] in ['+', '-']:
                                fields[i] = a[j][3]

                        if outfile is None:
                            print(line + '\t'
                                  + '(split.'
                                  + str(count)
                                  + ':' + ':'.join([str(i) for i in a[j-1]])
                                  + ')\t' +
                                  '\t'.join([str(i) for i in fields]))
                        else:
                            print(
                                '\t'.join([str(i) for i in fields]),
                                file=FILE_OUT)
            except:
                if outfile is None:
                    print(line + '\tFail')
                else:
                    print(line + '\tFail', file=UNMAP)
                continue

        # deal with bed12 and bed12+8 (genePred format). naive_bed_parsing must
        # be False (or we'd need an `elif` here)
        if ( len(fields) == 12 or len(fields) == 20 ) and not naive_bed_parsing:
            strand = fields[5]
            if strand not in ['+', '-']:
                raise Exception(
                    "Unknown strand: %s. Can only be '+' or '-'." % strand)
            fail_flag = False
            # [[chr,st,end],[chr,st,end],...]
            exons_old_pos = annoGene.getExonFromLine(line)
            # print exons_old_pos
            exons_new_pos = []
            for e_chr, e_start, e_end in exons_old_pos:
                a = map_coordinates(
                    mapping, e_chr, e_start, e_end, strand, chrom_style=cstyle)
                if a is None:
                    fail_flag = True
                    break

                if len(a) == 2:
                    exons_new_pos.append(a[1])
                else:
                    fail_flag = True
                    break

            if not fail_flag:
                # check if all exons were mapped to the same chromosome
                # and the same strand
                chr_id = set()
                exon_strand = set()

                for e_chr, e_start, e_end, e_strand in exons_new_pos:
                    chr_id.add(e_chr)
                    exon_strand.add(e_strand)
                if len(chr_id) != 1 or len(exon_strand) != 1:
                    fail_flag = True

                if not fail_flag:
                    # build new bed
                    cds_start_offset = int(fields[6]) - int(fields[1])
                    cds_end_offset = int(fields[2]) - int(fields[7])
                    new_chrom = exons_new_pos[0][0]
                    new_chrom_st = exons_new_pos[0][1]
                    new_chrom_end = exons_new_pos[-1][2]
                    new_name = fields[3]
                    new_score = fields[4]
                    new_strand = exons_new_pos[0][3]
                    new_thickStart = new_chrom_st + cds_start_offset
                    new_thickEnd = new_chrom_end - cds_end_offset
                    new_ittemRgb = fields[8]
                    new_blockCount = len(exons_new_pos)
                    new_blockSizes = ','.join(
                        [str(o - n) for m, n, o, p in exons_new_pos])
                    new_blockStarts = ','.join(
                        [str(n - new_chrom_st) for m, n, o, p in exons_new_pos])

                    new_bedline = '\t'.join(str(i) for i in (
                        new_chrom, new_chrom_st, new_chrom_end, new_name,
                        new_score, new_strand, new_thickStart, new_thickEnd,
                        new_ittemRgb, new_blockCount, new_blockSizes,
                        new_blockStarts))
                    if check_bed12(new_bedline) is False:
                        fail_flag = True
                    else:
                        if outfile is None:
                            print(line + '\t->\t' + new_bedline)
                        else:
                            print(new_bedline, file=FILE_OUT)

            if fail_flag:
                if outfile is None:
                    print(line + '\tFail')
                else:
                    print(line, file=UNMAP)
