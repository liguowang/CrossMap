import sys
import string
import random
from cmmodule  import ireader
from cmmodule.utils import map_coordinates


def crossmap_gff_file(mapping, ingff, outfile=None, cstyle='a'):
    '''
    Description
    -----------
    Convert genome coordinates (in GFF/GTF format) between assemblies.
    GFF (General Feature Format) lines have nine required fields that must
    be Tab-separated:

    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard
       feature types are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first
       base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute
       is set to 1 for this annotation data set, the score value will determine
       the level of gray in which this feature is displayed (higher numbers =
       darker gray). If there is no score value, enter ".".
    7. strand - Valid entries include '+', '-', or '.'
       (for don't know/don't care).
    8. frame - If the feature is a coding exon, frame should be a number
       between 0-2 that represents the reading frame of the first base. If the
       feature is not a coding exon, the value should be '.'.
    9. group - All lines with the same group are linked together into a
       single item.

    GTF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format4

    We do NOT check if features (exon, CDS, etc) originally belonging to the
    same gene  were converted into the same chromosome/strand.

    Parameters
    ----------
    mapping : dict
        Dictionary with source chrom name as key, IntervalTree object as value.

    ingff : file
        Input GFF/GTF file.

    outfile : str, optional
        Prefix of output files.

    cstyle : str, optional
        Chromosome ID style. Must be one of ['a', 's', 'l'], where
        'a' : as-is. The chromosome ID of the output file is in the same
            style of the input file.
        's' : short ID, such as "1", "2", "X.
        'l' : long ID, such as "chr1", "chr2", "chrX.
    '''

    if outfile is not None:
        rand_str = ''.join(
            random.choices(string.ascii_uppercase + string.digits, k=8))
        FILE_OUT = open(outfile, 'w')
        UNMAP = open(outfile + '.' + rand_str + '.unmap', 'w')

    for line in ireader.reader(ingff):
        if line.startswith(('#', 'track', 'browser', 'visibility')):
            continue
        if not line.strip():
            continue

        line = line.strip()
        fields = line.split('\t')
        try:
            # 0-based
            start = int(fields[3]) - 1
            end = int(fields[4])/1
            feature_size = end - start
        except:
            print(
                'Cannot recognize \"start\" and \"end\" coordinates. Skip '
                + line, file=sys.stderr)
            if outfile:
                print(line, file=UNMAP)
            continue
        if fields[6] not in ['+', '-', '.']:
            print('Cannot recognize \"strand\". Skip ' + line, file=sys.stderr)
            if outfile:
                print(line, file=UNMAP)
            continue

        strand = '-' if fields[6] == '-' else '+'

        chrom = fields[0]
        a = map_coordinates(
            mapping, chrom, start, end, strand, chrom_style=cstyle)

        if a is None:
            if outfile is None:
                print(line + '\tfail (no match to target assembly)')
            else:
                print(line, file=UNMAP)
            continue
        if len(a) != 2:
            if outfile is None:
                print(line + '\tfail (multpile match to target assembly)')
            else:
                print(line, file=UNMAP)
        else:
            # check if it is exact match
            if (int(a[1][2]) - int(a[1][1])) != feature_size:
                if outfile is None:
                    print(line + '\tfail (not exact match)')
                else:
                    print(line, file=UNMAP)
            fields[0] = a[1][0]            # chrom
            fields[3] = int(a[1][1]) + 1     # start, 1-based
            fields[4] = int(a[1][2])
            fields[6] = a[1][3]

            if outfile is None:
                print(line + '\t->\t' + '\t'.join([str(i) for i in fields]))
            else:
                print('\t'.join([str(i) for i in fields]), file=FILE_OUT)
