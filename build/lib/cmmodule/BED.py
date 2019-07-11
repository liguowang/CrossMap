#!/usr/bin/env python
'''
This python3 module was converted from python2.7 code using 2to3
'''
#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
from operator import itemgetter
from itertools import groupby


#import third-party modules
import psyco_full
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from itertools import *
from . import fasta
import pysam

#changes to the paths

#changing history to this module
#05/26/2011: suppport multiple spliced mapped reads

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="3.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


class ParseBED:
    '''manipulate BED (http://genome.ucsc.edu/FAQ/FAQformat.html) format file.'''
    
    def __init__(self,bedFile):
        '''This is constructor of ParseBED'''
        self.transtab = str.maketrans('ACGTNX', 'TGCANX')
        self.f=open(bedFile,'r')
        self.fileName=os.path.basename(bedFile)
        self.ABS_fileName=bedFile
        
    def groupingBED(infile,outfile=None,strand=True,boundary="utr"):
        '''Group overlapping bed entries together. When strand=Ture, overlapping bed entries within
        the same strand will be merged. When strand=False, all overlapping bed entries are merged,
        regardless of the strand. When boundary="utr" (recommended), UTR regions are considered in
        determining overlapping relationship. When boundary="cds", UTR resgions are NOT considered
        in determing overlapping relationship'''
        
        f=open(infile,'r')
        
        bed_line=re.compile(r'^\S+\s+\d+\s+\d+')
        groups = collections.defaultdict(list)  #key is internal ID, value is list of bed lines
        groups_boundary_st = {}     #key is internal ID, value is leftmost boundary
        groups_boundary_end = {}    #key is internal ID, value is rightmost boundary
        Orig_bedNum=0
        cluster_Num=0
        
        #check if the input bed file is properly sorted
        print("check if input bed file is sorted properly ...", end=' ', file=sys.stderr)
        i=0
        for line in f:
            line=line.rstrip('\r\n')
            line=line.lstrip()      
            if bed_line.match(line):    
                fields=line.split()
                i +=1
                if i==1:
                    chrom = fields[0]
                    txStart=int(fields[1])
                if i>1:
                    if ((fields[0] < chrom) or ((fields[0] == chrom) and (int(fields[1]) <txStart))):   #print first unsorted line if any
                        print("File not properly sorted:" + line, file=sys.stderr)
                        exit(1)
                    chrom = fields[0]
                    txStart=int(fields[1])
        else:   #well, file seems to be OK
            print("OK!", file=sys.stderr)
        f.seek(0)
        
        if strand:  #want to consider strand information. 
            for line in f:
                line=line.rstrip('\r\n')
                line=line.lstrip()  
                fields=line.split() 
                if bed_line.match(line) and len(fields)==12:                        
                    if fields[5] == '-':continue    #merge overlapping transcripts on plus strand               
                    Orig_bedNum +=1
                    if Orig_bedNum ==1: #this is first bed entry
                        overlap_flag=0  #has nothing to overlap because it's first entry
                        chrom = fields[0]
                        strand = fields[5]
                            
                    elif Orig_bedNum > 1:               #this is NOT first entry
                        overlap_flag=1                  #we suppose current line is overlapped with previous one. Unless it can prove it's NOT!
                        if fields[0] != chrom:          #not the same chromosome
                            overlap_flag=0
                        if (boundary == "utr") and  ( fields[0] == chrom) and (int(fields[1]) >= groups_boundary_end[cluster_Num]):
                            overlap_flag=0                                                                                          
                        if (boundary == "cds") and  ( fields[0] == chrom) and (int(fields[6]) > groups_boundary_end[cluster_Num]):
                            overlap_flag=0
                
                    if overlap_flag ==0:
                        cluster_Num +=1
                        groups[cluster_Num].append(line)
                        if (boundary == "utr"):
                            groups_boundary_st[cluster_Num]=int(fields[1])
                            groups_boundary_end[cluster_Num]=int(fields[2])
                        if (boundary == "cds"):
                            groups_boundary_st[cluster_Num]=int(fields[6])
                            groups_boundary_end[cluster_Num]=int(fields[7])
                        strand = fields[5]
                        chrom = fields[0]                       
                        
                    elif overlap_flag ==1:
                        groups[cluster_Num].append(line)                        
                        strand = fields[5]
                        chrom = fields[0]
                        if (boundary == "utr"):
                            groups_boundary_st[cluster_Num]=min(groups_boundary_st[cluster_Num],int(fields[1]))
                            groups_boundary_end[cluster_Num] = max(groups_boundary_end[cluster_Num],int(fields[2]))
                        if (boundary == "cds"):
                            groups_boundary_st[cluster_Num]=min(groups_boundary_st[cluster_Num],int(fields[6]))
                            groups_boundary_end[cluster_Num] = max(groups_boundary_end[cluster_Num],int(fields[7]))
                        
                else:
                    print("unknown line:" + line +'\n', file=sys.stderr)
            f.seek(0)
        
            #merge overlapping transcripts on minus strand
            Orig_bedNum=0
            for line in f:
                line=line.rstrip('\r\n')
                line=line.lstrip()  
                fields=line.split() 
                if bed_line.match(line) and len(fields)==12:                        
                    if fields[5] == '+':continue    #merge overlapping transcripts on plus strand               
                    Orig_bedNum +=1
                    if Orig_bedNum ==1: #this is first bed entry
                        overlap_flag=0  #has nothing to overlap because it's first entry
                        chrom = fields[0]
                        strand = fields[5]
                            
                    elif Orig_bedNum > 1:               #this is NOT first entry
                        overlap_flag=1                  #we suppose current line is overlapped with previous one. Unless it can prove it's NOT!
                        if fields[0] != chrom:          #not the same chromosome
                            overlap_flag=0
                        if (boundary == "utr") and  ( fields[0] == chrom) and (int(fields[1]) >= groups_boundary_end[cluster_Num]):
                            overlap_flag=0                                                                                          
                        if (boundary == "cds") and  ( fields[0] == chrom) and (int(fields[6]) > groups_boundary_end[cluster_Num]):
                            overlap_flag=0
                
                    if overlap_flag ==0:
                        cluster_Num +=1
                        groups[cluster_Num].append(line)
                        if (boundary == "utr"):
                            groups_boundary_st[cluster_Num]=int(fields[1])
                            groups_boundary_end[cluster_Num]=int(fields[2])
                        if (boundary == "cds"):
                            groups_boundary_st[cluster_Num]=int(fields[6])
                            groups_boundary_end[cluster_Num]=int(fields[7])
                        strand = fields[5]
                        chrom = fields[0]                       
                        
                    elif overlap_flag ==1:
                        groups[cluster_Num].append(line)                        
                        strand = fields[5]
                        chrom = fields[0]
                        if (boundary == "utr"):
                            groups_boundary_st[cluster_Num]=min(groups_boundary_st[cluster_Num],int(fields[1]))
                            groups_boundary_end[cluster_Num] = max(groups_boundary_end[cluster_Num],int(fields[2]))
                        if (boundary == "cds"):
                            groups_boundary_st[cluster_Num]=min(groups_boundary_st[cluster_Num],int(fields[6]))
                            groups_boundary_end[cluster_Num] = max(groups_boundary_end[cluster_Num],int(fields[7]))
                        
                else:
                    print("unknown line:" + line +'\n', file=sys.stderr)
            f.seek(0)


        else:       # NOT want to consider strand information. merge + and - together
            Orig_bedNum=0
            cluster_Num=0
            for line in f:
                line=line.rstrip('\r\n')
                line=line.lstrip()  
                fields=line.split() 
                if bed_line.match(line) and len(fields)==12:                        
                    #if fields[5] == '-':continue   #merge overlapping transcripts on plus strand               
                    Orig_bedNum +=1
                    if Orig_bedNum ==1: #this is first bed entry
                        overlap_flag=0  #has nothing to overlap because it's first entry
                        chrom = fields[0]
                        strand = fields[5]
                            
                    elif Orig_bedNum > 1:               #this is NOT first entry
                        overlap_flag=1                  #we suppose current line is overlapped with previous one. Unless it can prove it's NOT!
                        if fields[0] != chrom:          #not the same chromosome
                            overlap_flag=0
                        if (boundary == "utr") and  ( fields[0] == chrom) and (int(fields[1]) >= groups_boundary_end[cluster_Num]):
                            overlap_flag=0                                                                                          
                        if (boundary == "cds") and  ( fields[0] == chrom) and (int(fields[6]) > groups_boundary_end[cluster_Num]):
                            overlap_flag=0
                
                    if overlap_flag ==0:
                        cluster_Num +=1
                        groups[cluster_Num].append(line)
                        if (boundary == "utr"):
                            groups_boundary_st[cluster_Num]=int(fields[1])
                            groups_boundary_end[cluster_Num]=int(fields[2])
                        if (boundary == "cds"):
                            groups_boundary_st[cluster_Num]=int(fields[6])
                            groups_boundary_end[cluster_Num]=int(fields[7])
                        strand = fields[5]
                        chrom = fields[0]                       
                        
                    elif overlap_flag ==1:
                        groups[cluster_Num].append(line)                        
                        strand = fields[5]
                        chrom = fields[0]
                        if (boundary == "utr"):
                            groups_boundary_st[cluster_Num]=min(groups_boundary_st[cluster_Num],int(fields[1]))
                            groups_boundary_end[cluster_Num] = max(groups_boundary_end[cluster_Num],int(fields[2]))
                        if (boundary == "cds"):
                            groups_boundary_st[cluster_Num]=min(groups_boundary_st[cluster_Num],int(fields[6]))
                            groups_boundary_end[cluster_Num] = max(groups_boundary_end[cluster_Num],int(fields[7]))
                        
                else:
                    print("unknown line:" + line +'\n', file=sys.stderr)
            f.seek(0)
            
        try:
            FO=open(outfile,'w')
            print("Writing to " + outfile + "\n", file=sys.stderr)
            for k in sorted(groups.keys()):
                FO.write( "Group_"+str(k) + "\t" + str(groups_boundary_st[k]) +'\t' + str(groups_boundary_end[k]) + ":\n")
                for line in groups[k]:  
                    FO.write("\t"+line+"\n")
        except:
            #return (groups,groups_boundary_st,groups_boundary_end)
            return groups
    
    groupingBED=staticmethod(groupingBED)
            
    
    def complementBED(self,sizeFile,outfile=None):
        '''Return the complement regions of a bed file. Requires a genomeSizeFile that maps chromosome
        to sizes. In genomeSizeFile, each row contains two columns (separaed by white spaces): the 1st
        column is chromosome, and the 2nd column is size. An example of sizeFile: 
        chr1    247249719
        chr2    242951149
        chr3    199501827
        ...
        '''
        
        if outfile is None:
            outfile=self.fileName + ".compl.bed"
        FO=open(outfile,'w')
        
        #read sizeFile
        chrSize={}
        for line in open(sizeFile):
            line=line.rstrip(' \n')
            line=line.lstrip()
            fields = line.split()
            if ((len(fields)) != 2):continue
            chrSize[fields[0]] = int(fields[1])
        
        bitsets = binned_bitsets_from_file( open( self.ABS_fileName ) )
        
        for chrom in chrSize:
            if chrom in bitsets:
                bits = bitsets[chrom]
                bits.invert()
                length = chrSize[chrom]
                end=0
                while 1:
                    start = bits.next_set( end )
                    if start == bits.size: break
                    end = bits.next_clear( start )
                    if end > length: end = length
                    FO.write(chrom + "\t" + str(start) + "\t" + str(end) + "\n")
                    if end == length: break
            else:
                FO.write(chrom + "\t0\t" + str(chrSize[chrom]) + "\n")
        FO.close()
        self.f.seek(0)
                
    def bedToWig(self,outfile=None,log2scale=False,header=True):
        '''Transform bed file into wiggle format. Input bed must have at least 3 columns[chrom St End].
        For bed12 file, intron regions are automatically excluded.
        NOTE bed is 0-based and half-open while wiggle is 1-based.'''
        
        if outfile is None:
            outfile = self.fileName + ".wig"
        FO=open(outfile,'w')
        wig=collections.defaultdict(dict)
        headline="track type=wiggle_0 name=" + outfile + " track_label description='' visibility=full color=255,0,0"

        print("Writing wig file to \"",outfile,"\"...", file=sys.stderr)
        
        for line in self.f:
            line=line.rstrip('\r\n')
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue
            fields=line.rstrip('\n').split()
            coverReg=[]
            if len(fields)<12 and len(fields)>2 :
                chrom=fields[0]
                coverReg = list(range(int(fields[1])+1,int(fields[2])+1))
            elif len(fields)==12:
                chrom=fields[0]
                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                exon_starts = list(map(int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                for st,size in zip(exon_starts,exon_sizes):
                    coverReg.extend(list(range(int(fields[1]) + st + 1, int(fields[1]) + st + size +1)))                
            else:continue
            
            for i in coverReg:
                if i in wig[chrom]:
                    wig[chrom][i] +=1
                else:
                    wig[chrom][i]=1             
        if header:FO.write(headline + "\n")         
        for chr in sorted(wig.keys()):
            print("Writing ",chr, " ...", file=sys.stderr)
            FO.write('variableStep chrom='+chr+'\n')
            for coord in sorted(wig[chr]):
                if log2scale:
                    FO.write("%d\t%5.3f\n" % (coord,math.log(wig[chr][coord],2)))
                else:
                    FO.write("%d\t%d\n" % (coord,wig[chr][coord]))
        self.f.seek(0)
        FO.close()
            
        
    def bedToGFF(self,outfile=None):
        '''Transform bed file into GFF format. Borrowed from Galaxy with slight change'''
        input_name = self.ABS_fileName
        if outfile is None:
            output_name = self.fileName + ".GFF"
        else: output_name=outfile
        skipped_lines = 0
        first_skipped_line = 0
        out = open( output_name, 'w' )
        out.write( "##gff-version 2\n" )
        out.write( "##bed_to_gff_converter.py\n\n" )
        i = 0
        for i, line in enumerate( self.f ):
            complete_bed = False
            line = line.rstrip( '\r\n' )
            if line and not line.startswith( '#' ) and not line.startswith( 'track' ) and not line.startswith( 'browser' ):
                try:
                    elems = line.split( '\t' )
                    if len( elems ) == 12:
                        complete_bed = True
                    chrom = elems[0]
                    if complete_bed:
                        feature = "mRNA"
                    else:
                        try:
                            feature = elems[3]
                        except:
                            feature = 'feature%d' % ( i + 1 )
                    start = int( elems[1] ) + 1
                    end = int( elems[2] )
                    try:
                        score = elems[4]
                    except:
                        score = '0'
                    try:
                        strand = elems[5]
                    except:
                        strand = '+'
                    try:
                        group = elems[3]
                    except:
                        group = 'group%d' % ( i + 1 )
                    if complete_bed:
                        out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\t%s %s;\n' % ( chrom, feature, start, end, score, strand, feature, group  ) )
                    else:
                        out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\t%s;\n' % ( chrom, feature, start, end, score, strand, group  ) )
                    if complete_bed:
                        # We have all the info necessary to annotate exons for genes and mRNAs
                        block_count = int( elems[9] )
                        block_sizes = elems[10].split( ',' )
                        block_starts = elems[11].split( ',' )
                        for j in range( block_count ):
                            exon_start = int( start ) + int( block_starts[j] )
                            exon_end = exon_start + int( block_sizes[j] ) - 1
                            out.write( '%s\tbed2gff\texon\t%d\t%d\t%s\t%s\t.\texon %s;\n' % ( chrom, exon_start, exon_end, score, strand, group ) )
                except:
                    skipped_lines += 1
                    if not first_skipped_line:
                        first_skipped_line = i + 1
            else:
                skipped_lines += 1
                if not first_skipped_line:
                    first_skipped_line = i + 1
        out.close()
        info_msg = "%i lines converted to GFF version 2.  " % ( i + 1 - skipped_lines )
        if skipped_lines > 0:
            info_msg += "Skipped %d blank/comment/invalid lines starting with line #%d." %( skipped_lines, first_skipped_line )
        print(info_msg)
        
    def getUTR(self,utr=35):
        '''Extract UTR regions from input bed file (must be 12-column). output is 6-column bed format.
        When utr=35 [default], extract both 5' and 3' UTR. When utr=3, only extract 3' UTR. When utr=5,
        only extract 5' UTR'''
        
        ret_lst=[]
        for line in self.f:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue
            fields=line.rstrip('\r\n').split()
            chrom=fields[0]
            geneName=fields[3]
            strand=fields[5]
            txStart=int(fields[1])
            txEnd=int(fields[2])
            cdsStart=int(fields[6])
            cdsEnd=int(fields[7])       
            exon_start=list(map(int,fields[11].rstrip(',').split(',')))
            exon_start=list(map((lambda x: x + txStart),exon_start))
                
            exon_end=list(map(int,fields[10].rstrip(',').split(',')))
            exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
            
            if strand == '+':
                if (utr==35 or utr==5):
                    for st,end in zip(exon_start,exon_end):
                        if st < cdsStart:
                            utr_st = st
                            utr_end = min(end,cdsStart)
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])                  
                if (utr==35 or utr==3):
                    for st,end in zip(exon_start,exon_end):
                        if end > cdsEnd:
                            utr_st = max(st, cdsEnd)
                            utr_end = end
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])
            if strand == '-':
                if (utr==35 or utr==3):
                    for st,end in zip(exon_start,exon_end):
                        if st < cdsStart:
                            utr_st = st
                            utr_end = min(end,cdsStart)
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])                  
                if (utr==35 or utr==5):
                    for st,end in zip(exon_start,exon_end):
                        if end > cdsEnd:
                            utr_st = max(st, cdsEnd)
                            utr_end = end
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])
        self.f.seek(0)
        return ret_lst
            
    def getExon(self):
        '''Extract exon regions from input bed file (must be 12-column). output is 6-column Tab 
        separated bed file, each row represents one exon'''
        
        ret_lst=[]
        for f in self.f:
            f = f.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
            for base,offset in zip( blockStarts, blockSizes ):
                ret_lst.append((chrom, base, base+offset))
        self.f.seek(0)
        return ret_lst

    def getTranscriptRanges(self):
        '''Extract exon regions from input bed file (must be 12-column). Return ranges of
        transcript'''
        
        mRNA_ranges = []
        for line in self.f:
            try:
                if line.startswith('#'):continue
                if line.startswith('track'):continue
                if line.startswith('browser'):continue
                fields=line.rstrip('\r\n').split()
                txStart = int(fields[1])
                txEnd = int(fields[2])
                chrom = fields[0]
                strand = fields[5]
                geneName = fields[3]
                yield([chrom, txStart, txEnd, strand, geneName + ":" + chrom + ":" + str(txStart) + '-' + str(txEnd)])
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
        self.f.seek(0)

    def getCDSExon(self):
        
        '''Extract CDS exon regions from input bed file (must be 12-column).'''     
        ret_lst=[]
        for f in self.f:
            f = f.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
            # grab cdsStart - cdsEnd
            cds_exons = []
            cds_seq = ''
            genome_seq_index = []
            for base,offset in zip( blockStarts, blockSizes ):
                if (base + offset) < cdsStart: continue
                if base > cdsEnd: continue
                exon_start = max( base, cdsStart )
                exon_end = min( base+offset, cdsEnd ) 
                #cds_exons.append( (exon_start, exon_end) )
                ret_lst.append([chrom,exon_start,exon_end]) 
        self.f.seek(0)
        return ret_lst

    def truncate_bed(self, truncation_from=3, size = 250):
        
        '''
        truncate bed from either the 5' end or 3' end.
        truncation_from: 3 or 5
        size = truncation size (default 250 bp)
        '''     
        ret_lst=[]
        for f in self.f:
            gene_all_bases = []
            chose_bases = []
            f = f.strip().split()
            chrom = f[0]
            txStart  = int(f[1])
            txEnd    = int(f[2])
            score = f[4]
            name = f[3]
            strand = f[5]
            exon_start=list(map(int,f[11].rstrip(',').split(',')))
            exon_start=list(map((lambda x: x + txStart),exon_start))
            exon_end=list(map(int,f[10].rstrip(',').split(',')))
            exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
            for st,end in zip(exon_start,exon_end):
                gene_all_bases.extend(list(range(st+1,end+1)))  #1-based, closed
            gene_all_bases = sorted(gene_all_bases)
            
            new_chrom = chrom       #1
            new_name = name         #4
            new_score = score       #5
            new_strand = strand     #6

            
            if truncation_from == 3:
                if strand == '+':
                    chose_bases = gene_all_bases[-size:]
                elif strand == '-':
                    chose_bases = gene_all_bases[:size]
                else:
                    print('unknown strand ' + strand, file=sys.stderr)
            elif truncation_from == 5:
                if strand == '+':
                    chose_bases = gene_all_bases[:size]
                elif strand == '-':
                    chose_bases = gene_all_bases[-size:]
                else:
                    print('unknown strand ' + strand, file=sys.stderr)
            else:
                print("truncation_from takes '3' or '5'", file=sys.stderr)
            
            new_tx_start = chose_bases[0]-1 #2
            new_tx_end = chose_bases[-1]    #3
            new_thick_start = new_tx_start  #7
            new_thick_end = new_tx_end      #8
            new_igb = f[8]                  #9
            
            exon_boundarie = []
            for k, g in groupby(enumerate(chose_bases), lambda i_x:i_x[0]-i_x[1]):
                group = list(map(itemgetter(1),g))
                exon_boundarie.append((group[0],group[-1]))     #eg [(st1, end1), (st2, end2)]
            new_block_count = len(exon_boundarie)       #10
            
            tmp = []
            for i,j in (exon_boundarie):
                tmp.append(str(j-i+1))
            
            new_block_sizes = ','.join(tmp)     #11
            
            tmp = []
            for i,j in (exon_boundarie):
                tmp.append(str(i - 1 - new_tx_start))
            
            new_block_starts = ','.join(tmp)        #12
            
            ret_lst.append('\t'.join([str(i) for i in [new_chrom,new_tx_start,new_tx_end,new_name,new_score,new_strand,new_thick_start,new_thick_end,new_igb,new_block_count,new_block_sizes,new_block_starts]]))
        self.f.seek(0)
        return ret_lst

        
    def getIntron(self):
        '''Extract Intron regions from input bed file (must be 12-column).  output is 6-column Tab 
        separated bed file, each row represents one intron'''

        ret_lst=[]
        for line in self.f:
            try:
                if line.startswith('#'):continue
                if line.startswith('track'):continue
                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
                cds_start = int( fields[6] )
                cds_end   = int( fields[7] )
                if int(fields[9] ==1):
                    continue
            
                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                intron_start = exon_ends[:-1]
                intron_end=exon_starts[1:]
                
                if(strand == '-'):
                    intronNum=len(intron_start)
                    for st,end in zip(intron_start,intron_end):
                        #FO.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t0\t" + strand + '\n')
                        #intronNum -= 1
                        ret_lst.append([chrom,st,end])
                else:
                    intronNum=1
                    for st,end in zip(intron_start,intron_end):
                        #FO.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t0\t" + strand + '\n')
                        #intronNum += 1
                        ret_lst.append([chrom,st,end])
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
        self.f.seek(0)
        return ret_lst

    def getSpliceJunctions(self):
        '''
        Return splice junctions for each transcripts.
        Single exon gene will be skipped.
        '''

        for line in self.f:
            try:
                tmp = []
                if line.startswith('#'):continue
                if line.startswith('track'):continue
                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                if int(fields[9]) == 1: #skip single exon gene
                    yield ((geneName, chrom, tx_start, tx_end, tmp))
                    continue
            
                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                intron_start = exon_ends[:-1]
                intron_end=exon_starts[1:]
                
                
                for st,end in zip(intron_start,intron_end):
                    tmp.append(chrom + ':' + str(st) + '-' + str(end))
                yield ((geneName, chrom, tx_start, tx_end, tmp))
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue

    def getIntergenic(self,direction='up', size=1000):
        '''get intergenic regions. direction=up or down or both.'''
        
        ret_lst=[]
        for line in self.f:
            if line.startswith(('#','track','browser')):
                continue
            fields = line.split()
            chrom     = fields[0]
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            strand    = fields[5]
            if(direction=="up" or direction=="both"):
                if strand=='-':
                    region_st=tx_end
                    region_end=tx_end +size
                else:
                    region_st = max(tx_start-size,0)
                    region_end=tx_start
                ret_lst.append([chrom,region_st,region_end])
            if (direction=="down" or direction=="both"):
                if strand == '-':
                    region_st = max(0,tx_start-size)
                    region_end = tx_start
                else:
                    region_st = tx_end
                    region_end = tx_end+size
                ret_lst.append([chrom,region_st,region_end])
        self.f.seek(0)
        return ret_lst

        
    def getBedinfor(self,outfile=None, reffa=None):
        '''Extract information (such as exonNumber, exonSize (min,max,mean)etc,.) from bed entries.
        '''
        if reffa is not None:
            print("Indexing " + reffa + ' ...', end=' ', file=sys.stderr)
            pysam.faidx(reffa)
            print("Done!", file=sys.stderr)
        #transtab = maketrans("ACGTNX","TGCANX")
        transtab = str.maketrans('ACGTNX', 'TGCANX')
        if outfile is None:outfile = self.fileName + ".infor.xls"
        FO=open(outfile,'w')
        print("writing feature information to " + outfile + " ...", file=sys.stderr)                
        
        intron_sizes=[]
        FO.write("geneID\t" + "Exon_Num\t" + "Exon_Len_Sum\t" + "Exon_Len_Min\t" + "Exon_Len_Max\t" + "Exon_Len_Avg\t")
        FO.write("Intron_Len_Sum\t" + "Intron_Len_Min\t" + "Intron_Len_Max\t" + "Intron_Len_Avg\t" + "3'UTR\t" + "5'UTR\t" + "GC\t" + "\n")
        for line in self.f:
            intron_sizes=[]
            mRNA_seq = ''
            GC_content = 0.0
            utr5len = 0
            utr3len = 0
            try:
                if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
                cds_start = int( fields[6] )
                cds_end   = int( fields[7] )                
                exon_num= int(fields[9])
                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                exon_starts = list(map(int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                intron_start = exon_ends[:-1]
                intron_end=exon_starts[1:]
                
                if reffa is not None:
                    for st,end in zip(exon_starts, exon_ends):
                        exon_coord = chrom + ':' + str(st +1) + '-' + str(end)
                        tmp = pysam.faidx(reffa,exon_coord)
                        mRNA_seq += ''.join([i.rstrip('\n\r') for i in tmp[1:]])
                    if strand == '-':
                        mRNA_seq = mRNA_seq.upper().translate(transtab)[::-1]  
                    elif strand == '+':
                        mRNA_seq = mRNA_seq.upper()
                    GC_content = (mRNA_seq.count('C') + mRNA_seq.count('G'))*1.0/len(mRNA_seq)
                

                for st,end in zip(exon_starts,exon_ends):
                    if strand == '+':
                        if st < cds_start:  # 5' UTR
                            utr_st = st
                            utr_end = min(end,cds_start)
                            utr5len += (utr_end - utr_st)
                        if end > cds_end:   # 3' UTR
                            utr_st = max(st, cds_end)
                            utr_end = end
                            utr3len += (utr_end - utr_st)
                    if strand == '-':
                        if st < cds_start:  # 3' UTR
                            utr_st = st
                            utr_end = min(end,cds_start)
                            utr3len += (utr_end - utr_st)
                        if end > cds_end:   # 5' UTR
                            utr_st = max(st, cds_end)
                            utr_end = end
                            utr5len += (utr_end - utr_st)                               
                if cds_start == cds_end:
                    utr3len = 0
                    utr5len = 0
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
            for st,end in zip(intron_start,intron_end):
                intron_sizes.append(end-st)
                
            FO.write(chrom +':'+str(tx_start+1)+':'+str(tx_end)+ ':' + geneName +'\t')
            FO.write(str(exon_num) + "\t" + str(sum(exon_sizes)) +'\t' +str(min(exon_sizes)) + '\t' + str(max(exon_sizes)) +'\t' + str(sum(exon_sizes)/exon_num) +"\t")
            if intron_sizes:
                FO.write(str(sum(intron_sizes)) +'\t' + str(min(intron_sizes)) + '\t' + str(max(intron_sizes)) +'\t' + str(sum(intron_sizes)/len(intron_sizes)) + "\t")     
            else:
                FO.write('0' +'\t' + '0' + '\t' + '0' +'\t' + '0' + "\t")
            FO.write(str(utr3len) + '\t' + str(utr5len) + '\t' )
            if reffa is not None:
                FO.write(str(GC_content) + '\n')
            else:
                FO.write("NA" + '\n')

        self.f.seek(0)
        FO.close()  

    def filterBedbyIntronSize(self,outfile=None,min_intron=50):
        '''Filter bed files with intron size. Mamalian gene has minimum intron size
        '''
        if outfile is None:outfile = self.fileName + ".filterIntron.xls"
        FO=open(outfile,'w')
        print("writing feature information to " + outfile + " ...", file=sys.stderr)                
        
        intron_sizes=[]
        for line in self.f:
            intron_sizes=[]
            flag=0
            try:
                if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
                cds_start = int( fields[6] )
                cds_end   = int( fields[7] )                
                exon_num= int(fields[9])
                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                exon_starts = list(map(int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                intron_start = exon_ends[:-1]
                intron_end=exon_starts[1:]
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
            if exon_num <=1:    #intron size is 0
                continue
            for st,end in zip(intron_start,intron_end):
                if end-st <=min_intron:
                    flag=1
                    break
            if flag == 0:
                FO.write(line)              
            
        self.f.seek(0)
        FO.close()      


    def getAllUniqJunctions(self,outfile=None,flankSize=20):
        '''Extract unique (non-redundant) junctions from input bed file (must be 12-column).  use 
        flankSize to represent junction. Note that too large flankSize could exceed chromosome
        boundary and raise error when you convert bed to bigbed. output 12-column Tab separated bed
        file'''
        if outfile is None:outfile = self.fileName + ".uniqJunctions.bed"
        FO=open(outfile,'w')
        print("writing unique junctions " + outfile + " ...", file=sys.stderr)              
        uniqJunc = collections.defaultdict(int)
        for line in self.f:
            #try:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue   
            fields = line.split()
            chrom     = fields[0]
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            geneName      = fields[3]
            strand    = fields[5].replace(" ","_")
            cds_start = int( fields[6] )
            cds_end   = int( fields[7] )
            if int(fields[9] ==1):
                continue
            
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
                
            for st,end in zip(intron_start,intron_end):
                key=chrom + str(st) + str(end)
                if key not in uniqJunc:
                    FO.write(chrom + "\t" + str(st-flankSize) + "\t" + str(end+flankSize) + "\t" + geneName  + "\t0\t" + strand + '\t' + str(st-flankSize) + "\t" + str(end+flankSize) + '\t0,0,0\t2\t' + str(flankSize)+','+str(flankSize)+'\t' + '0,' + str(end-st+flankSize) +'\n')
                uniqJunc[key]+=1
            #except:
            #   print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
            #   continue
        self.f.seek(0)
        FO.close()  

    def collapseJunctionBed(self,outfile=None):
        '''Junctions spannig the same block will be merged. multiple spliced junctions will be reported
        as is without any changes'''
        
        if outfile is None:outfile = self.fileName + ".collapsed.bed"
        FO=open(outfile,'w')
        print("\tCollapse junctions for " +     self.fileName, end=' ', file=sys.stderr)    
        
        tss_start =collections.defaultdict(list)
        tss_end =collections.defaultdict(list)
        exon_size1 =collections.defaultdict(list)
        exon_size2 =collections.defaultdict(list)
        block_start1 =collections.defaultdict(list)
        block_start2 =collections.defaultdict(list)
        count = collections.defaultdict(int)
        chrm=dict()
        strands=dict()
        for line in self.f:
            skip_flag=0
            try:
                if line.startswith(('#','track','broser')):continue 
                fields = line.rstrip().split()
                exonSize = list(map(int, fields[10].rstrip(',\n').split(',')))
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName  = fields[3]
                score = int(fields[4])
                strand    = fields[5].replace(" ","_")
                if int(fields[9]) ==1:  #if bed has only 1 block
                    print(line, end=' ', file=FO)
                    continue
                if int(fields[9]) >=3:  #if bed has more than 3 blocks
                    print(line, end=' ', file=FO)
                    continue
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
            exonSize = list(map(int,fields[10].rstrip(',\n').split(',')))
            block_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
                
            intronKey = chrom + ":" + str(intron_start[0]+1) + "-" + str(intron_end[0]) + ":" + strand
            count[intronKey] += score
            chrm[intronKey] = chrom
            strands[intronKey] = strand
            tss_start[intronKey].append(tx_start)
            tss_end[intronKey].append(tx_end)
            exon_size1[intronKey].append(exonSize[0])
            exon_size2[intronKey].append(exonSize[1])
            block_start1[intronKey].append(block_starts[0])
            block_start2[intronKey].append(block_starts[1])
            
        print("Writing junctions to " + outfile + " ...", file=sys.stderr)
        for key in count:
            print(chrm[key] + '\t' + str(min(tss_start[key])) + '\t' + str(max(tss_end[key])) + '\t' + "SR=" + str(count[key]) + '\t' + str(count[key]), end=' ', file=FO)
            print('\t'+strands[key] + '\t' + str(min(tss_start[key])) + '\t' + str(max(tss_end[key])) + '\t', end=' ', file=FO)
            print('255,0,0' + '\t2\t' + str(max(exon_size1[key])) + ',' + str(max(exon_size2[key])) + '\t', end=' ', file=FO)
            print(str(min(block_start1[key])) + ',' + str(max(block_start2[key])), file=FO)
            
        self.f.seek(0)
        FO.close()
        
    def filterJunctionBed(self,outfile=None,overhang=8,min_intron=50,max_intron=1000000,cvg=1):
        '''filter junction bed file according to overhang size and supporting read'''
        
        if outfile is None:outfile = self.fileName + ".filter.bed"
        FO=open(outfile,'w')
        print("\tfilter junctions ... ", file=sys.stderr)               

        for line in self.f:
            skip_flag=0
            if line.startswith(('#','track','broser','\n')):continue 
            fields = line.rstrip().split()
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            exonSize = list(map(int,fields[10].rstrip(',\n').split(',')))
            block_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
            if fields[9] >=3:   #pass for multiple spliced read
                print(line, end=' ', file=FO)
                continue
            if fields[4]>cvg:       #pass for multi-read supporting junction
                print(line, end=' ', file=FO)
                continue                
            for i in exonSize:
                if i < overhang:    #filter out short overhang junctions
                    skip_flag=1
                    break
            for st,end in zip(intron_start,intron_end):
                if (end - st) > max_intron or (end -st) < min_intron:
                    skip_flag=1
                    break
            if skip_flag ==1:
                continue
            else:
                print(line, file=FO)
        self.f.seek(0)
        FO.close()      
        
    def unionBED(self,outfile=None,outNameFile=None,boundary="utr",stranded=True):
        '''Collapse bed entries through UNION all overlapping exons. Just like "dense" display mode 
        in UCSC genome browser. Bed entries will be merged if the following conditions are met:
        1) Coordinates overlapped
        2) on the same chromosome
        output two files, one bed file contains the merged bed entries, the other *.name.txt file 
        contains which genes were merged together.
        NOTE:
        1) input bed should 12-columns and sorted by chrom then by start position 
           (use "sort -k1,1 -k2,2n myfile.bed >myfile.sorted.bed")
        2) when stranded=True, only overlapping bed entries on the same strand will be merged. 
           Otherwise, all overlapping bed entries will be meraged without considering strand information.
        '''
        
        #open output file
        if outfile is None:outfile = self.fileName + ".merge.unionExon.bed"
        if outNameFile is None: outNameFile = self.fileName + ".name.txt"       
        FO=open(outfile,'w')
        FName=open(outNameFile,'w') 
        
        #some local variable within function
        bed_line=re.compile(r'^\S+\s+\d+\s+\d+')
        Merge_bed_TxStart=collections.defaultdict(list)
        Merge_bed_TxEnd=collections.defaultdict(list)
        Merge_bed_cdsStart=collections.defaultdict(list)
        Merge_bed_cdsEnd=collections.defaultdict(list)
        Merge_bed_ExonStart=collections.defaultdict(list)
        Merge_bed_ExonEnd=collections.defaultdict(list)
        Merge_bed_geneName=collections.defaultdict(list)
        Orig_bedNum=0   #number of bed entries before merging
        Merge_bedNum=0  #number of bed entries after merging
        Merge_bed_chr={}
        Merge_bed_strand={}
        final_exon_starts=[]
        final_exon_sizes=[]


        #check if the input bed file is properly sorted
        print("check if input bed file is sorted properly ...", end=' ', file=sys.stderr)
        i=0
        for line in self.f:
            line=line.strip()
            if bed_line.match(line):    
                fields=line.split()
                i +=1
                if i==1:
                    chrom = fields[0]
                    txStart=int(fields[1])
                if i>1:
                    if ((fields[0] < chrom) or ((fields[0] == chrom) and (int(fields[1]) <txStart))):   #print first unsorted line if any
                        print("File not properly sorted:" + line, file=sys.stderr)
                        exit(1)
                    chrom = fields[0]
                    txStart=int(fields[1])
        else:   #well, file seems to be OK
            print("OK!", file=sys.stderr)
        self.f.seek(0)


        
        print("merge bed file to " + outfile + " ...", file=sys.stderr) 
        FO.write("track name=" + self.fileName + " description=" + ' \"' + "Overlapping entries in " + self.fileName + " were merged" + '\"'+'\n') 
        if stranded:    #want to consider strand information
            for line in self.f:
                line=line.strip()
                if bed_line.match(line):
                    fields=line.split() 
                    if(len(fields)!=12):
                        print("bed file must be 12 columns speparated by Tab(whilte space)" +line, file=sys.stderr)
                    if fields[5] == '-':continue    #merge overlapping transcripts on plus strand               
                    Orig_bedNum +=1
                    if Orig_bedNum ==1: #this is first bed entry
                        overlap_flag=0  #has nothing to overlap because it's first entry
                        chrom = fields[0]
                        strand = fields[5]
                        txStart= int(fields[1])
                        txEnd= int(fields[2])
                        cdsStart= int(fields[6])
                        cdsEnd= int(fields[7])      
                        geneName=fields[3]
                        score=fields[4] 
                        txEnd_float=txEnd
                        cdsEnd_float=cdsEnd
                    elif Orig_bedNum > 1:               #this is NOT first entry
                        overlap_flag=1                  #we suppose current line is overlapped with previous one. Unless it can prove it's NOT!
                        if fields[0] != chrom:          #not the same chromosome
                            overlap_flag=0
                            txEnd_float = int(fields[2])
                        if (boundary == "utr") and  ( fields[0] == chrom) and (int(fields[1]) >= txEnd_float):
                            overlap_flag=0                                                                                          
                        if ((boundary == "cds") and  ( fields[0] == chrom) and (int(fields[6]) > cdsEnd_float)):
                            overlap_flag=0
                        txStart= int(fields[1])
                        txEnd= int(fields[2])
                        strand = fields[5]
                        chrom = fields[0]               
                        cdsStart= int(fields[6])
                        cdsEnd= int(fields[7])      
                        geneName=fields[3]
                        score=fields[4] 
                        txEnd_float=max(txEnd_float,int(fields[2]))
                        cdsEnd_float=max(cdsEnd_float,int(fields[7]))
                    exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                    exon_starts = list(map((lambda x: x + txStart ), exon_starts))  #0-based half open [)
                    exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                    exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));           
                    #intron_start = exon_ends[:-1]  #0-based half open [)
                    #intron_end=exon_starts[1:]
                
                    if overlap_flag ==0:
                        Merge_bedNum +=1
                        Merge_bed_TxStart[Merge_bedNum].append(txStart)
                        Merge_bed_TxEnd[Merge_bedNum].append(txEnd)
                        Merge_bed_cdsStart[Merge_bedNum].append(cdsStart)
                        Merge_bed_cdsEnd[Merge_bedNum].append(cdsEnd)                   
                        Merge_bed_ExonStart[Merge_bedNum].extend(exon_starts)
                        Merge_bed_ExonEnd[Merge_bedNum].extend(exon_ends)
                        Merge_bed_geneName[Merge_bedNum].append(geneName)
                        Merge_bed_chr[Merge_bedNum]=chrom
                        Merge_bed_strand[Merge_bedNum]=strand
                    elif overlap_flag ==1:
                        Merge_bed_TxStart[Merge_bedNum].append(txStart)
                        Merge_bed_TxEnd[Merge_bedNum].append(txEnd)
                        Merge_bed_cdsStart[Merge_bedNum].append(cdsStart)
                        Merge_bed_cdsEnd[Merge_bedNum].append(cdsEnd)                   
                        Merge_bed_ExonStart[Merge_bedNum].extend(exon_starts)
                        Merge_bed_ExonEnd[Merge_bedNum].extend(exon_ends)
                        Merge_bed_geneName[Merge_bedNum].append(geneName)
                        Merge_bed_chr[Merge_bedNum]=chrom       
                        Merge_bed_strand[Merge_bedNum]=strand
            self.f.seek(0)
        
        #merge overlapping transcripts on minus strand
            Orig_bedNum=0
            for line in self.f:
                line=line.strip()       
                if bed_line.match(line):
                    fields=line.split()
                    if(len(fields)!=12):
                        print("bed file must be 12 columns speparated by Tab(whilte space)" +line, file=sys.stderr)                 
                    if fields[5] == '+':continue                    
                    Orig_bedNum +=1
                    if Orig_bedNum ==1: #this is first bed entry
                        overlap_flag=0  #has nothing to overlap because it's first entry
                        chrom = fields[0]
                        strand = fields[5]
                        txStart= int(fields[1])
                        txEnd= int(fields[2])
                        cdsStart= int(fields[6])
                        cdsEnd= int(fields[7])      
                        geneName=fields[3]
                        score=fields[4] 
                        txEnd_float=txEnd
                        cdsEnd_float=cdsEnd
                    elif Orig_bedNum > 1:               #this is NOT first entry
                        overlap_flag=1              #we suppose current line is overlapped with previous one. Unless it can prove it's NOT!
                        if fields[0] != chrom:      #not the same chromosome
                            overlap_flag=0
                            txEnd_float = int(fields[2])
                        if (boundary == "utr") and  ( fields[0] == chrom) and (int(fields[1]) >= txEnd_float):
                            overlap_flag=0                                                                                          
                        if ((boundary == "cds") and  ( fields[0] == chrom) and (int(fields[6]) > cdsEdn_float)):
                            overlap_flag=0
                        txStart= int(fields[1])
                        txEnd= int(fields[2])
                        strand = fields[5]
                        chrom = fields[0]               
                        cdsStart= int(fields[6])
                        cdsEnd= int(fields[7])      
                        geneName=fields[3]
                        score=fields[4] 
                        txEnd_float=max(txEnd_float,int(fields[2]))
                        cdsEnd_float=max(cdsEnd_float,int(fields[7]))
                    exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                    exon_starts = list(map((lambda x: x + txStart ), exon_starts))  #0-based half open [)
                    exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                    exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));           
                    #intron_start = exon_ends[:-1]  #0-based half open [)
                    #intron_end=exon_starts[1:]
                
                    if overlap_flag ==0:
                        Merge_bedNum +=1
                        Merge_bed_TxStart[Merge_bedNum].append(txStart)
                        Merge_bed_TxEnd[Merge_bedNum].append(txEnd)
                        Merge_bed_cdsStart[Merge_bedNum].append(cdsStart)
                        Merge_bed_cdsEnd[Merge_bedNum].append(cdsEnd)                   
                        Merge_bed_ExonStart[Merge_bedNum].extend(exon_starts)
                        Merge_bed_ExonEnd[Merge_bedNum].extend(exon_ends)
                        Merge_bed_geneName[Merge_bedNum].append(geneName)
                        Merge_bed_chr[Merge_bedNum]=chrom
                        Merge_bed_strand[Merge_bedNum]=strand
                    elif overlap_flag ==1:
                        Merge_bed_TxStart[Merge_bedNum].append(txStart)
                        Merge_bed_TxEnd[Merge_bedNum].append(txEnd)
                        Merge_bed_cdsStart[Merge_bedNum].append(cdsStart)
                        Merge_bed_cdsEnd[Merge_bedNum].append(cdsEnd)                   
                        Merge_bed_ExonStart[Merge_bedNum].extend(exon_starts)
                        Merge_bed_ExonEnd[Merge_bedNum].extend(exon_ends)
                        Merge_bed_geneName[Merge_bedNum].append(geneName)
                        Merge_bed_chr[Merge_bedNum]=chrom       
                        Merge_bed_strand[Merge_bedNum]=strand
            self.f.seek(0)
        else:       # NOT want to consider strand information. merge + and - together
            Orig_bedNum=0
            Merge_bedNum=0
            for line in self.f:
                line=line.rstrip('\r\n')
                line=line.lstrip()      
                if bed_line.match(line):
                    fields=line.split()     
                    Orig_bedNum +=1
                    if Orig_bedNum ==1: #this is first bed entry
                        overlap_flag=0  #has nothing to overlap because it's first entry
                        chrom = fields[0]
                        #strand = fields[5]
                        txStart= int(fields[1])
                        txEnd= int(fields[2])
                        cdsStart= int(fields[6])
                        cdsEnd= int(fields[7])      
                        geneName=fields[3]
                        score=fields[4] 
                        txEnd_float=txEnd
                        cdsEnd_float=cdsEnd
                    elif Orig_bedNum > 1:               #this is NOT first entry
                        overlap_flag=1              #we suppose current line is overlapped with previous one. Unless it can prove it's NOT!
                        if fields[0] != chrom:      #not the same chromosome
                            overlap_flag=0
                            txEnd_float = int(fields[2])
                        if (boundary == "utr") and  ( fields[0] == chrom) and (int(fields[1]) >= txEnd_float):
                            overlap_flag=0                                                                                          
                        if ((boundary == "cds") and  ( fields[0] == chrom) and (int(fields[6]) > cdsEdn_float)):
                            overlap_flag=0
                        txStart= int(fields[1])
                        txEnd= int(fields[2])
                        #strand = fields[5]
                        chrom = fields[0]               
                        cdsStart= int(fields[6])
                        cdsEnd= int(fields[7])      
                        geneName=fields[3]
                        score=fields[4] 
                        txEnd_float=max(txEnd_float,int(fields[2]))
                        cdsEnd_float=max(cdsEnd_float,int(fields[7]))
                    exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                    exon_starts = list(map((lambda x: x + txStart ), exon_starts))  #0-based half open [)
                    exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                    exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))    
                    #intron_start = exon_ends[:-1]  #0-based half open [)
                    #intron_end=exon_starts[1:]
                
                    if overlap_flag ==0:
                        Merge_bedNum +=1
                        Merge_bed_TxStart[Merge_bedNum].append(txStart)
                        Merge_bed_TxEnd[Merge_bedNum].append(txEnd)
                        Merge_bed_cdsStart[Merge_bedNum].append(cdsStart)
                        Merge_bed_cdsEnd[Merge_bedNum].append(cdsEnd)                   
                        Merge_bed_ExonStart[Merge_bedNum].extend(exon_starts)
                        Merge_bed_ExonEnd[Merge_bedNum].extend(exon_ends)
                        Merge_bed_geneName[Merge_bedNum].append(geneName)
                        Merge_bed_chr[Merge_bedNum]=chrom
                        Merge_bed_strand[Merge_bedNum]='+'
                    elif overlap_flag ==1:
                        Merge_bed_TxStart[Merge_bedNum].append(txStart)
                        Merge_bed_TxEnd[Merge_bedNum].append(txEnd)
                        Merge_bed_cdsStart[Merge_bedNum].append(cdsStart)
                        Merge_bed_cdsEnd[Merge_bedNum].append(cdsEnd)                   
                        Merge_bed_ExonStart[Merge_bedNum].extend(exon_starts)
                        Merge_bed_ExonEnd[Merge_bedNum].extend(exon_ends)
                        Merge_bed_geneName[Merge_bedNum].append(geneName)
                        Merge_bed_chr[Merge_bedNum]=chrom       
                        Merge_bed_strand[Merge_bedNum]='+'      

        
        #meged exon and output
        for id in sorted(Merge_bed_TxStart.keys()):
            final_exon = {}     #key is exon start, value is exon end
            final_exon_sizes=[]
            final_exon_starts=[]
            FName.write("BED_box_" + str(id) + "\t")
            FName.write(','.join(Merge_bed_geneName[id]) + "\n")
            for E_st,E_end in sorted(zip(Merge_bed_ExonStart[id],Merge_bed_ExonEnd[id])):
                startOfFirstExon = E_st
                endOfFirstExon = E_end
                final_exon[E_st] = E_end
                break
            for E_st,E_end in sorted(zip(Merge_bed_ExonStart[id][1:],Merge_bed_ExonEnd[id][1:])):       
                if E_st in final_exon:  #the current start position already there. We only need to compare end position
                    if E_end > endOfFirstExon:
                        final_exon[startOfFirstExon] = E_end
                        endOfFirstExon = E_end
                    else:
                        continue
                else:                       #the current start postion is different
                    if E_st <=endOfFirstExon:
                        if E_end <= endOfFirstExon:
                            continue
                        elif E_end > endOfFirstExon:
                            final_exon[startOfFirstExon] = E_end
                            endOfFirstExon = E_end
                    else:
                        final_exon[E_st]=E_end
                        startOfFirstExon=E_st
                        endOfFirstExon=E_end
                
            for k in sorted(final_exon.keys()):
                final_exon_sizes.append(final_exon[k]-k)
                final_exon_starts.append(k-(min(Merge_bed_TxStart[id])))
            
            FO.write(Merge_bed_chr[id] + "\t" + str(min(Merge_bed_TxStart[id])) + "\t" + str(max(Merge_bed_TxEnd[id])) +"\t")   #column 1,2,3
            FO.write("BED_box_" + str(id) + "\t0\t")                                                                            #column 4,5
            FO.write(Merge_bed_strand[id] +"\t")                                                                                #column 6
            FO.write(str(min(Merge_bed_cdsStart[id])) + "\t" + str(max(Merge_bed_cdsEnd[id])) + "\t" + "255,0,0"+"\t")          #column 7,8,9
            FO.write(str(len(final_exon_starts)) +'\t')                                                                         #column 10
            FO.write(','.join(map(str,final_exon_sizes))+"\t")                                                                  #column 11
            FO.write(','.join(map(str,final_exon_starts)) + "\n")                                                               #column 12

        self.f.seek(0)
        FO.close()      

    def correctSplicingBed(self,genome,outfile,sp="GTAG,GCAG,ATAC"):
        '''input should be bed12 file representing splicing junctions. The function will compare
        the splicing motif to genome. To see if the direcion is correct or not. Only consider GT/AG
        GC/AG, AT/AC motifs. Multiple spliced reads are accepted'''
        
        fout=open(outfile,'w')
        motif=sp.upper().split(',')
        motif_rev = [m.translate(self.transtab)[::-1] for m in motif]
        print("\tloading "+genome+'...', file=sys.stderr)
        tmp=fasta.Fasta(genome)
        for line in self.f:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue
            line=line.strip()
            fields=line.split()
            if (len(fields)<12):
                print(line, file=fout)
                continue
            
            chrom     = fields[0]
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            geneName      = fields[3]
            strand    = fields[5]
            cds_start = int( fields[6] )
            cds_end   = int( fields[7] )                
            exon_num= int(fields[9])
            exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
            exon_starts = list(map(int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
            splice_strand=[]
            for st,end in zip(intron_start,intron_end):
                splice_motif = tmp.fetchSeq(chrom, st, st+2) + tmp.fetchSeq(chrom, end-2, end)
                if splice_motif in motif:
                    splice_strand.append('+')
                elif splice_motif in motif_rev:
                    splice_strand.append('-')
                else:
                    splice_strand.append('.')

            real_strand = set(splice_strand)
            if (len(real_strand) ==1):
                print("\t".join((fields[0],fields[1],fields[2],fields[3],fields[4],real_strand.pop(),fields[6],fields[7],fields[8],fields[9],fields[10],fields[11])), file=fout)
            else:
                print("\t".join((fields[0],fields[1],fields[2],fields[3],fields[4],'.',fields[6],fields[7],fields[8],fields[9],fields[10],fields[11])), file=fout)
        self.f.seek(0)
        fout.close()

    def nrBED(self,outfile=None):
        '''redundant bed entries (exactly the same gene structure) in bed12 file will be merged.'''
        if outfile is None:
            outfile = self.fileName + ".nr.bed"
        FO=open(outfile,'w')
        
        mergeGene=collections.defaultdict(list)
        print("Removing redundcany from " + self.fileName + "  ...", file=sys.stderr)
        for line in self.f:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue  
            fields=line.rstrip().split()
            chrom = fields[0]       #
            txStart= fields[1]  #
            txEnd= fields[2]    #
            geneName=fields[3]
            score=fields[4]
            strand = fields[5]      #
            cdsStart= fields[6] #
            cdsEnd= fields[7]       #
            blockCount = fields[9]  #
            blockSize = fields[10]  #
            blockStart = fields[11] #
            
            key = ":".join((chrom,txStart,txEnd,strand,cdsStart,cdsEnd,blockCount,blockSize,blockStart))
            mergeGene[key].append(geneName)
        for k in mergeGene:
            fields=k.split(":")
            name=";".join(mergeGene[k])
            FO.write(fields[0] +'\t' + fields[1] + '\t' + fields[2] + '\t' + name + '\t0\t' + fields[3] + '\t' + fields[4] +'\t'+fields[5] +'\t0,0,0\t'+fields[6] +'\t' + fields[7] +'\t'+fields[8] +'\n')
        self.f.seek(0)
        FO.close()
    
        
class CompareBED:
    '''Compare two bed fies. Standard BED file has 12 fields. (http://genome.ucsc.edu/FAQ/FAQformat.html)'''
    
    def __init__(self,bedFileA,bedFileB):
        '''This is constructor of ParseBED. Must provide two bed files for comprison. 1st bed file is
        user input bed, while 2nd bed file is usually a reference gene model'''
        self.A_fh=open(bedFileA,'r')
        self.A_full_Name=bedFileA
        self.A_base_Name=os.path.basename(bedFileA)
        self.B_fh=open(bedFileB,'r')
        self.B_full_Name=bedFileB
        self.B_base_Name=os.path.basename(bedFileB)


    def annotateEvents(self,outfile=None):
        '''Compare bed file A to bed file B (usually a bed file for reference gene model). 
        NOTE that only intron boundaries are compared. This function is useful if one want to 
        find if a junctions is novel or not. bed file A will divided into two files: *.known.bed
        and *.novel.bed'''
        
        if outfile is None:
            KnownBed = self.A_base_Name + ".known.bed"
            NovelBed = self.A_base_Name + ".novel.bed"
        else:
            KnownBed = outfile + ".known.bed"
            NovelBed = outfile + ".novel.bed"       
        KNO=open(KnownBed,'w')
        NOV=open(NovelBed,'w')
        
        #print >>sys.stderr,"writing intron to " + outfile + " ..."             
        ref_blocks=collections.defaultdict(dict)
        
        print("reading reference bed file",self.B_full_Name, " ... ", end=' ', file=sys.stderr)
        for line in self.B_fh:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue   
                 # Parse fields from gene tabls
            fields = line.split()
            if(len(fields)<12):
                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                continue
            chrom     = fields[0]
            tx_start = int( fields[1] )
            tx_end   = int( fields[2] )
            if int(fields[9] ==1):
                continue        
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
            for i_st,i_end in zip (intron_start, intron_end):
                key_str=str(i_st)  + '_' + str(i_end)
                ref_blocks[chrom][key_str]=i_end
                
        print("Done", file=sys.stderr)
        
        print("processing",self.A_full_Name, "...", end=' ', file=sys.stderr)
        for line in self.A_fh:
            found=0 
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue   
            # Parse fields from gene tabls
            fields = line.split()
            if(len(fields)<12):
                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                continue
            chrom     = fields[0]
            tx_start = int( fields[1] )
            tx_end   = int( fields[2] )
            if int(fields[9] ==1):
                continue        
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]

            for i_st, i_end in sorted (zip (intron_start,intron_end)):
                key_str=str(i_st)  + '_' + str(i_end)
                if (key_str in ref_blocks[chrom]):
                    #found +=1
            #if (found == len(intron_start)):
                    print(line, end=' ', file=KNO)
                else:
                    print(line, end=' ', file=NOV)
        print("Done", file=sys.stderr)              
                    
    def annotateSplicingSites(self,outfile=None):
        '''Compare bed file A to bed file B (usually a bed file for reference gene model). NOTE that
        only intron boundaries are compared. This function is useful if one want to find if a spilcing
        site is novel or not. bed file A will divided into two files: *.known.bed, *.35novel.bed,
        *.3novel.bed and *.5novel.bed
        
        NOTE: splicing junctions with multiple introns will be split'''
        
        if outfile is None:
            KnownBed = self.A_base_Name + ".known.bed"
            NovelBed5 = self.A_base_Name + ".5novel.bed"        
            NovelBed3 = self.A_base_Name + ".3novel.bed"    
            NovelBed35 = self.A_base_Name + ".35novel.bed"
        else:
            KnownBed = outfile + ".known.bed"
            NovelBed5 = outfile + ".5novel.bed"     
            NovelBed3 = outfile + ".3novel.bed" 
            NovelBed35 = outfile + ".35novel.bed"   
        KNO=open(KnownBed,'w')
        N5=open(NovelBed5,'w')
        N3=open(NovelBed3,'w')
        N35=open(NovelBed35,'w')
        
        #print >>sys.stderr,"writing intron to " + outfile + " ..."             
        refIntronStarts=collections.defaultdict(dict)
        refIntronEnds=collections.defaultdict(dict)
        
        print("\treading reference bed file",self.B_full_Name, " ... ", end=' ', file=sys.stderr)
        for line in self.B_fh:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue   
                 # Parse fields from gene tabls
            fields = line.split()
            if(len(fields)<12):
                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                continue
            chrom     = fields[0]
            tx_start = int( fields[1] )
            tx_end   = int( fields[2] )
            if int(fields[9] ==1):
                continue        
            
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
            for i_st,i_end in zip (intron_start, intron_end):
                refIntronStarts[chrom][i_st] =i_st
                refIntronEnds[chrom][i_end] =i_end
                
        print("Done", file=sys.stderr)
        
        print("\tprocessing",self.A_full_Name, "...", end=' ', file=sys.stderr)
        for line in self.A_fh:
            found=0 
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue   
            # Parse fields from gene tabls
            fields = line.split()
            if(len(fields)<12):
                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                continue
            chrom     = fields[0]
            tx_start = int( fields[1] )
            tx_end   = int( fields[2] )
            geneName = fields[3]
            score = fields[4]
            strand = fields[5]
            if int(fields[9] ==1):
                continue        
            exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end = exon_starts[1:]        
            counter=0
            for i_st, i_end in  zip (intron_start,intron_end):
                counter +=1
                if(strand == '+' or strand == '.'):
                    if (i_st in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                        found =2                                                                        #known both
                    elif (i_st in refIntronStarts[chrom] and i_end not in refIntronEnds[chrom]):
                        found=5                                                                         # 5' splice site known, 3' splice site unkonwn
                    elif (i_st not in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                        found=3                                                                         # 5' splice site uknown, 3' splice site konwn
                    else:
                        found=10
                elif(strand == '-'):
                    if (i_st in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                        found =2                                                                        #known
                    elif (i_st in refIntronStarts[chrom] and i_end not in refIntronEnds[chrom]):
                        found=3                                                                     # 5' splice site uknown, 3' splice site konwn
                    elif (i_st not in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                        found=5                                                                     # 5' splice site known, 3' splice site unkonwn
                    else:
                        found=10
                else:
                    continue
                
                if (found == 2):
                    print('\t'.join((chrom,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),geneName + '_intron' + str(counter),score,strand,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),'0,255,0','2',str(exon_sizes[counter-1]) +','+str(exon_sizes[counter]),'0'+','+str(exon_sizes[counter-1] + i_end - i_st))), file=KNO)
                elif(found ==3):
                    print('\t'.join((chrom,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),geneName + '_intron' + str(counter),score,strand,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),'0,255,0','2',str(exon_sizes[counter-1]) +','+str(exon_sizes[counter]),'0'+','+str(exon_sizes[counter-1] + i_end - i_st))), file=N5)
                elif(found ==5):
                    print('\t'.join((chrom,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),geneName + '_intron' + str(counter),score,strand,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),'0,255,0','2',str(exon_sizes[counter-1]) +','+str(exon_sizes[counter]),'0'+','+str(exon_sizes[counter-1] + i_end - i_st))), file=N3)
                elif(found ==10):
                    print('\t'.join((chrom,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),geneName + '_intron' + str(counter),score,strand,str(i_st - exon_sizes[counter-1]),str(i_end + exon_sizes[counter]),'0,255,0','2',str(exon_sizes[counter-1]) +','+str(exon_sizes[counter]),'0'+','+str(exon_sizes[counter-1] + i_end - i_st))), file=N35)
        print("Done", file=sys.stderr)              
        
    def distribBed(self,outfile=None):
        '''Compare bed file A (usually a bed file of reads mapping results) to bed file B 
        (usually a bed file for reference gene model). For each exon/intron of a gene, calculate
        how many reads mapped into it.'''

        if outfile is None:
            exon_count = self.B_base_Name + "_exon.count.bed"
            intron_count = self.B_base_Name + "_intron.count.bed"
            rscript=self.B_base_Name + ".piechart.r"
            rpdf=self.B_base_Name + ".piechart.pdf"
        else:
            exon_count = outfile + "_exon.count.bed"
            intron_count = outfile +  "_intron.count.bed"
            rscript=outfile + ".piechart.r"
            rpdf=outfile + ".piechart.pdf"
        
        EXON_OUT = open(exon_count,'w')
        INTRON_OUT =open(intron_count,'w')
        R_OUT = open(rscript,'w')
        
        ranges={}
        intronReads=0
        exonReads=0
        intergenicReads=0
        totalReads=0
        splicedReads=0
        
        #read SAM 
        print("reading "+ self.A_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.A_fh:
            if line.startswith("track"):continue
            if line.startswith("#"):continue
            if line.startswith('browser'):continue   
            fields=line.rstrip('\n ').split()
            totalReads +=1
            if(len(fields)==12 and fields[9]>1):
                splicedReads +=1    
                continue
            else:
                chrom=fields[0].upper()
                mid = int(fields[1]) + int((int(fields[2])-int(fields[1]))/2)
                if chrom not in ranges:
                    ranges[chrom] = Intersecter()
                else:
                    ranges[chrom].add_interval( Interval( mid, mid ) )
        self.A_fh.seek(0)
        print("Done", file=sys.stderr)
        
        #read refbed file
        print("Assign reads to "+ self.B_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.B_fh:
            try:
                if line.startswith('#'):continue
                if line.startswith('track'):continue
                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0].upper()
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
                
                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                intron_starts = exon_ends[:-1]
                intron_ends=exon_starts[1:]
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue

                # assign reads to intron                
            if(strand == '-'):
                intronNum=len(intron_starts)
                exonNum=len(exon_starts)
                for st,end in zip(intron_starts,intron_ends):
                    if chrom in ranges:
                        hits= len(ranges[chrom].find(st,end))
                        intronReads += hits
                        INTRON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\n')
                        intronNum -= 1
                            
                for st,end in zip(exon_starts,exon_ends):
                    if chrom in ranges:
                        hits= len(ranges[chrom].find(st,end))
                        exonReads += hits
                        EXON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\n')
                        exonNum -= 1
            elif(strand == '+'):
                intronNum=1
                exonNum=1
                for st,end in zip(intron_starts,intron_ends):
                    if chrom in ranges:
                        hits= len(ranges[chrom].find(st,end))
                        intronReads += hits
                        INTRON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\n')
                        intronNum += 1    
                for st,end in zip(exon_starts,exon_ends):
                    if chrom in ranges:
                        hits= len(ranges[chrom].find(st,end))
                        exonReads += hits
                        EXON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\n')
                        exonNum += 1        
        intergenicReads=totalReads-exonReads-intronReads-splicedReads       
        print("Done." + '\n', file=sys.stderr)
        print("Total reads:\t" + str(totalReads), file=sys.stderr)
        print("Exonic reads:\t" + str(exonReads), file=sys.stderr) 
        print("Intronic reads:\t" + str(intronReads), file=sys.stderr) 
        print("Splicing reads:\t" + str(splicedReads), file=sys.stderr)
        print("Intergenic reads:\t" + str(intergenicReads), file=sys.stderr)
        
        print("writing R script ...", end=' ', file=sys.stderr)
        totalReads=float(totalReads)
        print("pdf('%s')" % rpdf, file=R_OUT)
        print("dat=c(%d,%d,%d,%d)" % (exonReads,splicedReads,intronReads,intergenicReads), file=R_OUT)
        print("lb=c('exon(%.2f)','junction(%.2f)','intron(%.2f)','intergenic(%.2f)')" % (exonReads/totalReads,splicedReads/totalReads,intronReads/totalReads,intergenicReads/totalReads), file=R_OUT)
        print("pie(dat,labels=lb,col=rainbow(4),clockwise=TRUE,main='Total reads = %d')" % int(totalReads), file=R_OUT)
        print("dev.off()", file=R_OUT)
        print("Done.", file=sys.stderr)
        self.B_fh.seek(0)
                
    def distribBedWithStrand(self,outfile=None,output=True):
        '''Compare bed file A (usually a bed file of reads mapping results) to bed file B 
        (usually a bed file for reference gene model). 
        NOTE: When assigning reads (from bed fileA) to gene (bed file B), program will consider
        strand information'''
        
        if output:
            if outfile is None:
                read_count = self.B_base_Name + "_count.xls"
            else:
                read_count = outfile + "_count.xls"
        
            READ_OUT = open(read_count,'w')
        
        Minus_ranges={}
        Plus_ranges={}
        unknown_ranges={}
        redat = collections.defaultdict(list)
        #read bed
        print("\treading "+ self.A_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.A_fh:
            if line.startswith(("track","#","browser")):continue 
            fields=line.rstrip('\n ').split()
            #juntionID = fields[0] + ":" + fields[1] + ":" fields[2] + ":" + fields[5];
            if(len(fields)<6):
                print("[NOTE:input bed must be at least 6 columns] skipped this line: " + line, end=' ', file=sys.stderr)
            if (fields[5] == '-'):
                chrom=fields[0].upper()
                mid = int(fields[1]) + int((int(fields[2])-int(fields[1]))/2)
                if chrom not in Minus_ranges:
                    Minus_ranges[chrom] = IntervalTree()
                else:
                    Minus_ranges[chrom].insert(mid, mid,fields[4])
            elif (fields[5] == '+'):
                chrom=fields[0].upper()
                mid = int(fields[1]) + int((int(fields[2])-int(fields[1]))/2)
                if chrom not in Plus_ranges:
                    Plus_ranges[chrom] = IntervalTree()
                else:
                    Plus_ranges[chrom].insert(mid, mid,fields[4])
            else:
                chrom=fields[0].upper()
                mid = int(fields[1]) + int((int(fields[2])-int(fields[1]))/2)
                if chrom not in unknown_ranges:
                    unknown_ranges[chrom] = IntervalTree()
                else:
                    unknown_ranges[chrom].insert(mid, mid,fields[4])
                
        self.A_fh.seek(0)
        print("Done", file=sys.stderr)
        
        #read refbed file
        print("\tAssign reads to "+ self.B_base_Name + '...', file=sys.stderr)
        for line in self.B_fh:
            try:
                if line.startswith(("track","#","browser")):continue   
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0].upper()
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue

            if output:  
                READ_OUT.write(line.rstrip() + '\t')
                if chrom in Plus_ranges:
                    tmp=Plus_ranges[chrom].find(tx_start,tx_end)
                    if(len(tmp)>0):
                        READ_OUT.write(','.join(tmp) + '\t')
                    else:
                        READ_OUT.write('0\t')
                if chrom in Minus_ranges:
                    tmp=Minus_ranges[chrom].find(tx_start,tx_end)
                    if (len(tmp)>0):
                        READ_OUT.write(','.join(tmp) + '\t')
                    else:
                        READ_OUT.write('0\t')
                        
                if chrom in unknown_ranges:
                    tmp=unknown_ranges[chrom].find(tx_start,tx_end)
                    if (len(tmp)>0):
                        READ_OUT.write(','.join(tmp) + '\n')
                    else:
                        READ_OUT.write('0\n')
            else:
                key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                if chrom in Plus_ranges and len(Plus_ranges[chrom].find(tx_start,tx_end))>0:
                    redat[key].append(','.join(Plus_ranges[chrom].find(tx_start,tx_end)))
                else:
                    redat[key].append('0')
                        
                if chrom in Minus_ranges and len(Minus_ranges[chrom].find(tx_start,tx_end))>0:
                    redat[key].append(','.join(Minus_ranges[chrom].find(tx_start,tx_end)))
                else:
                    redat[key].append('0')
                if chrom in unknown_ranges and len(unknown_ranges[chrom].find(tx_start,tx_end))>0:
                    redat[key].append(','.join(unknown_ranges[chrom].find(tx_start,tx_end)))
                else:
                    redat[key].append('0')
        if output is not True:
            return redat
        print("Done.", file=sys.stderr)
        self.B_fh.seek(0)
                
    def distribSpliceSites(self,outfile=None,output=True):
        '''Compare bed file A (usually a bed file of splicing junctions) to bed file B 
        (usually a bed file for reference gene model). 
        NOTE: When assigning reads (from bed fileA) to gene (bed file B), program will consider
        strand information'''
        
        if output:
            if outfile is None:
                read_count = self.B_base_Name + "_count.xls"
            else:
                read_count = outfile + "_count.xls"
        
            READ_OUT = open(read_count,'w')
        
        Minus_ranges={}
        Plus_ranges={}
        unknown_ranges={}
        redat = collections.defaultdict(list)
        #read bed
        print("\treading "+ self.A_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.A_fh:
            if line.startswith(("track","#","browser")):continue 
            fields=line.rstrip('\n ').split()
            juntionID = fields[0] + ":" + fields[1] + ":" + fields[2] + ":" + fields[5]
            if(len(fields)<12):
                print("[NOTE:input bed must be at least 12 columns] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
            
            
            chrom     = fields[0].upper()
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_starts = exon_ends[:-1]
            intron_ends=exon_starts[1:]
            
            for st,end in zip(intron_starts,intron_ends):
                junctionID=chrom + ":" + str(st) + ":" + str(end)
                if (fields[5] == '-'):
                    if chrom not in Minus_ranges:
                        Minus_ranges[chrom] = IntervalTree()
                    else:
                        Minus_ranges[chrom].insert(st, st,(junctionID,fields[4]))
                        Minus_ranges[chrom].insert(end, end,(junctionID,fields[4]))
                elif (fields[5] == '+'):
                    if chrom not in Plus_ranges:
                        Plus_ranges[chrom] = IntervalTree()
                    else:
                        Plus_ranges[chrom].insert(st, st,(junctionID,fields[4]))
                        Plus_ranges[chrom].insert(end, end,(junctionID,fields[4]))
                else:
                    if chrom not in unknown_ranges:
                        unknown_ranges[chrom] = IntervalTree()
                    else:
                        unknown_ranges[chrom].insert(st, st,(junctionID,fields[4]))
                        unknown_ranges[chrom].insert(end, end,(junctionID,fields[4]))
                
        self.A_fh.seek(0)
        print("Done", file=sys.stderr)
        
        #read refbed file
        print("\tAssign reads to "+ self.B_base_Name + '...', file=sys.stderr)
        for line in self.B_fh:
            try:
                if line.startswith(("track","#","browser")):continue   
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0].upper()
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue

            if output:  
                READ_OUT.write(line.rstrip() + '\t')
                if chrom in Plus_ranges:
                    tmp=Plus_ranges[chrom].find(tx_start,tx_end)
                    locus=[]
                    junction_list=[]
                    if(len(tmp)>0):
                        for j in tmp:
                            if j[0] not in junction_list:
                                locus.append(j[1])
                                junction_list.append(j[0])
                        
                        READ_OUT.write(','.join(locus) + '\t')
                    else:READ_OUT.write('0\t')
                else:READ_OUT.write('0\t')
                
                if chrom in Minus_ranges:
                    tmp=Minus_ranges[chrom].find(tx_start,tx_end)
                    locus=[]
                    junction_list=[]
                    if(len(tmp)>0):
                        for j in tmp:
                            if j[0] not in junction_list:
                                locus.append(j[1])
                                junction_list.append(j[0])
                        
                        READ_OUT.write(','.join(locus) + '\t')
                    else:READ_OUT.write('0\t')
                else:READ_OUT.write('0\t')
                
                if chrom in unknown_ranges:
                    tmp=unknown_ranges[chrom].find(tx_start,tx_end)
                    locus=[]
                    junction_list=[]
                    if(len(tmp)>0):
                        for j in tmp:
                            if j[0] not in junction_list:
                                locus.append(j[1])
                                junction_list.append(j[0])
                                
                        READ_OUT.write(','.join(locus) + '\n')
                    else:READ_OUT.write('0\n')
                else:READ_OUT.write('0\n')

            else:
                key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                if chrom in Plus_ranges:
                    tmp = Plus_ranges[chrom].find(tx_start,tx_end)
                    locus=[]
                    junction_list=[]
                    if(len(tmp)>0):
                        for j in tmp:
                            if j[0] not in junction_list:
                                locus.append(j[1])
                                junction_list.append(j[0])
                    
                        redat[key].append(','.join(locus))
                    else:redat[key].append('0')
                else:redat[key].append('0')
                        
                        
                if chrom in Minus_ranges:
                    tmp = Minus_ranges[chrom].find(tx_start,tx_end)
                    locus=[]
                    junction_list=[]
                    if(len(tmp)>0):
                        for j in tmp:
                            if j[0] not in junction_list:
                                locus.append(j[1])
                                junction_list.append(j[0])
                    
                        redat[key].append(','.join(locus))
                    else:redat[key].append('0')
                else:redat[key].append('0')
                
                
                if chrom in unknown_ranges:
                    tmp = unknown_ranges[chrom].find(tx_start,tx_end)
                    locus=[]
                    junction_list=[]
                    if(len(tmp)>0):
                        for j in tmp:
                            if j[0] not in junction_list:
                                locus.append(j[1])
                                junction_list.append(j[0])
                    
                        redat[key].append(','.join(locus))
                    else:redat[key].append('0')
                else:redat[key].append('0')

        if output is not True:
            return redat
        print("Done.", file=sys.stderr)
        self.B_fh.seek(0)
        
    def findClosestTSS(self,outfile=None,downStream=50000,upStream=50000):
        '''For each entry in input bed file (1st bed file), find the nearest gene (2nd bed file) based
        on TSS. Genes shared the same TSS will be grouped together.
        NOTE: gene is represented by its TSS, input bed entry is represented by its middle point.
        NOTE: this is peak centered, for each peak find the nearest gene(s)'''
        
        if outfile is None:
            outfileName = self.A_base_Name + ".nearestTSS.xls"
        else:
            outfileName = outfile + ".nearestTSS.xls"
        OUT=open(outfileName,'w')
        
        ranges={}
        tss_group={}
        tss_group_num=0
        #read reference bed file
        print("Reading "+ self.B_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.B_fh:
            if line.startswith(("#","track","browser")):continue
            fields=line.rstrip('\n ').split()
            
            chrom = fields[0]
            geneName = fields[3]
            if (len(fields) >=6):
                if (fields[5] == '-'):
                    tss_st = int(fields[2]) - 1
                    tss_end = int(fields[2])
                else:
                    tss_st = int(fields[1])
                    tss_end =int(fields[1]) + 1
            elif (len(fields) < 6 and len(fields) >=3 ):
                tss_st = int(fields[1])
                tss_end =int(fields[1]) + 1
            else:
                print("reference bed file must be at least 3 columns", file=sys.stderr)
                sys.exit()
            
            key=chrom + ":" + str(tss_st) + ':' + str(tss_end)
            if key not in tss_group:
                tss_group[key] = geneName + ";"
            else:
                tss_group[key] += geneName + ";"
            
        for key in tss_group:
            tss_group_num +=1
            chrom,tss_st,tss_end = key.split(":")
            #print chrom + '\t'+ tss_st +'\t' + tss_end
            if chrom not in ranges:
                ranges[chrom] = IntervalTree()
            ranges[chrom].insert_interval( Interval( int(tss_st), int(tss_end),value=tss_group[key] + '\t' + key) )
            
        self.B_fh.seek(0)
        print("Done. Total " + str(tss_group_num) + " TSS groups", file=sys.stderr)

        #a=ranges['chr1'].find(1000000,1200000)
        #print a
        
        
        
        #read input bed file
        print("Find nearest TSS(s) for "+ self.A_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.A_fh:
            if line.startswith(("#","track","browser")):continue
            fields=line.rstrip('\n ').split()

            if (len(fields)>=6):
                chain = fields[5]
            elif (len(fields) < 6 and len(fields) >=3 ):
                chain = '+'
            else:
                print("Input bed file must be at least 3 columns", file=sys.stderr)
                sys.exit()  
                
            chrom = fields[0]
            bed_st =  int(fields[1]) + int((int(fields[2]) - int(fields[1]))/2)
            bed_end = bed_st +1
            
            up=ranges[chrom].upstream_of_interval(Interval(bed_st,bed_end,strand=chain),num_intervals=1,max_dist=upStream)
            down=ranges[chrom].downstream_of_interval(Interval(bed_st,bed_end,strand=chain),num_intervals=1,max_dist=downStream)
            
            if (len(up) >0):
                up_name=up[0].value
                up_dist=abs(up[0].end - bed_end)
            else:
                up_name="NA\tNA"
                up_dist="NA"
            if (len(down) >0):
                down_name=down[0].value
                down_dist=abs(down[0].end - bed_end)
            else:
                down_name="NA\tNA"
                down_dist="NA"          
            print(line.rstrip() + "\t" + up_name + "\t" + str(up_dist) + '\t' + down_name + "\t" + str(down_dist), file=OUT)
            
        self.A_fh.seek(0)
        print("Done.", file=sys.stderr)


    def findClosestTTS(self,outfile=None,downStream=50000,upStream=50000):
        '''For each entry in input bed file (1st bed file), find the nearest gene (2nd bed file) based
        on TTS. Genes shared the same TTS will be grouped together.
        NOTE: gene is represented by its TTS, input bed entry is represented by its middle point.'''
        
        if outfile is None:
            outfileName = self.A_base_Name + ".nearestTTS.xls"
        else:
            outfileName = outfile + ".nearestTTS.xls"
        OUT=open(outfileName,'w')
        
        ranges={}
        tts_group={}
        tts_group_num=0
        #read reference bed file
        print("Reading "+ self.B_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.B_fh:
            if line.startswith(("#","track","browser")):continue
            fields=line.rstrip('\n ').split()
            
            chrom = fields[0]
            geneName = fields[3]
            if (len(fields) >=6):
                if (fields[5] == '-'):
                    tts_st = int(fields[1])
                    tts_end = int(fields[1]) + 1
                else:
                    tts_st = int(fields[2]) - 1
                    tts_end =int(fields[2])
            elif (len(fields) < 6 and len(fields) >=3 ):
                tts_st = int(fields[2]) - 1
                tts_end =int(fields[2])
            else:
                print("reference bed file must be at least 3 columns", file=sys.stderr)
                sys.exit()
            
            key=chrom + ":" + str(tts_st) + ':' + str(tts_end)
            if key not in tts_group:
                tts_group[key] = geneName + ";"
            else:
                tts_group[key] += geneName + ";"
            
        for key in tts_group:
            tts_group_num +=1
            chrom,tts_st,tts_end = key.split(":")
            #print chrom + '\t'+ tts_st +'\t' + tts_end
            if chrom not in ranges:
                ranges[chrom] = IntervalTree()
            ranges[chrom].insert_interval( Interval( int(tts_st), int(tts_end),value=tts_group[key] + '\t' + key) )
            
        self.B_fh.seek(0)
        print("Done. Total " + str(tts_group_num) + " TTS groups", file=sys.stderr)

        #a=ranges['chr1'].find(1000000,1200000)
        #print a
        
        
        
        #read input bed file
        print("Find nearest TTS(s) for "+ self.A_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.A_fh:
            if line.startswith(("#","track","browser")):continue
            fields=line.rstrip('\n ').split()

            if (len(fields)>=6):
                chain = fields[5]
            elif (len(fields) < 6 and len(fields) >=3 ):
                chain = '+'
            else:
                print("Inut bed file must be at least 3 columns", file=sys.stderr)
                sys.exit()  
                
            chrom = fields[0]
            bed_st =  int(fields[1]) + int((int(fields[2]) - int(fields[1]))/2)
            bed_end = bed_st +1
            
            up=ranges[chrom].upstream_of_interval(Interval(bed_st,bed_end,strand=chain),num_intervals=1,max_dist=upStream)
            down=ranges[chrom].downstream_of_interval(Interval(bed_st,bed_end,strand=chain),num_intervals=1,max_dist=downStream)
            
            if (len(up) >0):
                up_name=up[0].value
                up_dist=abs(up[0].end - bed_end)
            else:
                up_name="NA\tNA"
                up_dist="NA"
            if (len(down) >0):
                down_name=down[0].value
                down_dist=abs(down[0].end - bed_end)
            else:
                down_name="NA\tNA"
                down_dist="NA"          
            print(line.rstrip() + "\t" + up_name + "\t" + str(up_dist) + '\t' + down_name + "\t" + str(down_dist), file=OUT)
            
        self.A_fh.seek(0)
        print("Done.", file=sys.stderr)     


    def findClosestPeak(self,mod, downStream=50000,upStream=50000):
        '''For each entry in second bed file (reference gene model) find the closest peak defined in the fist bed file'''
        
        mode={
        0:'TSS-up, TSS-down',
        1:'TSS-up, TES-down',
        2:'TES-up, TES-down',
        3:'CDSS-up, CDSS-down',
        4:'CDSS-up, CDSE-down',
        5:'CDSE-up, CDSE-down'
        }
        ranges={}
        
        #read peak bed file
        print("\tReading "+ self.A_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.A_fh:
            if line.startswith(("#","track","browser")):continue
            if not line.strip(): continue
            line = line.strip('\n')
            fields=line.split()
            
            if len(fields)>=3 and len(fields) < 6:
                chrom = fields[0]
                start = fields[1]
                end = fields[2]
                peak_id = "ID=" + chrom + ":" + start + '-' + end + ":" + 'Score=0'
                #middle = int(int(start) + (int(start) + int(end))/2)
            if len(fields)>=6:
                chrom = fields[0]
                start = fields[1]
                end = fields[2]
                peak_id = "ID=" + chrom + ":" + start + '-' + end + ":" + 'Score=' + fields[4]
                #middle = int(int(start) + (int(start) + int(end))/2)
            if chrom not in ranges:
                ranges[chrom] = IntervalTree()
            ranges[chrom].insert_interval( Interval(int(start), int(end), value = peak_id) )
        print("Done", file=sys.stderr)

        #read reference gene model bed file
        line_id=0
        nearest_peak={}
        print("\tFind nearest peaks for "+ self.B_base_Name + '...', end=' ', file=sys.stderr)
        for line in self.B_fh:
            if line.startswith(("#","track","browser")):continue
            if not line.strip():continue
            fields=line.rstrip('\n ').split()
            line_id +=1
            try:
                strand = fields[5]      
                chrom = fields[0]
                if strand == '-':
                    TSS = int(fields[2])    
                    TES = int(fields[1])    
                    CDSS = int(fields[7])   
                    CDSE = int(fields[6])
                elif strand != '-':
                    strand='+'
                    TSS = int(fields[1])    
                    TES = int(fields[2])    
                    CDSS = int(fields[6])   
                    CDSE = int(fields[7])   
                geneName=fields[3]
                score=fields[4]
                blockCount = fields[9]  
                blockSize = fields[10]  
                blockStart = fields[11] 
                #line_id = geneName
            except:
                print("Reference gene model must 12 column BED files", file=sys.stderr)
                sys.exit(1)
            if chrom in ranges:
                if mod==0:  #TSS-up, TSS-down
                    hits = ranges[chrom].find(TSS - upStream, TSS + downStream )
                elif mod==1:    # TSS-up, TES-down
                    if strand == '+': hits = ranges[chrom].find(TSS - upStream, TES + downStream)
                    if strand == '-': hits = ranges[chrom].find(TSS + upStream, TES - downStream)
                elif mod==2:    #TES-up, TES-down
                    hits = ranges[chrom].find(TES - upStream, TES + downStream )
                elif mod==3:    #CDSS-up, CDSS-down
                    hits = ranges[chrom].find(CDSS - upStream, CDSS + downStream )
                elif mod==4:    #CDSS-up, CDSE-down
                    if strand == '+': hits = ranges[chrom].find(CDSS - upStream, CDSE + downStream)
                    if strand == '-': hits = ranges[chrom].find(CDSS + upStream, CDSE - downStream)
                elif mod==5:    #CDSE-up, CDSE-down
                    hits = ranges[chrom].find(CDSE - upStream, CDSE + downStream )
                else:
                    print("unknow arguments for 'mod'", file=sys.stderr)
                    return None
                if (len(hits) >0):
                    val = ','.join([i.value for i in hits])
                    nearest_peak[line_id]=val
                else:
                    nearest_peak[line_id] ="NA"
            else: nearest_peak[line_id] = "NA"
        print("Done.", file=sys.stderr)
        return nearest_peak 
        #FO=open("aaa",'w')
        #for k,v in nearest_peak.items():
        #   print >>FO,  str(k) + '\t' + str(v)
                        

    def bestMatch(self):
        '''Exon chain comparison. Given a bed entry in bedFileA, find the best-matched gene from bedFileB.
        If multiple genes from bedFileB matched equally good. Randomly report one.'''
        
        #read reference gene model

        ref_ranges={}
        ref_cdsRange={}
        print("\tReading reference gene model " + self.B_base_Name + '...', file=sys.stderr)
        for line in self.B_fh:
            boundaries=set()
            #try:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue
            fields=line.rstrip('\r\n').split()
            chrom = fields[0]
            txStart = int(fields[1])
            txEnd = fields[2]
            geneName = fields[3]
            score = fields[4]
            strand = fields[5]
            cdsStart = fields[6]
            cdsEnd = fields[7]
            geneID = chrom + ":" + str(txStart) + "-" + txEnd + ":" + strand + ':' + geneName
            
            #boundaries.add(int(cdsStart))
            #boundaries.add(int(cdsEnd))
            
            ref_cdsRange[geneID] = [int(cdsStart),int(cdsEnd)]
            
            exon_start=list(map(int,fields[11].rstrip(',').split(',')))
            exon_start=list(map((lambda x: x + txStart),exon_start))
            exon_end=list(map(int,fields[10].rstrip(',').split(',')))
            exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
            #except:
            #   print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
            #   continue
            for st,end in zip(exon_start,exon_end):
                boundaries.add(st)
                boundaries.add(end)
            chrom = chrom + ':' + strand
            if chrom not in ref_ranges:
                ref_ranges[chrom] = IntervalTree()
            ref_ranges[chrom].insert(int(txStart), int(txEnd),{geneID : boundaries})
        
        print("\tReading " + self.A_base_Name + '...', file=sys.stderr)
        #read bed file
        utr_diff_score = -0.5
        cds_diff_score = -1
        common_score = 8
        line_id=0
        ret_dict={}
        for line in self.A_fh:
            boundaries=set()
            try:
                if line.startswith('#'):continue
                if line.startswith('track'):continue
                if line.startswith('browser'):continue
                if not line.strip():continue
                
                line_id += 1
                fields=line.rstrip('\r\n').split()
                chrom = fields[0]
                txStart = int(fields[1])
                txEnd = fields[2]
                geneName = fields[3]
                score = fields[4]
                strand = fields[5]
                cdsStart = fields[6]
                cdsEnd = fields[7]
                geneID = chrom + ":" + str(txStart) + "-" + txEnd + ":" + strand + ':' + geneName
                
                #boundaries.add(int(cdsStart))
                #boundaries.add(int(cdsEnd))
                
                exon_start=list(map(int,fields[11].rstrip(',').split(',')))
                exon_start=list(map((lambda x: x + txStart),exon_start))
                exon_end=list(map(int,fields[10].rstrip(',').split(',')))
                exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue

            for st,end in zip(exon_start,exon_end):
                boundaries.add(st)
                boundaries.add(end)

            chrom = chrom + ':' + strand
            #print line.strip() + '\t',
            if chrom in ref_ranges:
                overlap_genes = ref_ranges[chrom].find(int(txStart),int(txEnd))
                if len(overlap_genes)==0:               #not overlap with any gene
                    status = "novel_gene(non-overlap with known gene)"
                else:                                   #overlap with known gene
                    diff_dict=collections.defaultdict(int)      #note the higher score, the better similarity

                    #which ref gene is best match
                    for i in overlap_genes:
                        status = ''
                        ref_gene_marks =   list(i.values())[0]
                        ref_gene_id = list(i.keys())[0]
                        input_gene_marks = boundaries
                        
                        ref_gene_uniq = ref_gene_marks.difference(input_gene_marks)
                        input_gene_uniq = input_gene_marks.difference(ref_gene_marks)
                        common = input_gene_marks.intersection(ref_gene_marks)
                        if ref_gene_marks.__eq__(input_gene_marks):     #boundaries 100% match to known exon chain
                            status = 'complete_match' + '_' + ref_gene_id
                            break                           
                        
                        for bd in common:
                            diff_dict[ref_gene_id] += common_score
                        for bd in input_gene_uniq:
                            if bd >= int(cdsStart)  or bd <= int(cdsEnd):
                                diff_dict[ref_gene_id] += cds_diff_score                
                            else:
                                diff_dict[ref_gene_id] += utr_diff_score
                                
                        for bd in ref_gene_uniq:
                            if bd >= ref_cdsRange[ref_gene_id][0] or bd <= ref_cdsRange[ref_gene_id][1]:
                                diff_dict[ref_gene_id] += cds_diff_score
                            else:
                                diff_dict[ref_gene_id] += utr_diff_score
                    
                    if status.find('complete_match')==-1:
                        bestID = max(diff_dict,key=diff_dict.get)
                        for i in overlap_genes:
                            if list(i.keys())[0] == bestID:
                                ref_best_match_set = list(i.values())[0]
                                break

                        #we found a ref gene best match to input gene
                        input_gene_marks = boundaries
                        i_union = ref_best_match_set.union(input_gene_marks)
                        i_intersection = ref_best_match_set.intersection(input_gene_marks)  
                        ref_gene_uniq = ref_best_match_set.difference(input_gene_marks)
                        input_gene_uniq = input_gene_marks.difference(ref_best_match_set)           
                    
                        if len(i_intersection) == 0:
                            status = "novel_gene(overlap with known gene)" + bestID
                        else:
                            tmp=set()
                            for bd in ref_gene_uniq:
                                if bd >= ref_cdsRange[bestID][0] or bd <= ref_cdsRange[bestID][1]:
                                    tmp.add('different CDS')
                                    #status = "partial_match (different CDS)" + bestID
                                else:
                                    tmp.add('different UTR')
                                    #status = "partial_match (different UTR)" + bestID
                            status = 'partial_match (' + ';'.join(tmp) + ')'
            else:
                status = 'unknown'
            ret_dict[line_id] = status
        return ret_dict

def unionBed3(lst):
    '''Take the union of 3 column bed files. return a new list'''
    bitsets = binned_bitsets_from_list(lst)
    ret_lst=[]
    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while 1:
            start = bits.next_set( end )
            if start == bits.size: break
            end = bits.next_clear( start )
            ret_lst.append([chrom, start, end])
    bitsets=dict()
    return ret_lst

def intersectBed3(lst1,lst2):
    '''Take the intersection of two bed files (3 column bed files)'''
    bits1 = binned_bitsets_from_list(lst1)
    bits2 = binned_bitsets_from_list(lst2)

    bitsets = dict()
    ret_lst = []
    for key in bits1:
        if key in bits2:
            bits1[key].iand( bits2[key] )
            bitsets[key] = bits1[key]

    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while 1:
            start = bits.next_set( end )
            if start == bits.size: break
            end = bits.next_clear( start )
            ret_lst.append([chrom, start, end])
    bits1.clear()
    bits2.clear()
    bitsets.clear()
    return ret_lst

def subtractBed3(lst1,lst2):
    '''subtrack lst2 from lst1'''
    bitsets1 = binned_bitsets_from_list(lst1)
    bitsets2 = binned_bitsets_from_list(lst2)
    
    ret_lst=[]
    for chrom in bitsets1:  
        if chrom not in bitsets1:
            continue
        bits1 = bitsets1[chrom]
        if chrom in bitsets2:
            bits2 = bitsets2[chrom]
            bits2.invert()
            bits1.iand( bits2 )
        end=0
        while 1:
            start = bits1.next_set( end )
            if start == bits1.size: break
            end = bits1.next_clear( start )
            ret_lst.append([chrom,start,end])
    bitsets1 = dict()
    bitsets2 = dict()
    return ret_lst

def tillingBed(chrName,chrSize,stepSize=10000):
    '''tilling whome genome into small sizes'''
    #tilling genome
    for start in range(0,chrSize,stepSize):
        end = start + stepSize
        if end < chrSize:
            yield (chrName,start,end)
        else:
            yield (chrName,start,chrSize)
        
