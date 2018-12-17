#!/usr/bin/env python
'''
manipulate fasta for fastq format files.
This python3 module was converted from python2.7 code using 2to3
'''

#import built-in modules
import numpy
import re
import sys
from optparse import OptionParser
import collections
#import third-party modules

#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="3.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


class Fasta:
    '''manipulate fasta or fastq format file
    '''
    
    def __init__(self,fastafile=None):
        '''initialize object, lowercase in sequence is automatically converted into uppercase'''
        self.seqs={}
        self.IDs=[]
        self.transtab = str.maketrans('ACGTNX', 'TGCANX')
        self.filename = fastafile
        tmpseq=''
        if fastafile is not None:
            for line in open(fastafile,'r'):
                line=line.strip(' \n')
                if line.startswith('>'):
                    if(tmpseq):
                        self.seqs[name]=tmpseq
                    name=line[1:]
                    tmpseq =''
                    self.IDs.append(name)
                    print("\tloading "+name+' ...', file=sys.stderr)
                else:
                    tmpseq += line.upper()
            self.seqs[name]=tmpseq
                
    def addSeq(self,id,seq):
        '''add sequence to current data'''
        if id in self.seqs:
            print(id +" already exists!", file=sys.stderr)
            return
        else:
            self.seqs[id]=seq.upper()
            self.IDs.append(id)
            
    def getNames (self,file=None):
        '''return all sequence IDs'''
        return self.IDs
        
    def getSeq(self,seqID=None):
        '''return sequence for sepcified seqID, otherwise all sequences are returned'''
        if seqID is None:
            return list(self.seqs.values())
        else:
            return self.seqs[seqID]

    def printSeqs(self,n=50):
        '''print all seqs '''
        for k,v in list(self.seqs.items()):
            print('>' + k)
            for i in range(0, len(v), n):
                print(v[i:i+n])

            
    def getSeqLen(self,seqID=None):
        seqlen=collections.defaultdict(dict)
        if seqID is None:
            for (k,v) in list(self.seqs.items()):
                seqlen[k]=len(v)
        else:
            try:
                seqlen[seqID]=len(self.seqs[seqID])
            except:
                print("Not found", file=sys.stderr)
        return seqlen
            
    def countBase(self,pattern=None):
        '''count occurence of substring (defined by pattern), otherwise count A,C,G,T,N,X
        NOTE: pattern is counted non-overlappingly'''
        if pattern is None:
            print("ID\tTotal\tA\tC\tG\tT\tN\tX")
            for (k,v) in list(self.seqs.items()):
                print(k+"\t", end=' ')
                print(len(v),"\t", end=' ')
                print(str(v.count('A'))+"\t", end=' ')
                print(str(v.count('C'))+"\t", end=' ')
                print(str(v.count('G'))+"\t", end=' ')
                print(str(v.count('T'))+"\t", end=' ')
                print(str(v.count('N'))+"\t", end=' ')
                print(v.count('X'))
        else:
            for (k,v) in list(self.seqs.items()):
                print(k+"\t", end=' ')
                print(str(len(v))+"\t", end=' ')
                print(v.count(pattern))
                
    def revComp(self,seqID=None):
        '''return reverse-complemented sequence for sepcified seqID, otherwise all sequences are 
        reverse-complemented'''
        if seqID is None:
            for (k,v) in list(self.seqs.items()):
                print(">" + k + "_rev")
                tmp = v.translate(self.transtab)
                return tmp[::-1]            
        else:
            return self.seqs[seqID].translate(self.transtab)[::-1]

            
    def getUniqSeqs(self):
        '''remove redundancy from original fasta files.
        duplicated sequences will be only report once'''

        seq2Name={}
        seq2Count={}
        for (key,value) in list(self.seqs.items()):
            seq2Name[value]=key
            if value in seq2Count:
                seq2Count[value]+=1
            else:
                seq2Count[value]=1
        for value in list(seq2Name.keys()):
                print('>'+ str(seq2Name[value]) + '_' + str(seq2Count[value]))
                print(value)


    def findPattern(self,pat,outfile,seqID=None,rev=True):
        ''' find pattern in all sequence unless seqID is specified, coordinates will be returned as bed format file'''
        
        fout=open(outfile,'w')
        length=len(pat) 

        Pat=pat.upper()
        start=0
        
        
        if seqID is None:
            for (k,v) in list(self.seqs.items()):
                loopSwitch=0
                start=0
                while loopSwitch !=-1:
                    loopSwitch = v.find(Pat,start)
                    if loopSwitch !=-1:
                        print(k + "\t" + str(loopSwitch) + "\t" + str(loopSwitch + length) + "\t" + Pat + "\t0\t+", file=fout) 
                        start = loopSwitch +1
                    
        else:
            loopSwitch=0
            start=0
            while loopSwitch !=-1:
                loopSwitch = self.seqs[seqID].find(Pat,start)
                print(seqID + "\t" + str(loopSwitch) + "\t" + str(loopSwitch + length) + "\t" + Pat + "\t0\t+", file=fout) 
                start = loopSwitch +1

        if rev==True:
            Pat_rev=Pat.translate(self.transtab)[::-1]
            if seqID is None:
                for (k,v) in list(self.seqs.items()):
                    loopSwitch=0
                    start=0
                    while loopSwitch !=-1:
                        loopSwitch = v.find(Pat_rev,start)
                        if loopSwitch !=-1:
                            print(k + "\t" + str(loopSwitch) + "\t" + str(loopSwitch + length) + "\t" + Pat + "\t0\t-", file=fout) 
                            start = loopSwitch +1
                        
            else:
                loopSwitch=0
                start=0
                while loopSwitch !=-1:
                    loopSwitch = self.seqs[seqID].find(Pat_rev,start)
                    print(seqID + "\t" + str(loopSwitch) + "\t" + str(loopSwitch + length) + "\t" + Pat + "\t0\t-", file=fout) 
                    start = loopSwitch +1       

    def fetchSeq(self,chr=None,st=None,end=None,infile=None,outfile=None):
        ''' Fetching sequence based on chrName (should be exactly the same as fasta file), St, End. 
        NOTE: the coordinate is 0-based,half-open. use infile to specify multiple coordinates. infile
        should be bed3, bed6 or bed12'''

        if (infile is not None) and (outfile is not None):
            fout=open(outfile,'w')
            for line in open(infile):
                fields=line.strip().split()
                if (len(fields)==3):
                    print(fields[0]+":"+fields[1]+"-"+fields[2]+"\t"+"strand=+", file=fout)
                    print(self.seqs[fields[0]][int(fields[1]):int(fields[2])].upper(), file=fout)
                elif (len(fields)>3):
                    if fields[5]=='-':
                        print(fields[0]+":"+fields[1]+"-"+fields[2]+"\t"+"strand=-", file=fout)
                        print(self.seqs[fields[0]][int(fields[1]):int(fields[2])].translate(self.transtab)[::-1].upper(), file=fout)
                    else:
                        print(fields[0]+":"+fields[1]+"-"+fields[2]+"\t"+"strand=+", file=fout)
                        print(self.seqs[fields[0]][int(fields[1]):int(fields[2])].upper(), file=fout)
        else:
            try:
                return self.seqs[chr][st:end].upper()
            except:
                print("cannot fetch sequence from " + self.filename + " for " + chr + ":" + str(st) + "-" + str(end), file=sys.stderr)
                return ''
                #print >>sys.stderr, chr + "\t" + str(st) +'\t' + str(end) + "  Please input chr,st,end"
                        
def main():
    parser = OptionParser()
    parser.add_option("-i","--input_file",dest="in_file",help="input file name")
    (options,args)=parser.parse_args()
    obj=Fasta(options.in_file)
    obj.printSeqs(n=80)

if __name__ == "__main__":
    main()
    
