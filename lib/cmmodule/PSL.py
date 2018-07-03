#!/usr/bin/env python
'''manipulate blat PSL file.'''

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings

#import third-party modules

#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2010, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production


class PSL:
    '''manipulate PSL format file (blat output file)'''
    
    def __init__(self,inputfile_name,score=20,block=2):
        '''initialize this class
        arg1: inputfile psl file
        arg2: score cutoff. only report alignments with score >= this value. default=20
        arg3: block count cutoff. only report alignments with block number <= this value. default=2
        eg: a=PSL("filename",30,2) or a=PSL.PSL("filename",30,2)
        '''
        self.__inputfile=open(inputfile_name,'r')
        self.__scoreCutoff=score        #matched score smaller than this value will be removed
        self.__blockCutoff=block        #matched blcok bigger than this value will be removed
        self.__pslLine=re.compile(r'^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+[+-]')   #use this re to remove head lines
        self.__pslSplit=re.compile(r'\s+')
        self.__blankLine=re.compile(r'\s*\n')
        self.__data=[]
        totalLine=0
        usedLine=0
        nonPslLine=0
        field=[]
        while True:
            self.__line=self.__inputfile.readline()
            if self.__blankLine.match(self.__line):     #skip blank line
                continue
            elif self.__pslLine.match(self.__line):
                totalLine=totalLine+1
                field=self.__pslSplit.split(self.__line)
                if string.atoi(field[0])< self.__scoreCutoff:
                    continue
                elif string.atoi(field[17]) > self.__blockCutoff:
                    continue
                usedLine = usedLine+1
                self.__data.append(self.__line)
                
            else:
                nonPslLine = nonPslLine+1
            if not self.__line:break    #end of file
            
        print  "\nTotal: ", totalLine, "lines"
        print  "Used: ", usedLine, "lines"
        print  "Non-PSL: ", nonPslLine,"lines","\n"
          
    def head(self,limit=10):
        '''print out header lines of PSL file, default first 10 lines
        eg: a.head(50)
        '''
        count=0 #count how many lines have been printed
        field=[]
        for line in self.__data:
            line=line.rstrip("\n")
            print line
            count=count+1
            if count >= limit:break
            
    def psl2bedFile(self,output_file=None):
        '''transform psl format into bed format. creat 6 column bed files.
        col1: target name, typically a chromosome
        col2: start coordinate
        col3: end coordinate
        col4: name of query. if a query is split. sequential number will be added
             at the end of the name
        col5: score. "matchScore.mismatch"
        col6: strand
        '''
        
        if output_file is None:
        	sys.stdout=sys.__stdout__
        else:
        	outfile=open(output_file,'w')
        	sys.stdout=outfile
        field=[]
        blockSize=[]
        blockStart=[]
        for line in self.__data:
            line=line.rstrip('\n')
            field=self.__pslSplit.split(line)
            if string.atoi(field[17]) == 1:
                print field[13],"\t",field[15],"\t",field[16],"\t",field[9],"\t",field[0]+'.'+field[1],"\t",field[8]
            else:
                blockSize=field[18].split(',')
                blockSize.pop(-1)
                blockStart=field[20].split(',')
                blockStart.pop(-1)
                for i in range(0,len(blockSize)):
                    print field[13],"\t",blockStart[i],"\t",string.atoi(blockStart[i])+string.atoi(blockSize[i]),"\t",field[9]+'.'+str(i+1),"\t",field[0]+'.'+field[1],"\t",field[8]

def main():
	parser=OptionParser()
	parser.add_option('-i','--input_file',dest="inputFileName",help="Input file name")
	parser.add_option('-o','--output_file',dest="outputFileName",help="Output file name")
	parser.add_option('-s','--score',dest="scoreCutoff",help="blat mapping score cutoff",type="int",default=20)
	parser.add_option('-b','--block_num',dest="blockCutoff",help="block number cutoff",type="int",default=2)
	(options,args)=parser.parse_args()
	
	obj=PSL(options.inputFileName,options.scoreCutoff,options.blockCutoff)
	obj.psl2bedFile(options.outputFileName)
            
#if __name__== '__main__':
#	main()
#else:
#    print >>sys.stderr, "module " + __name__ + " imported!"

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
