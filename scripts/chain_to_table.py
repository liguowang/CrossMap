#!/usr/bin/env python

import os,sys
import optparse
from textwrap import wrap
from time import strftime

import pyBigWig
import pysam
from bx.intervals import *
import numpy as np
import datetime
from cmmodule  import ireader
from cmmodule  import BED
from cmmodule  import annoGene
from cmmodule  import sam_header
from cmmodule  import bgrMerge

import collections


def printlog (mesg_lst):
    """
    print progress into stderr
    """
    if len(mesg_lst)==1:
        msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " +  mesg_lst[0]
    else:
        msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + ' '.join(mesg_lst)
    print(msg, file=sys.stderr)



def read_chain_file (chain_file, print_table=False):
    '''
    input chain_file could be either plain text, compressed file (".gz", ".Z", ".z", ".bz", ".bz2", ".bzip2"),
    or a URL pointing to the chain file ("http://", "https://", "ftp://"). If url was used, chain file must be plain text
    '''
    
    printlog(["Read chain_file: ", chain_file]),
    maps={}
    target_chromSize={}
    source_chromSize={}
    if print_table:
        blocks=[]
    
    for line in ireader.reader(chain_file):
        # Example: chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
        if not line.strip():
            continue
        line=line.strip()
        if line.startswith(('#',' ')):continue
        fields = line.split()
        
        if fields[0] == 'chain' and len(fields) in [12, 13]: 
            score = int(fields[1])        # Alignment score
            source_name = fields[2]       # E.g. chrY
            source_size = int(fields[3])  # Full length of the chromosome
            source_strand = fields[4]     # Must be +
            if source_strand != '+':
                raise Exception("Source strand in an .over.chain file must be +. (%s)" % line)
            source_start = int(fields[5]) # Start of source region
            source_end = int(fields[6])   # End of source region
        
            target_name = fields[7]       # E.g. chr5
            target_size = int(fields[8])  # Full length of the chromosome
            target_strand = fields[9]     # + or -
            target_start = int(fields[10])
            target_end = int(fields[11])
            target_chromSize[target_name]= target_size
            source_chromSize[source_name] = source_size
            
            if target_strand not in ['+', '-']:
                raise Exception("Target strand must be - or +. (%s)" % line)
            #if target_strand == '+':
            #   target_start = int(fields[10])
            #   target_end = int(fields[11])
            #if target_strand == '-':
            #   target_start = target_size - target_end
            #   target_end = target_size - target_start
            chain_id = None if len(fields) == 12 else fields[12]
            if source_name not in maps:
                maps[source_name] = Intersecter()           
            
            sfrom, tfrom = source_start, target_start
            
        # Now read the alignment chain from the file and store it as a list (source_from, source_to) -> (target_from, target_to)        
        elif fields[0] != 'chain' and len(fields) == 3: 
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
            if print_table:
                if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
                elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))
                
            if target_strand == '+':
                maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
            elif  target_strand == '-':
                maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
                
            sfrom += size + sgap
            tfrom += size + tgap
        
        elif fields[0] != 'chain' and len(fields) == 1: 
            size = int(fields[0])
            if print_table:
                if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
                elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))
                    
            if target_strand == '+':        
                maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
            elif target_strand == '-':
                maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
        else:
                raise Exception("Invalid chain format. (%s)" % line)
    if (sfrom + size) != source_end  or (tfrom + size) != target_end:
        raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)       
    
    if print_table:
        for i in blocks:
            print('\t'.join([str(n) for n in i]))
    
    return (maps,target_chromSize, source_chromSize)    
    
if len(sys.argv) != 2:
    print("Usage: python    chain_to_table.py   input_chain_file\n", file=sys.stderr)
else:   
    read_chain_file(sys.argv[1], print_table=True)  
