#!/usr/bin/env python
'''manipulate fasta for fastq format files.'''

#import built-in modules
import re
import sys
from string import maketrans
from random import shuffle
from heapq import nlargest
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


def S_diff(lst):
	'''Given a list of int or float, calculate S_diff and S_point'''
	
	S_avg = sum(lst) / len(lst)
	S_dist = [i-S_avg for i in lst]	#distance to average
	S_cum=[]	#list of cumulative sum
	S_cum.append(0)
	for i in range(0,len(S_dist)):
		S_cum.append(S_cum[i] + S_dist[i])
	return [nlargest(1,range(0,len(S_cum)),key=lambda i: S_cum[i]),(max(S_cum) - min(S_cum))]
	#return the index of maximum_diff index, and maximum_diff

def bootstrap(lst,obs,rep=1000):
	'''Given a list of int or float (lst) and an observation value(obs). calcualte the chance (pvalue) 
	of getting this observation through bootstrapping.'''
	
	shuffled_diff=[]
	count=0
	tmp=lst
	for i in range(0,rep):
		shuffle(tmp)
		shuffled_diff.append(S_diff(tmp))
	
	for i in sorted(shuffled_diff):
		if (i>=obs):
			count += 1 
	if count/rep <0.5:
		return count/rep
	else:
		return 1- count/rep
