#!/usr/bin/env python
'''manipulate fastq files'''

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

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



class FQ:
	'''provides method to processing fastaq files'''
	
	def __init__(self,fqFile):
		'''This is constructor of FQ'''
		self.f=open(fqFile,'r')
		self.fileName=os.path.basename(fqFile)
		self.ABS_fileName=fqFile
