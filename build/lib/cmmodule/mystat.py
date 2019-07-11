#!/usr/bin/env python

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

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="3.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def RSS(arg):
    '''calculate Square root of sum of square. Input is ',' separated numbers'''
    lst=arg.split(',')
    lst_sum=0
    for i in [ int(i)**2 for i in lst]:
        lst_sum += i 
    #nsr=10*math.log10((1+noi_sum**0.5)/(1+sig_sum**0.5))
    return lst_sum**0.5
    
def H_mean(arg):
    '''calculate harmornic mean. Input is ',' separated numbers'''
    lst=[1/float(i) for i in arg.split(',') if float(i) !=0]
    if len(lst) == 0:
        return "NA"
    else:
        return len(lst)/(sum(lst))

def shannon_entropy(arg):
    '''calculate shannon's entropy (or Shannon-Wiener index).'''
    lst=arg
    lst=[float(i) for i in lst if float(i)>0]
    entropy=0.0
    for i in lst:
        entropy += (i/sum(lst)) * math.log((i/sum(lst)))
    if entropy == 0:
        return 0
    else:
        return -entropy

    
def shannon_entropy_es(arg):
    '''calculate estimator of shannon's entropy (Chao & Shen, 2003)'''
    lst=arg
    lst=[float(i) for i in lst if float(i)>0]
    if sum(lst)<=0 or min(lst)<0:return "NA"    #if there is no fragmental splicing
    if (len(lst)==1): return 0                  #if there is only 1 fragmental splicing
    lst.append(2)
    
    #estimate C_bar
    singleton=0
    entropy=0.0
    for i in lst:
        if i ==1:singleton +=1
    
    C_bar = 1- (singleton/sum(lst))
    for i in lst:entropy += ( (C_bar*i/sum(lst)) * math.log((C_bar*i/sum(lst))) )/(1-(1-C_bar*i/sum(lst))**sum(lst))
    if entropy == 0:
        return 0
    else:
        return -entropy

def shannon_entropy_ht(arg):
    '''calculate estimator of shannon's entropy based on Horzitz-Thompson'''
    lst=arg.split(',')
    lst=[float(i) for i in lst if float(i)>0]
    if sum(lst)<=0 or min(lst)<0:return "NA"    #if there is no fragmental splicing
    if (len(lst)==1): return 0                  #if there is only 1 fragmental splicing
    
    #estimate C_bar
    entropy=0.0
    for i in lst:
        entropy += ( (i/sum(lst)) * math.log((i/sum(lst))) )/(1-(1-i/sum(lst))**sum(lst))
    return -entropy
    
def simpson_index(arg):
    '''calculate Gini-Simpson's index. Input is ',' separated numbers'''
    lst=arg.split(',')
    lst=[float(i) for i in lst if float(i)>0]
    simpson=0.0
    
    try:
        for i in lst:
            simpson = simpson + (i/sum(lst))**2
        return 1-simpson
    except: return 0
    
def simpson_index_es(arg):
    '''calculate estimator Gini-Simpson's index. Input is ',' separated numbers'''
    lst=arg.split(',')
    lst=[float(i) for i in lst if float(i)>0]
    simpson=0.0
    
    try:
        for i in lst:
            simpson = simpson + i*(i-1)
        return 1- (simpson/(sum(lst)*(sum(lst)-1)))
    except: return 0
    
def Hill_number(arg,qvalue=1):
    '''Calculate real diversity (Hill's number). Input is ',' separated numbers. qvalue is the only
    parameter for Hill's function. When q=1, it return exp(H) which is the effective number of junctions
    calculated by Shannon's entropy. When q<1, Hill's function was favors low frequency junctions. 
    When q>1, Hill's function was favors high frequency junctions (common junctions). Simpon's Index
    is particular case of Hill's function as q=2'''
    
    lst=arg.split(',')
    lst=[float(i) for i in lst if float(i)>0]
    freq=[(i/sum(lst))**qvalue for i in lst]
    try:
        return (sum(freq))**(1/(1-qvalue))
    except:
        return math.exp(shannon_entropy(arg))
import math
import functools

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0 to 100.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent/100.0
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1
    
def percentile_list(N):
    """
    Find the percentile of a list of values.
    @parameter N - is a list of values. Note N MUST BE already sorted.
    @return - the list of percentile of the values
    """
    if not N:return None
    if len(N) <100: return N
    per_list=[]
    for i in range(1,101):
        k = (len(N)-1) * i/100.0
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            per_list.append( int(N[int(k)])  )
        else:
            d0 = N[int(f)] * (c-k)
            d1 = N[int(c)] * (k-f)
            per_list.append(int(round(d0+d1)))  
    return per_list
