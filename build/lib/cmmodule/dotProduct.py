import timeit # module with timing subroutines                                                              
import random # module to generate random numnbers                                                          
from itertools import starmap
from operator import mul


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="3.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"

def v(N=50,min=-10,max=10):
    """Generates a random vector (in an array) of dimension N; the                                          
    values are integers in the range [min,max]."""
    out = []
    for k in range(N):
        out.append(random.randint(min,max))
    return out

def check(v1,v2):
    if len(v1)!=len(v2):
        raise ValueError("the lenght of both arrays must be the same")
    pass

def d0(v1,v2):
    """                                                                                                     
    d0 is Nominal approach:                                                                                 
    multiply/add in a loop                                                                                  
    """
    check(v1,v2)
    out = 0
    for k in range(len(v1)):
        out += v1[k] * v2[k]
    return out

def d1(v1,v2):
    """                                                                                                     
    d1 uses an imap (from itertools)                                                                        
    """
    check(v1,v2)
    return sum(map(mul,v1,v2))

def d2(v1,v2):
    """                                                                                                     
    d2 uses a conventional map                                                                              
    """
    check(v1,v2)
    return sum(map(mul,v1,v2))

def d3(v1,v2):
    """                                                                                                     
    d3 uses a starmap (itertools) to apply the mul operator on an izipped (v1,v2)                           
    """
    check(v1,v2)
    return starmap(mul,zip(v1,v2))
