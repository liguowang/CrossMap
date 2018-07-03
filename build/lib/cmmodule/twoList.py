'''manipulate ndarray list'''                                                       
from itertools import imap,starmap,izip
from operator import mul,add,sub

def check_list(v1,v2):
	'''check if the length of two list is same'''
	if v1.size != v2.size:
		raise ValueError,"the lenght of both arrays must be the same"
	pass

def Add(v1,v2):
	'''add two list'''                                                                                 
	check_list(v1,v2)
	return v1.__add__(v2)

def Subtract(v1,v2):
	'''subtract v2 from v1'''                                                                                 
	check_list(v1,v2)
	return v1.__sub__(v2)

def Product(v1,v2):
	'''return product of two list'''                                                                                 
	check_list(v1,v2)
	return v1.__mul__(v2)

def Division(v1,v2):
	'''return divide v1 by v2. add 1 to both v1 and v2'''                                                                                 
	check_list(v1,v2)
	return (v1+1).__div__(v2+1)

def Average(v1,v2):
	'''return arithmetic mean of two list'''                                                                                 
	check_list(v1,v2)
	return v1.__add__(v2)/2

def geometricMean(v1,v2):
	'''return geometric mean of two list'''                                                                                 
	check_list(v1,v2)
	return (v1.__mul__(v2))**0.5

def Max(v1,v2):
	'''pairwise comparison two list. return  the max one between two paried number'''                                                                                 
	check_list(v1,v2)
	return imap(max,izip(v1,v2))

def Min(v1,v2):
	'''pairwise comparison two list. return  the max one between two paried number'''                                                                                 
	check_list(v1,v2)
	return imap(min,izip(v1,v2))
def euclidean_distance(v1,v2):
	'''return euclidean distance'''                                                                                 
	check_list(v1,v2)
	return (sum((v1.__sub__(v2))**2) / v1.size)**0.5