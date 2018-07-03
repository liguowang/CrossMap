import math
def point_poip(actual, mean):
	'''give poisson pvalue. P[obs ==mean]'''
	# naive:   math.exp(-mean) * mean**actual / factorial(actual)
	# iterative, to keep the components from getting too large or small:
	p = math.exp(-mean)
	for i in xrange(actual):
		p *= mean
		p /= i+1
	return p

def cumu_poip(num, mean,logp=False):
	'''give poisson pvalue P[obs >=mean]'''
	s=0.0
	for i in range(0,num+1):
		s += point_poip(i,mean)
	if logp is True:
		try:
			return -10*math.log10(1-s)
		except:
			return 3000
	else:
		return 1-s
	
