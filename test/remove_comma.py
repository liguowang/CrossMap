import os, sys
infile = sys.argv[1]

for l in open(infile,'r'):
	l = l.strip()
	if l.startswith('#'):
		print (l)
	else:
		f = l.split()
		if ',' in f[4]:
			continue
		else:
			print (l)