import sys,re,operator,os
from string import *
from shutil import *
from math import *

inputSHANN=sys.argv[1]
finalout=sys.argv[2]

fh_msa = open(inputSHANN, 'r')
msa_lines = fh_msa.readlines()
if len(msa_lines)<1:
	fh_finalout = open(finalout, 'w')
	fh_finalout.write('A BLAST ERROR HAS OCCURED! NO BLAST HITS HAVE BEEN FOUND FOR YOUR QUERY SEQUENCE!')
	
	fh_finalout.write('\n')
	fh_finalout.close()
	#print "exit"

else:
	print "MSA OK!"

