import sys,re,operator,os
from string import *
from shutil import *
from math import *

inputPSI =sys.argv[1]
inputTOPO=sys.argv[2]

fh_psi  = open(inputPSI , 'r')
fh_topo = open(inputTOPO, 'r')

psi_lst   =[]
psi_lines = fh_psi.readlines()

for lines in psi_lines:
	if re.compile('\s*\d').match(lines):
	    lines = lines.strip()
            items = lines.split()
            #print items [2:22]
            psi_lst.append(items[1])

topo_lines = fh_topo.readlines()
fh_topo.close()
if len(topo_lines)<1:
	fh_topo = open(inputTOPO, 'w')
	fh_topo.write('>query'+'\n')
	for i in xrange(len(psi_lst)):
		fh_topo.write('o')

	fh_topo.write('\n')

fh_topo.close()
