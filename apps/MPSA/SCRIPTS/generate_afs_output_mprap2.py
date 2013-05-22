import sys,re,operator,os
from string import *
from shutil import *
from math import *

#########
#opend input files from command line
outSVMclassify=sys.argv[1] # SVM classify output
AA_file=sys.argv[2]
finalOUT2=sys.argv[3]

##################################################################
##parsing the input files and generate lists used from the final output
#parse the classify output, showing the predicted accessibility real values and the converted binary (E=exposed > 25.0 B=buried otherwise) 
fh4=open(outSVMclassify,'r')
list_lin4=fh4.readlines()
column1=[]
column1_b=[]
for aline in list_lin4:
	aline= aline.strip()
	items= aline.split()
	column1.append(float(items[0]))
for i in range(0, len(column1),1):
	if float(column1[i])<=25.0:
		column1_b.append('B')
	else:
		column1_b.append('E')



fh3=open(AA_file,'r')
AA_lines = fh3.readlines()


fh12=open(finalOUT2,'w') # final afs termilnal output

fh12.write('>mpSA_bin'+'\n'+'XXXX')
for i in range(0,len(column1_b),1):
	fh12.write(column1_b[i])
fh12.write('XXXX'+'\n')
#write the predicted accessibility real values 
fh12.write('>mpSA_real'+'\n'+'X_X_X_X_')
for s in range(0,len(column1),1):
	fh12.write('%.1f' %(column1[s])) #number of decimal
	fh12.write('_')
fh12.write('X_X_X_X'+'\n')
# write aa in a fasta format
fh12.write('>AA\n')
fh12.write(AA_lines[1] )
fh12.close()
#####################################################################################################################################################################
