import sys,re,operator,os
from string import *
from shutil import *
from math import *

#########
#opend input files from command line
outSVMclassify=sys.argv[1] # SVM calssify output
inputPSI=sys.argv[2]
inputZPRED=sys.argv[3]
inputSHANN=sys.argv[4]
inputTOPO=sys.argv[5]
#
inputBLAST=sys.argv[6] # the blast output that is shown to the user of mprap webserver
# 
#output files used as webserver output
finalOUT2=sys.argv[7]

##################################################################
##parsing the input files and generate lists used from the final output
#parse the classify output, showing the predicted accessibility real values and the converted binary (E=exposed > 25.0 B=buried otherwise) 
fh4=open(outSVMclassify,'r')
list_lin4=fh4.readlines()
column1=[]
column1_b=[]

first=list_lin4[0];
last=list_lin4[len(list_lin4)-1]
list_lin4.insert(0,first)
list_lin4.insert(0,first)
list_lin4.insert(0,first)
list_lin4.insert(0,first)
list_lin4.append(last)
list_lin4.append(last)
list_lin4.append(last)
list_lin4.append(last)

for aline in list_lin4:
	aline= aline.strip()
	items= aline.split()
	column1.append(float(items[0]))
for i in range(0, len(column1),1):
	if float(column1[i])<=25.0:
		column1_b.append('B')
	else:
		column1_b.append('E')
#parse shannon output and list
fh3=open(inputSHANN,'r')
sha_lines = fh3.readlines()
sha_lst=[]
abs_sha_lst=[]
for i in xrange(len(sha_lines)):
	key = sha_lines[i]
	if key.startswith('>NORM_Shannon_noindel'):
		sha_lst = sha_lines[i+1].replace('\n','').split('_')
	if key.startswith('>ABS_Shannon_noindel'):
		abs_sha_lst = sha_lines[i+1].replace('\n','').split('_')
#parse topo output and list
fh13=open(inputTOPO,'r')
topo_lines = fh13.readlines()
topo_lst=[]
for i in xrange(len(topo_lines)):
	key = topo_lines[i]
	if key.startswith('>query'):
		topo_lst = topo_lines[i+1].replace('\n','').split()
topo_lst2=[]
for i in xrange(len(topo_lst)):
	topo_lst2=topo_lst[i]
#
fh3.close()
fh4.close()

###################################################################################################################################################################
########### this is the part for generating the output for the mpRAP webserver
######## output sample for terminal afs file
#	>mpRAP_bin
#	XXXXEEEEEEEEEEEEEBEEBBEEBBEBEEEEEEEEXXXX
#	>mpRAP_real
#	X_X_X_X_47.6_76.3_34.9_69.8_80.8_27.2_52.3_38.7_34.9_35.5_26.3_36.5_38.1_16.7_26.0_39.0_19.7_15.0_31.2_28.8_6.5_17.2_49.5_23.3_35.2_32.5_35.4_49.9_39.8_44.7_51.4_56.3_X_X_X_X
#	>AA_Len=40
#	VQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK
#######

fh6=open(inputPSI,'r') #input PSI for AA
fh12=open(finalOUT2,'w') # final afs termilnal output

listlin6=fh6.readlines()
list6=[]
length=repr(len(list6))
#take aa from psi.file
for aline in listlin6:
	if re.compile('\s*\d').match(aline):
		aline=aline.strip()
		items=aline.split()
		list6.append(items[0:2])

length=repr(len(list6))
#write the binary conversion like B E
#fh12.write('>mpRAP_bin'+'\n') #+'XXXX')
#for i in range(0,len(column1_b),1):
#	fh12.write(column1_b[i])
#fh12.write('\n')
#fh12.write('XXXX'+'\n')
#write the predicted accessibility real values 
#fh12.write('>mpRAP_real'+'\n') #+'X_X_X_X_')
for s in range(0,len(column1),1):
#	fh12.write('%1s %.1f', list6[s][1],column1[s]) #number of decimal
	fh12.write(list6[s][1] + ' ' + column1_b[s] + ' ' + str(column1[s]) + ' ' +sha_lst[s] + ' ' + abs_sha_lst[s]) #number of decimal

	fh12.write('\n')
#fh12.write('X_X_X_X'+'\n')
#fh12.write('\n')
# write aa in a fasta format
#fh12.write('>AA_Len='+length+'\n')
#for y in range(0,len(list6),1):
#	fh12.write(list6[y][1])
fh6.close()
fh12.close()
#####################################################################################################################################################################
