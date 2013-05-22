import sys,re,operator,os
from string import *
from shutil import *
from math import *

#mprap output directory
path_out=('/var/www/html/mprap/outdir/')

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
finalOUT=sys.argv[7]
finalOUT2=sys.argv[8]

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
#parse shannon output and list
fh3=open(inputSHANN,'r')
sha_lines = fh3.readlines()
sha_lst=[]
for i in xrange(len(sha_lines)):
	key = sha_lines[i]
	if key.startswith('>NORM_Shannon_noindel'):
		sha_lst = sha_lines[i+1].replace('\n','').split('_')
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
fh6=open(inputPSI,'r')
fh9=open(inputZPRED,'r')
fh11=open(finalOUT,'w')
listlin6=fh6.readlines()
listlin9=fh9.readlines()
list6=[]
list9=[]
list9_b=[]
####### WRITE mprap table
#zpred output
for aline in listlin9:
	aline=aline.strip()
	items=aline.split()
	list9.append(float(items[0]))
#PSI score andd aa output
for aline in listlin6:
	if re.compile('\s*\d').match(aline):
		aline=aline.strip()
		items=aline.split()
		list6.append(items[0:2])

length=repr(len(list6))

fh11.write("#RN"+'_'+'AA'+'_'+'Zcoord'+'_'+'Conservation'+'_'+'Accessibility'+'\n')
for i in range(0,len(column1),1):
	fh11.write(list6[i+4][0]+'_'+list6[i+4][1]+'_'+'%.1f' %(list9[i+4])+'_'+'%.3f' %(float(sha_lst[i+4]))+'_'+'%.1f' %(column1[i])+'\n')
##### WRITE mprap file output2 the table
fh12=open(finalOUT2,'w')
fh12.write('>AA_Len='+length+'\n')
for y in range(0,len(list6),1):
	fh12.write(list6[y][1])
fh12.write('\n'+'\n'+'>mpRAP'+'\n'+'XXXX')
for i in range(0,len(column1_b),1):
	fh12.write(column1_b[i])
fh12.write('XXXX'+'\n'+'\n')
fh12.write("#RN"+'\t'+'AA'+'\t'+'Topology'+'\t'+'Zcoord'+'\t'+'Conservation'+'\t'+'Accessibility'+'\n')
for k in range(0,len(column1),1):
	fh12.write(list6[k+4][0]+'\t'+list6[k+4][1]+'\t'+topo_lst2[k+4]+'\t'+'%.1f' %(list9[k+4])+'\t'+'%.3f' %(float(sha_lst[k+4]))+'\t'+'%.1f' %(column1[k])+'\n')
fh6.close()
fh9.close()
fh11.close()		
fh12.close()
######### put the files needed for the webserver output in a temporany directory
argsystem = os.system('cp'+' '+finalOUT2+' '+path_out)
argsystem = os.system('cp'+' '+inputBLAST+' '+path_out)
############################################################################################################################################################################

