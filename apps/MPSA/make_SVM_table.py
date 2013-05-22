#!/usr/bin/python 
import sys,re,operator,os
from string import *
from shutil import *
from math import *

# input files needed to generate the input vector for SVM classify
PSIblast=sys.argv[1]
SHANNON=sys.argv[2]
ZPRED=sys.argv[3]
FILE1=sys.argv[4]
FILE2=sys.argv[5]

########################################
# reading the files
fh_sha = open(SHANNON, 'r')
fh_zpred = open(ZPRED, 'r')
fh_psi = open(PSIblast, 'r')
 # writing two temporany output files
fh = open(FILE1, 'w')
fh2 = open(FILE2, 'w')
########################################

aa_lst=[]
sha_lst=[]
zpred_lst=[]
psi_lst=[]
sha_lines = fh_sha.readlines()
zpred_lines = fh_zpred.readlines()
psi_lines = fh_psi.readlines()

#######################################################################
# parsing and generating the needed list from the input files (psi_blast, scampi and zpred outputs)
#shannon
for i in xrange(len(sha_lines)):
	key = sha_lines[i]
        if key.startswith('>NORM_Shannon_noindel'):
            sha_lst = sha_lines[i+1].replace('\n','').split('_')
        elif key.startswith('>AA'):
            aa_lst = sha_lines[i+1].strip()
#zpred
for lines in zpred_lines:
	lines = lines.strip()
        items = lines.split()
        zpred_lst.append(items[0])        
#psi_balst
for lines in psi_lines:
	if re.compile('\s*\d').match(lines):
	    lines = lines.strip()
            items = lines.split()
            psi_lst.append(items[2:22])
#aa
#HACK
#i=0
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#
#
#ii=0
for i in xrange(len(aa_lst)): 
	fh.write('0')
	for j in xrange(len(psi_lst[i])):
            fh.write('\t' + psi_lst[i][j] )
        fh.write('\t'+ sha_lst[i])
	fh.write('\t'+ zpred_lst[i]+'\n')
	ii=i


#i=ii
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')
#fh.write('0')
#for j in xrange(len(psi_lst[i])):
#	fh.write('\t' + '0' )
#fh.write('\t'+ '0')
#fh.write('\t'+ zpred_lst[i]+'\n')

fh.close()




###################################################################
##Creating SVM vector of symmetric windows size 9 from the list files above
fh = open(FILE1, 'r')
listLin=[]
listLin=fh.readlines()
list_items=[]
nlist=[]
t_list=[]
xlist=[]
TF_item=[]
fh2.write("#"+FILE2+'\n')
for aline in listLin:
	TF_item=[]
	if not re.compile("#").match(aline):
		aline=aline.strip()
		items=aline.split()
		list_items.append(items[1:])
for i in range(0,len(list_items),1):
	if i==0:
		h=1
		values=[]
		xlist=list_items[0]+list_items[0]+list_items[0]+list_items[0]+list_items[0]+list_items[1]+list_items[2]+list_items[3]+list_items[4]
		for j in range(0,len(xlist),1):
			h+=1
	elif i==1:
		h=1
		values=[]
		xlist=list_items[1]+ list_items[0]+list_items[0]+list_items[0]+list_items[0]+list_items[2]+list_items[3]+list_items[4]+list_items[5]
		for j in range(0,len(xlist),1):
			h+=1
	elif i==2:
		h=1
		values=[]
		xlist=list_items[2]+ list_items[1]+list_items[0]+list_items[0]+list_items[0]+list_items[3]+list_items[4]+list_items[5]+list_items[6]
		for j in range(0,len(xlist),1):
			h+=1
	elif i==3:
		h=1
		values=[]
		xlist=list_items[3]+ list_items[2]+list_items[1]+list_items[0]+list_items[0]+list_items[4]+list_items[5]+list_items[6]+list_items[7]			
		for j in range(0,len(xlist),1):
			h+=1
	elif i==(len(list_items)-4):
		h=1
		values=[]
		xlist=list_items[i]+ list_items[i-5]+list_items[i-6]+ list_items[i-7]+list_items[i-8]+list_items[i+1]+list_items[i+2]+list_items[i+3]+list_items[0]
		for j in range(0,len(xlist),1):
			h+=1
	elif i==(len(list_items)-3):
		h=1
		values=[]
		xlist=list_items[i]+ list_items[i-4]+list_items[i-5]+ list_items[i-6]+list_items[i-7]+list_items[i+1]+list_items[i+2]+list_items[0]+list_items[0]		
		for j in range(0,len(xlist),1):
			h+=1
	elif i==(len(list_items)-2):
		h=1
		values=[]
		xlist=list_items[i]+ list_items[i-3]+list_items[i-4]+ list_items[i-5]+list_items[i-6]+list_items[i+1]+list_items[0]+list_items[0]+list_items[0]
		for j in range(0,len(xlist),1):
			h+=1             
	elif i==(len(list_items)-1):
		h=1
		values=[]
		xlist=list_items[i]+list_items[i-2]+list_items[i-3]+ list_items[i-4]+list_items[i-5]+list_items[0]+list_items[0]+list_items[0]+list_items[0]
		for j in range(0,len(xlist),1):
			h+=1
	else:
		h=1
		values=[]
		nlist=[]
		fh2.write('0'+"\t")
		nlist=list_items[i]+list_items[i-1]+list_items[i-2]+ list_items[i-3]+list_items[i-4]+list_items[i+1]+list_items[i+2]+list_items[i+3]+list_items[i+4]
		for j in range(0,len(nlist),1):
			fh2.write(repr(h)+':'+nlist[j]+"\t")
			h+=1
		fh2.write("\n")
#########################################################################################################################################

# the final output is a vactor of windows size9 that include the PSI_blast score, the Shannon Entropy and the |zpred coordinates| (ZPRED 1.0)
 
