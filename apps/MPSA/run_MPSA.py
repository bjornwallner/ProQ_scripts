#!/usr/bin/python 
import sys,re,operator,os, commands
from string import *
from shutil import *
from math import *

### Don't blame me this is messy master student code.... :-) /BW

#install_path = '/local/www/services/ProQM/'
pathname = os.path.dirname(sys.argv[0]) + "/../"
install_path = os.path.abspath(pathname) + "/"
#print pathname + '\n'
#print install_path + '\n'

path_scripts    = install_path + 'MPSA/SCRIPTS/'
path_scripts_new    = install_path + 'MPSA/' # ('/Library/WebServer/ProQM/MPSA/')
#path_bin_scripts=('/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/MPSA/SCRIPTS/bin/')#('/home/simsuede/mpRAP/bin/') #/var/www/html/mprap/SCRIPTS/bin/')
#path_scampi     =('/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/scampi-seq/')
path_SVM        =install_path + 'MPSA/bin/' #('/Library/WebServer/ProQM/bin/') #/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/MPSA/SCRIPTS/SVM/')
path_SVMmodel   =install_path + 'MPSA/'  #('/Library/WebServer/ProQM/MPSA/') #('/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/MPSA/SCRIPTS/SVM/')
#uniref_database =install_path + 'DB/uniref90.fasta' #('/Library/WebServer/ProQM/DB/uniref90.fasta') #/afs/pdc.kth.se/home/s/simsuede/Public/uniref90.fasta') #/afs/pdc.kth.se/home/y/yogesh/disk3/tmp/uniref90.fasta')

#path_input      =('/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/MPSA/temp_input/')
#path_output     =('/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/MPSA/temp_output/')

temp_dir=commands.getoutput('mktemp -d /tmp/MPSA.XXXXXXXXX')
temp_dir=temp_dir+'/'

path_blastpgp   =install_path +'blast-2.2.18_x86_64/bin/' #('/Library/WebServer/ProQM/blast-2.2.16/bin/')
#database        =install_path + 'DB/uniref90.fasta' #('/Library/WebServer/ProQM/DB/uniref90.fasta') #/afs/pdc.kth.se/home/b/bjornw/.vol/bjornw31/DB/uniref90.fasta')
#bindir          =('/afs/pdc.kth.se/home/k/kriil/D_vol/SCRIPTS/bin/')
formatdb        =install_path + 'blast-2.2.18_x86_64/bin/formatdb'#('/Library/WebServer/ProQM/blast-2.2.16/bin/formatdb')
zpred           =install_path + 'zpred/bin/zpred.pl' #('/Library/WebServer/ProQM/zpred/bin/zpred.pl')
infile_AA  = sys.argv[1] # 
path_output= sys.argv[2] # 
uniref_database = sys.argv[3]
database = uniref_database


#for filename in dirlist:

path_lst = infile_AA.split('/')
filename = path_lst[len(path_lst)-1:][0]
key      = filename.replace('.fasta','')


#query1=path_input+filename.replace('.fa','')
#print infile_AA
#print filename + " " + key

query2=re.sub('.fasta$','',temp_dir+filename)
#zpred_file=infile_AA.replace('.fasta$',".zpred")
#topo_file=infile_AA.replace('.fasta$',".topcons.fa")
zpred_file=re.sub('.fasta$',".zpred",infile_AA)
topo_file=re.sub('.fasta$',".topcons.fa",infile_AA)
#print zpred_file + '\n' + topo_file + '\n'
#sys.exit()


########################################################################################
# RUN PSI_BLAST AS in TOPCONS
print path_blastpgp+'blastpgp'+' '+'-i'+' '+infile_AA+' '+'-d'+' '+database+' '+'-e 1e-3 -v 500 -b 500 -T F -m 6 -a 2 >'+' '+query2+'.blast'

argsystem= os.system(path_blastpgp+'blastpgp'+' '+'-i'+' '+infile_AA+' '+'-d'+' '+database+' '+'-e 1e-3 -T F -v 500 -b 500 -m 6 -a 8 >'+' '+query2+'.blast')
	
#print path_scripts+'msa62fasta_oneround.pl'+' '+query2+'.blast'+' '+'>'+' '+query2+'.hits.db'
print path_scripts+'msa62fasta_oneround.pl'+' '+query2+'.blast'+' '+'>'+' '+query2+'.hits.db'
argsystem = os.system(path_scripts+'msa62fasta_oneround.pl'+' '+query2+'.blast'+' '+'>'+' '+query2+'.hits.db')

#print formatdb+' ' +'-i'+' '+query2+'.hits.db'+' '+'-l /dev/null'
print formatdb+' ' +'-i'+' '+query2+'.hits.db'+' '+'-l /dev/null'
argsystem= os.system(formatdb+' ' +'-i'+' '+query2+'.hits.db'+' '+'-l /dev/null')	

#print path_blastpgp+'blastpgp'+' '+'-j 2'+' '+'-i'+' '+infile_AA+' '+'-d'+' '+query2+'.hits.db'+' '+'-e 1000000'+' '+'-v 1000000'+' '+'-b 1000000'+' '+'-a 4'+' '+'-C'+' '+query2+'.chk'+' '+'-Q'+' '+query2+'.psi'+' '+'>/dev/null'
print path_blastpgp+'blastpgp'+' '+'-j 2'+' '+'-i'+' '+infile_AA+' '+'-d'+' '+query2+'.hits.db'+' '+'-e 1000000'+' '+'-v 1000000'+' '+'-b 1000000'+' '+'-a 8'+' '+'-C'+' '+query2+'.chk'+' '+'-Q'+' '+query2+'.psi'+' '+'>/dev/null'
argsystem= os.system(path_blastpgp+'blastpgp'+' '+'-j 2'+' '+'-i'+' '+infile_AA+' '+'-d'+' '+query2+'.hits.db'+' '+'-e 1000000'+' '+'-v 1000000'+' '+'-b 1000000'+' '+'-a 8 -T F'+' '+'-C'+' '+query2+'.chk'+' '+'-Q'+' '+query2+'.psi'+' '+'>/dev/null')


#######################################################################################


################   #######################################################################################
################   # USE the input sequence to generate SCAMPI TOPOLOGY (needed for zpred and web server output)
################   argsystem= os.system('perl'+' '+path_scampi+'SCAMPI_run.pl'+' '+infile_AA)
################   
################   argsystem= os.system('mv'+' '+infile_AA+'.res'+' '+path_output)
################   
################   argsystem= os.system('python'+' '+path_scampi+'compacttopo2top.py'+' '+query2+'.fa.res'+' '+'>'+' '+query2+'.topo')
################   # #script needed for checking and printing possible error messages
################   
################   argsystem= os.system('python'+' '+path_scripts+'topo_checker.py'+' '+query2+'.psi'+' '+query2+'.topo')
################   #####################################################################################
################   
################   profile_file=query2+'.psi'
################   topo_file   =query2+'.topo'
################   zpred_file  =infile_AA+'.zpred'
################   
################   argsystem= os.system(zpred +' -mode profile_topology -profile '+ profile_file +' -topology '+topo_file+' -prediction modhmm -out '+ zpred_file+' -tmpdir '+temp_dir )
################   
################   #print zpred +' -mode profile_topology -profile '+ profile_file +' -topology '+topo_file+' -prediction modhmm -out '+ zpred_file+' -tmpdir '+temp_dir 

########################################################################################
# USE PSI- AND TOPO-FILES AS INPUT TO RUN ZPRED.pl from TOPCONS
#argsystem= os.system('perl'+' '+'/afs/pdc.kth.se/home/k/kriil/vol_03/Programs/zpred/bin/zpred.pl'+' '+' '+'-mode profile'+' '+'-profile'+' '+query2+'.psi'+' '+'-topology'+' '+query2+'.topo'+' '+'-prediction modhmm'+' '+'-out'+' '+query2+'.zpred'+' '+'/afs/pdc.kth.se/home/k/kriil/D_vol/Programs/MPSA/tmp/')


########################################################################################

#######################################################################################
#########SCRIPTS FOR MSA AND SHANNON ENTROPY CALCULATION##############################
# CONVERT BLAST OUTPUT TO MSA IN FASTA-FORMAT

#print path_scripts+'blast_m62fasta_msa_S.pl'+' '+query2+'.blast'+' '+infile_AA+' '+'>'+' '+query2+'.msa.fa'
print path_scripts+'blast_m62fasta_msa_S.pl'+' '+query2+'.blast'+' '+infile_AA+' '+'>'+' '+query2+'.msa.fa'
argsystem= os.system(path_scripts+'blast_m62fasta_msa_S.pl'+' '+query2+'.blast'+' '+infile_AA+' '+'>'+' '+query2+'.msa.fa')


# FIX SO MSA ARE CORRECTLY ALIGNED	
print 'perl'+' '+path_scripts+'align_msa2pdb_S2.pl'+' '+query2+'.msa.fa'+' '+infile_AA+' '+'>'+' '+query2+'.2_msa.fa'
argsystem= os.system('perl'+' '+path_scripts+'align_msa2pdb_S2.pl'+' '+query2+'.msa.fa'+' '+infile_AA+' '+'>'+' '+query2+'.2_msa.fa')

#print "get shannon"

# CALCULATE SHANNON (CONSERVATION) SCORE FROM MSA-FILE	
print 'perl'+' '+path_scripts+'get_Shannon_S.pl'+' '+query2+'.2_msa.fa'+' '+infile_AA+' '+'>'+' '+query2+'.shann'
argsystem= os.system('perl'+' '+path_scripts+'get_Shannon_S.pl'+' '+query2+'.2_msa.fa'+' '+infile_AA+' '+'>'+' '+query2+'.shann')
# SCRIPT THAT PRINTS ERROR for the web server in case of missing blast hits or MSA	
print 'python'+' '+path_scripts+'MSA_checker.py'+' '+query2+'.shann'+' '+query2+'.finalout'
argsystem= os.system('python'+' '+path_scripts+'MSA_checker.py'+' '+query2+'.shann'+' '+query2+'.finalout')
#####################################################################################

##############################################################################
#MAKE SVM TABLE using psi blast score, zpred, and shannon (conservation) windows size 9 
print 'python'+' '+path_scripts_new+'make_SVM_table.py'+' '+query2+'.psi'+' '+query2+'.shann'+' ' + zpred_file +' '+query2+'.table.1'+' '+query2+'.table.2'
argsystem= os.system('python'+' '+path_scripts_new+'make_SVM_table.py'+' '+query2+'.psi'+' '+query2+'.shann'+' ' + zpred_file +' '+query2+'.table.1'+' '+query2+'.table.2')
#RUN SVM CLASSIFY using the previous generated table as input and an SVM model file
print 'python'+' '+path_scripts_new+'mprap_SVM_classify.py'+' '+query2+'.table.2'+' '+query2+'.class'+' '+path_SVM+' '+path_SVMmodel
argsystem= os.system('python'+' '+path_scripts_new+'mprap_SVM_classify.py'+' '+query2+'.table.2'+' '+query2+'.class'+' '+path_SVM+' '+path_SVMmodel)

############################################################################
#GENERATE THE OUTPUT FOR THE WEB SERVER #!!!!!!!!!!!!! #COMMENT THIS IF SIMPLE AFS output	
#argsystem= os.system('python'+' '+path_scripts+g'enerate_webserver_output_mprap.py'+' '+query2+'.class'+' '+query2+'.psi'+' '+query2+'.zpred'+' '+query2+'.shann'+' '+query2+'.topo'+' '+query2+'.blast'+' '+query2+'.finalout'+' '+query2+'.mpRAP')

#or...

#GENERATE THE OUTPUT FOR local AFS #!!!!!!!!!!!!! #COMMENT THIS from WEB server
print 'python'+' '+path_scripts_new+'generate_afs_output_mprap.py'+' '+query2+'.class'+' '+query2+'.psi'+' '+query2+'.zpred'+' '+query2+'.shann'+' '+topo_file+' '+query2+'.blast'+' '+query2+'.mpSA'
argsystem= os.system('python'+' '+path_scripts_new+'generate_afs_output_mprap.py'+' '+query2+'.class'+' '+query2+'.psi'+' '+query2+'.zpred'+' '+query2+'.shann'+' '+topo_file+' '+query2+'.blast'+' '+query2+'.mpSA')
argsystem= os.system('cat'+' '+query2+'.mpSA')
###########################################################################

# remove all temporany files

#query2+'.finalout'

#argsystem= os.system('rm'+' '+query2+'.blast'+' '+query2+'.topo'+' '+query2+'.class'+' '+query2+'.shann'+' '+query2+'.zpred'+' '+query2+'.psi'+' '+query2+'.table.2'+' '+query2+'.table.1'+' '+query2+'.2_msa.fa'+' '+query2+'.msa.fa'+' '+query2+'.fa.res'+' '+query2+'.chk'+' '+query2+'.hits.db'+' '+query2+'.hits.db.phr'+' '+query2+'.hits.db.pin'+' '+query2+'.hits.db.psq')


#argsystem= os.system('rm -fr'+' '+query2+'.blast'+' '+query2+'.topo'+' '+query2+'.class'+' '+query2+'.zpred'+' '+query2+'.psi'+' '+query2+'.table.2'+' '+query2+'.table.1'+' '+query2+'.2_msa.fa'+' '+query2+'.msa.fa'+' '+query2+'.fa.res'+' '+query2+'.chk'+' '+query2+'.hits.db'+' '+query2+'.hits.db.phr'+' '+query2+'.hits.db.pin'+' '+query2+'.hits.db.psq')
#
#
#argsystem= os.system('cp'+' '+query2+'.mpSA'+' '+path_output)
#argsystem= os.system('rm'+' '+query2+'.mpSA') #+' '+path_output)
os.system('rm -fr'+' '+temp_dir);#+' '+path_output)


