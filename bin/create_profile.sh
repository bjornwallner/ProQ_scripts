#!/bin/bash
#module add blast/2.0.11
# This is a simple script which will carry out all of the basic steps
# required to make a PSIPRED V2 prediction. Note that it assumes that the
# following programs are in the appropriate directories:
# blastpgp - PSIBLAST executable (from NCBI toolkit)
# makemat - IMPALA utility (from NCBI toolkit)
# psipred - PSIPRED V2 program
# psipass2 - PSIPRED V2 program

dirname=`dirname $0`;
dirname=`( cd $dirname && pwd )`
echo $dirname
basename=${1%.*}
outfile=$basename.mtx
echo $basename
#exit
# The name of the BLAST data bank
dbname=$2 # /local/www/services/ProQ2/DB/uniref90.fasta
cpu=$3
# Where the NCBI programs have been installed
ncbidir=$dirname/../apps/blast-2.2.18_x86_64/bin/
if [ ! -f ".ncbirc" ] 
    then
echo "[ncbi]" > .ncbirc
echo "Data=$dirname/../apps/blast-2.2.18_x86_64/data/" >> .ncbirc
fi
#set ncbidir = /usr/ebiotools/bin/



# Where the PSIPRED V2 programs have been installed
execdir=$dirname/../apps/psipred25/bin/

# Where the PSIPRED V2 data files have been installed
datadir=$dirname/../apps/psipred25/data/

if [ -f  $outfile ] 
    then
echo Already exist $outfile ...
exit
#else
#echo "Let's go!!!!"
fi

#echo $outfile
#echo $rootname
#echo $basename
#exit;


#\cp -f $1 psitmp.fasta
echo cp -f $1 psitmp.$$.fasta
cp -f $1 psitmp.$$.fasta

echo "Running PSI-BLAST with sequence" $1 "..."
echo $ncbidir/blastpgp -a $cpu -j 3 -h 0.001 -d $dbname -F F -i psitmp.$$.fasta -C psitmp.$$.chk -Q psitmp.$$.psi
#cpu=`grep -c CPU /proc/cpuinfo`;
$ncbidir/blastpgp -a $cpu -j 3 -h 0.001 -d $dbname -F F -i psitmp.$$.fasta -C psitmp.$$.chk -o psitmp.$$.blastpgp -Q psitmp.$$.psi > /dev/null #& $rootname.blast


echo "Running Makemat..."
echo psitmp.$$.chk > psitmp.$$.pn
echo psitmp.$$.fasta > psitmp.$$.sn
$ncbidir/makemat -P psitmp.$$
#cp psitmp.$$.mtx $rootname.mtx
#cp psitmp.$$.chk $rootname.chk
#cp psitmp.$$.psi $rootname.psi
#echo Cleaning up ...
#rm -f psitmp.$$.* error.log 
#exit;
echo "Predicting secondary structure..."
echo Pass1 ...

$execdir/psipred psitmp.$$.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > psitmp.$$.ss

echo Pass2 ...

$execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.0 psitmp.$$.ss2 psitmp.$$.ss 
cp psitmp.$$.ss2 $basename.ss2
cp psitmp.$$.mtx $basename.mtx
cp psitmp.$$.chk $basename.chk
$dirname/../bin/check_psiblast_matrix.py psitmp.$$.psi psitmp.$$.fasta psitmp.$$.psi.corrected
cp psitmp.$$.psi.corrected $basename.psi
cp psitmp.$$.blastpgp $basename.fasta.blastpgp # this is for accpro will skip the psiblast runs
# Remove temporary files
#exit
echo Cleaning up ...
rm -f psitmp.$$.* error.log
exit

#rm $rootname.ss $rootname.blast $rootname.ss2 #$rootname.horiz

#echo "Final output files:" $rootname.ss2 $rootname.horiz
echo "Final output files:" $rootname.horiz $rootname.ss2 
echo "Finished."
