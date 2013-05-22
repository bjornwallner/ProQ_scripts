#!/bin/sh
#predict ss for a single sequence from scratch using ab-initio neural network.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file." 
	exit 1
fi
#output: a file with predicted ss and sa.
/home/bjornw/BAR99/ProQ_scripts/apps/sspro4/script/predict_ssa_ab.pl /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/blast2.2.8/ /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/data/big/big_98_X /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/data/nr/nr /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/server/predict_seq_ss.sh /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/script/ $1 $2 
