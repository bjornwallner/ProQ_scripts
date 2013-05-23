#!/bin/sh
#predict sa for a single sequence from scratch.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file" 
	exit 1
fi
#output: a file with predicted ss and sa.
/nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/script/predict_acc_ab.pl /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/blast2.2.8/ /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/data/big/big_98_X /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/data/nr/nr /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/server/predict_seq_sa.sh /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/script/ $1 $2 
