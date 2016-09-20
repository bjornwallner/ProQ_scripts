#!/bin/sh
#predict ss, sa for a single sequence from scratch.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file." 
	exit 1
fi
#output: a file with predicted ss and sa.
/proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/script/predict_ss_sa.pl /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/blast2.2.8/ /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/data/big/big_98_X /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/data/nr/nr /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/server/predict_seq_ss.sh /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/server/predict_seq_sa.sh /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/script/ $1 $2 
