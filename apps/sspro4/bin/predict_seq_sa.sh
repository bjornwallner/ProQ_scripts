#!/bin/sh
if [ $# -ne 3 ]
then
	echo "need three parameters:seq_file(name,seq in compact format), align_dir, output file." 
	exit 1
fi
#assumption: alignment_file = align_dir + name
/proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/script/predict_seq_sa.pl /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/server/predict_seq_sa.sh $1 $2 $3 
