#!/bin/sh
if [ $# -ne 4 ]
then
	echo "need 4 parameters:seq_file(name,seq in compact format), align_dir, output file, solvent acc threshold index(0:0%,1:5%,..., 19:95%)." 
	exit 1
fi
#assumption: alignment_file = align_dir + name
/nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/script/predict_seq_sa_multi.pl /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/server/predict_seq_sa_multi.sh $1 $2 $3 $4 
