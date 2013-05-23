#!/bin/sh
#predict the secondary stx for one sequence.
if [ $# -ne 3 ]
then
	echo "need three parameters:seq_file, ali_dir, data_format."
	exit 1
fi
#format: 1:9-line, 2:3-line, all seq must have a title line(1 20 3).

/nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/server/predict_seq_ss /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/model/sspro.model $1 $2 $3 
