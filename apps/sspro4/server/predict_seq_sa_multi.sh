#!/bin/sh
#predict the solvent acc for one sequence.
if [ $# -ne 4 ]
then
	echo "need four parameters:seq_file, ali_dir, data_format, solvent accessibility threshold index(0:0%,1:5%,...,5:25%,....19:95%)."
	exit 1
fi
#format: 1:9-line, 2:3-line, all seq must have a title line(1 20 3).

/nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/server/predict_seq_sa /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/model/accpro.model $1 $2 $3 $4 
