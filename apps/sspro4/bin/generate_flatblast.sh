#!/bin/sh
if [ $# -ne 2 ]
then
	echo "need two parameters:seq_file, output_file." 
	exit 1
fi
/nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/script/generate_flatblast.pl /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/blast2.2.8/ /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/script/ /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/data/big/big_98_X /nfs/bjornw/Research/git/ProQ_scripts/apps/sspro4/data/nr/nr $1 $2 
