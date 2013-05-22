#!/bin/sh
if [ $# -ne 2 ]
then
	echo "need two parameters:seq_file, output_file." 
	exit 1
fi
/home/bjornw/BAR99/ProQ_scripts/apps/sspro4/script/generate_flatblast.pl /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/blast2.2.8/ /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/script/ /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/data/big/big_98_X /home/bjornw/BAR99/ProQ_scripts/apps/sspro4/data/nr/nr $1 $2 
