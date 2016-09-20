#!/bin/sh
if [ $# -ne 2 ]
then
	echo "need two parameters:seq_file, output_file." 
	exit 1
fi
/proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/script/generate_flatblast.pl /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/blast2.2.8/ /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/script/ /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/data/big/big_98_X /proj/wallner/users/x_bjowa/github/ProQ_scripts/apps/sspro4/data/nr/nr $1 $2 
