#!/bin/bash -x 

num=$$

mkdir test.$num
cp tests_clean/1e12A_0001.pdb test.$num/
bin/run_all_external.pl -pdb test.$num/1e12A_0001.pdb 
cp tests_clean/1e12A_0001.subset.pdb test.$num/
bin/copy_features_from_master.pl test.$num/1e12A_0001.subset.pdb test.$num/1e12A_0001.pdb

#echo Membrane
#bin/run_all_external.pl -pdb test.$$/1e12A_0001.pdb -membrane 1
