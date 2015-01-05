#!/bin/bash -x 

mkdir test.$$
cp tests_clean/1e12A_0001.pdb test.$$/
bin/run_all_external.pl -pdb test.$$/1e12A_0001.pdb 
cp tests_clean/1e12A_0001.subset.pdb test.$$/
bin/copy_features_from_master.pl test.$$/1e12A_0001.subset.pdb test.$$/1e12A_0001.pdb

#echo Membrane
#bin/run_all_external.pl -pdb test.$$/1e12A_0001.pdb -membrane 1
