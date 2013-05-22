
SSpro/ACCpro 4.01: Secondary Structure and Solvent Accessibility Predictors

Reference:
J. Cheng, A. Randall, M. Sweredoski, P. Baldi, SCRATCH: a Protein Structure and
 Structural Feature Prediction Server, Nucleic Acids Research, 
Web Server Issue, vol. 33, 72-76, 2005

Contact:
Jianlin Cheng 
Ph.D. Candidate
School of Information and Computer Sciences
University of California Irvine
jianlinc@uci.edu

################################################################################
Install SSpro4/ACCpro on Linux (Tested on Redhat Linux) 
################################################################################

1) unzip the tarball 
   e.g. tar xzf sspro4.tar.gz

2) change the directory to unzipped sspro4 directory 
   e.g. cd sspro4

3) open configure.pl and set the $install_dir to the sspro4 installation dir.
   e.g. /home/your_home_dir/sspro4/

4) execute configure.pl to configure and install the package. 
   e.g. ./configure.pl

   Installaltion is done!  

###################################################################
Test Installation
###################################################################

1) test secondary structure predictions on sample sequences (1aqta or 1s12)
   with homology analysis:

   cd test
   ../bin/predict_ssa.sh 1aqta.fasta 1aqta.test

   The prediction of secondary structure is stored in 1aqta.test file. 
   Compare 1aqta.test with 1aqta.ss in the
   test directory. They should be identical. 

2) test secondary structure predictions using neural network only (ab-initio):

   cd test
   ../bin/predict_ssa_ab.sh 1aqta.fasta 1aqta.test

   The prediction of secondary structure is stored in 1aqta.test file. 
   Compare 1aqta.test with 1aqta.ss.ab in the
   test directory. They should be identical. 

3) test solvent accessibility predictions on a sample sequence with homology
   analysis:

   cd test
   ../bin/predict_acc.sh 1aqta.fasta test.acc

   Compare test.acc with 1aqta.acc. They should be identical. 
   predict_acc.sh combines results from neural networks and homology.

4) test solvent accessibility predictions using neural networks only
  (ab-initio)

   cd test
   ../bin/predict_acc.sh 3chya out.acc

   out.acc should be the same as 3chya.acc.ab


##################################################################
Descriptions of sub directories
##################################################################

bin: shell scripts
model: model files of neural network
script: perl script (homology script and main script)
data: the database (big and non-redudandant database and pdb database)
server: the sspro and accpro executables
blast2.2.8: the blast tool (version 2.2.8 is used)
test: test directory. 

###################################################################
Usage of SSpro and ACCpro
###################################################################

Commands: 
For secondary structure prediction:
  path/sspro4/bin/predict_ssa.sh    sequence.fasta output_file

For solvent accessibility prediction at 25% threshold using both neural
networks and homology information:
  path/sspro4/bin/predict_acc.sh    sequence.fasta output_file

  The output file format: 
  name
  sequence
  predicted secondary structure or solvent accessibility(e: exposed, -: buried
  at 25% threshold)

  See 1aqta.fasta, 1aqta.ss, 1aqta.acc  examples in test directory. 

  NOTICE: The sequence in the fasta file should be one single line. The 
  sequence can be up to a few thousand residues long.  

For solvent accessibility prediction at 25% threshold using neural networks
only (ab-initio approach)

  path/sspro4/bin/predict_acc_ab.sh    sequence.fasta output_file

#####################################################################
Predict solvent accessibility at the thresholds other than 25%
#####################################################################

To predict solvent accessibility at other thresholds, use script
bin/predict_acc_multi.sh

predict_acc_multi.sh input_fasta_file  output_file threshold

threshold is an integer index between 0 and 19. 
0  ->  0%
1  ->  5%
...
5  -> 25%
...
19 -> 95%.

e.g. if you want 30% threshold, change it to 6.

The threshold = integer  * 5%. 

#####################################################################
		Release notes
#####################################################################

4.03: released on 1/19/2006
Add ab-initio prediction script (using neural network only, no 
homology) for secondary structure

4.02: released on 12/17/2005
Add ab-initio prediction script (using neural network only, no
homology)  for solvent accessibility

4.01: released on 8/12/2005

Fix a bug about solvent accessibility prediction for thresholds other than 25%.

Add a new script predict_acc_multi.pl to predict solvent accessibility 
at different thresholds. Don't need to change parameters manually anymore. 

I'd like to thank Dan Klose for the question leading to the bug fix.


4.0: released on 5/1/2005
First released version

-----------------------------------------------------------------------




