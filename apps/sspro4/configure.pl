#!/usr/bin/perl
###########################################################
#
#SSpro and ACCpro: Protein Secondary Structure/Solvent Accessibility Prediction Programs
#configure.pl: to configure the installation of PSpro
#
#Author: Jianlin cheng
#Date: April 8, 2004
#Institute for Genomics and Bioinformatics
#School of Information and Computer Science
#University of California, Irvine
#
# Modified: April 15, 2005

#References:
#J. Cheng, A. Randall, M. Sweredoski, P. Baldi, SCRATCH: a Protein Structure and Structural Feature Prediction Server, Nucleic Acids Research, Special Issue on Web Servers, in press, 2005.
##########################################################


use Cwd 'abs_path';
use File::Basename;
#######Customize settings here###############################################

#set installation directory of SSpro to your unzipped sspro4 directory 
#$install_dir = "/scratch/sspro4-bjornw/";
#$install_dir="/local/www/services/ProQ2/apps/sspro4/";
$install_dir=dirname(abs_path($0));
#exit;

##########End of Customizaton################################################





################Don't Change the code below##############

if (! -d $install_dir)
{
	die "can't find installation directory.\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
	$install_dir .= "/"; 
}

#set the fullpath of nr database
$nr_db_dir = "${install_dir}data/nr/";
$nr_db = "nr"; 
#set the fullpath of big database
$big_db_dir = "${install_dir}data/big/";
#set the big database name
$big_db = "big_98_X";

#if (! -d $nr_db_dir)
#{
#	die "can't find nr database directory.\n";
#}
#if ( substr($nr_db_dir, length($nr_db_dir) - 1, 1) ne "/" )
#{
#	$nr_db_dir .= "/"; 
#}
#if (! -d $big_db_dir)
#{
#	die "can't find big directory.\n";
#}
#if ( substr($big_db_dir, length($big_db_dir) - 1, 1) ne "/" )
#{
#	$big_db_dir .= "/"; 
#}

#check if the database are there. 
$nr_db_f = $nr_db_dir . $nr_db; 
$big_db_f = $big_db_dir . $big_db; 

#if (! -f "${big_db_f}.pal") {
#if (! -f "${big_db_f}.phr" || ! -f "${big_db_f}.pin" || ! -f "${big_db_f}.psq")
#{
#	die "can't find the big coil database.\n"; 	
#} 
#}
#
#if (! -f "${nr_db_f}.pal") {
#if (! -f "${nr_db_f}.phr" || ! -f "${nr_db_f}.pin" || ! -f "${nr_db_f}.psq" || ! -f "${nr_db_f}.pnd"
#    || ! -f "${nr_db_f}.pni" || ! -f "${nr_db_f}.psd" || ! -f "${nr_db_f}.psi" )
#{
#	die "can't find the non-redundant database.\n"; 	
#} 
#}
 
#check if the installation directory is right
#the configuration file must run in the installation directory
$cur_dir = `pwd`;  
chomp $cur_dir; 
$configure_file = "$cur_dir/configure.pl";
if (! -f $configure_file || $install_dir ne "$cur_dir/")
{
	die "Please check the installation directory setting and run the configure program there.\n";
}

$bin_dir = "${install_dir}bin/";
$data_dir = "${install_dir}data/";
$model_dir = "${install_dir}model/";
$script_dir = "${install_dir}script/";
$server_dir = "${install_dir}server/";
$test_dir = "${install_dir}test/";
$blast_dir = "${install_dir}blast2.2.8/"; 

if ( ! -d $bin_dir || ! -d $data_dir || ! -d $model_dir || ! -d $script_dir
   || ! -d $server_dir || ! -d $test_dir || ! -d $blast_dir )
{
	die "some sub directories don't exist. check the installation tar ball.\n";
}

$nr_db = ${nr_db_dir} . $nr_db;
$big_db = ${big_db_dir} . $big_db;

#generate model definition file, server script for sspro
$sspro_exe = "${server_dir}predict_seq_ss";
$sspro_sh = "${server_dir}predict_seq_ss.sh";
$sspro_model_dir = "${model_dir}sspro/";
$sspro_model_def = "${model_dir}sspro.model";

$accpro_exe = "${server_dir}predict_seq_sa";
$accpro_sh = "${server_dir}predict_seq_sa.sh";
$accpro_sh_multi = "${server_dir}predict_seq_sa_multi.sh";
$accpro_model_dir = "${model_dir}accpro/";
$accpro_model_def = "${model_dir}accpro.model";

print "generate sspro server script...\n";
open(SERVER_SH, ">$sspro_sh") || die "can't write sspro shell script.\n";
print SERVER_SH "#!/bin/sh\n#predict the secondary stx for one sequence.\n";
print SERVER_SH "if [ \$# -ne 3 ]\n";
print SERVER_SH "then\n\techo \"need three parameters:seq_file, ali_dir, data_format.\"\n\texit 1\nfi\n";
print SERVER_SH "#format: 1:9-line, 2:3-line, all seq must have a title line(1 20 3).\n\n";
print SERVER_SH "$sspro_exe $sspro_model_def \$1 \$2 \$3 \n"; 
close SERVER_SH;

print "generate accpro server script...\n";
open(SERVER_SH, ">$accpro_sh") || die "can't write accpro shell script.\n";
print SERVER_SH "#!/bin/sh\n#predict the solvent acc for one sequence. Don't change this script.\n";
print SERVER_SH "if [ \$# -ne 3 ]\n";
print SERVER_SH "then\n\techo \"need three parameters:seq_file, ali_dir, data_format.\"\n\texit 1\nfi\n";
print SERVER_SH "#format: 1:9-line, 2:3-line, all seq must have a title line(1 20 3).\n\n";
print SERVER_SH "$accpro_exe $accpro_model_def \$1 \$2 \$3 5\n"; 
close SERVER_SH;

print "generate accpro multi-threshold server script...\n";
open(SERVER_SH, ">$accpro_sh_multi") || die "can't write accpro shell script.\n";
print SERVER_SH "#!/bin/sh\n#predict the solvent acc for one sequence.\n";
print SERVER_SH "if [ \$# -ne 4 ]\n";
print SERVER_SH "then\n\techo \"need four parameters:seq_file, ali_dir, data_format, solvent accessibility threshold index(0:0%,1:5%,...,5:25%,....19:95%).\"\n\texit 1\nfi\n";
print SERVER_SH "#format: 1:9-line, 2:3-line, all seq must have a title line(1 20 3).\n\n";
print SERVER_SH "$accpro_exe $accpro_model_def \$1 \$2 \$3 \$4 \n"; 
close SERVER_SH;

print "generate sspro model definition file...\n";
opendir(MODEL_DIR, "$sspro_model_dir") || die "can't open the sspro model dir.\n";
open(MODEL_DEF, ">$sspro_model_def") || die "can't create sspro model def file.\n";
@file_list = readdir(MODEL_DIR);
closedir(MODEL_DIR); 
if (@file_list < 2)
{
	die "can't find sspro model file. check intallation tar ball.\n";
}
$model_num = @file_list;
$model_num -= 2; 
print MODEL_DEF "$model_num 3\n";
while (@file_list)
{
	$model_file = shift @file_list;
	if ($model_file ne "." && $model_file ne ".." && $model_file ne "CVS")
	{
		print MODEL_DEF "$sspro_model_dir$model_file\n";
	}
}
close MODEL_DEF; 

print "generate accpro model definition file...\n";
opendir(MODEL_DIR, "$accpro_model_dir") || die "can't open the accpro model dir.\n";
open(MODEL_DEF, ">$accpro_model_def") || die "can't create sspro model def file.\n";
@file_list = readdir(MODEL_DIR);
closedir(MODEL_DIR); 
if (@file_list < 2)
{
	die "can't find accpro model file. check intallation tar ball.\n";
}
$model_num = @file_list;
$model_num -= 2; 
print MODEL_DEF "$model_num 20\n";
while (@file_list)
{
	$model_file = shift @file_list;
	if ($model_file ne "." && $model_file ne ".." && $model_file ne "CVS")
	{
		print MODEL_DEF "$accpro_model_dir$model_file\n";
	}
}
close MODEL_DEF; 

############################Generate Shell Script For user##############################


#generate script of creating alignment for one sequence
open(SHELL, ">${bin_dir}generate_flatblast.sh") || die "can't create genereate_flatblast.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "if [ \$# -ne 2 ]\n";
print SHELL "then\n\techo \"need two parameters:seq_file, output_file.\" \n\texit 1\nfi\n";
print SHELL "${script_dir}generate_flatblast.pl $blast_dir $script_dir $big_db $nr_db \$1 \$2 \n";
close SHELL; 

#generate script of predicting ss for one sequence
open(SHELL, ">${bin_dir}predict_seq_ss.sh") || die "can't create predict_seq_ss.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "if [ \$# -ne 3 ]\n";
print SHELL "then\n\techo \"need three parameters:seq_file(name,seq in compact format), align_dir, output file.\" \n\texit 1\nfi\n";
print SHELL "#assumption: alignment_file = align_dir + name\n";
print SHELL "${script_dir}predict_seq_ss.pl ${server_dir}predict_seq_ss.sh \$1 \$2 \$3 \n";
close SHELL; 

#generate script of predicting sa for one sequence
open(SHELL, ">${bin_dir}predict_seq_sa.sh") || die "can't create predict_seq_ss.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "if [ \$# -ne 3 ]\n";
print SHELL "then\n\techo \"need three parameters:seq_file(name,seq in compact format), align_dir, output file.\" \n\texit 1\nfi\n";
print SHELL "#assumption: alignment_file = align_dir + name\n";
print SHELL "${script_dir}predict_seq_sa.pl ${server_dir}predict_seq_sa.sh \$1 \$2 \$3 \n";
close SHELL; 

#generate script of predicting sa for one sequence
open(SHELL, ">${bin_dir}predict_seq_sa_multi.sh") || die "can't create predict_seq_ss_multi.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "if [ \$# -ne 4 ]\n";
print SHELL "then\n\techo \"need 4 parameters:seq_file(name,seq in compact format), align_dir, output file, solvent acc threshold index(0:0%,1:5%,..., 19:95%).\" \n\texit 1\nfi\n";
print SHELL "#assumption: alignment_file = align_dir + name\n";
print SHELL "${script_dir}predict_seq_sa_multi.pl ${server_dir}predict_seq_sa_multi.sh \$1 \$2 \$3 \$4 \n";
close SHELL; 

#generate script of predicting ss one single sequence from scratch(generate alignments too) 
#use homology ss and sa predictor
open(SHELL, ">${bin_dir}predict_ssa.sh") || die "can't create predict_ssa.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "#predict ss, sa for a single sequence from scratch.\n";
print SHELL "if [ \$# -ne 2 ]\n";
print SHELL "then\n\techo \"need 2 parameters:seq_file(in fasta format), output_file.\" \n\texit 1\nfi\n";
print SHELL "#output: a file with predicted ss and sa.\n"; 
print SHELL "${script_dir}predict_ssa.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_ss.sh $script_dir \$1 \$2 \n";
close SHELL; 

#generate script of predicting ss one single sequence from scratch(generate alignments too) 
#use ab-initio predictor
open(SHELL, ">${bin_dir}predict_ssa_ab.sh") || die "can't create predict_ssa_ab.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "#predict ss for a single sequence from scratch using ab-initio neural network.\n";
print SHELL "if [ \$# -ne 2 ]\n";
print SHELL "then\n\techo \"need 2 parameters:seq_file(in fasta format), output_file.\" \n\texit 1\nfi\n";
print SHELL "#output: a file with predicted ss and sa.\n"; 
print SHELL "${script_dir}predict_ssa_ab.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_ss.sh $script_dir \$1 \$2 \n";
close SHELL; 

#generate script of predicting ss  and sa one single sequence from scratch(generate alignments too) 
#use homology ss and sa predictor
open(SHELL, ">${bin_dir}predict_ss_sa.sh") || die "can't create predict_ssa.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "#predict ss, sa for a single sequence from scratch.\n";
print SHELL "if [ \$# -ne 2 ]\n";
print SHELL "then\n\techo \"need 2 parameters:seq_file(in fasta format), output_file.\" \n\texit 1\nfi\n";
print SHELL "#output: a file with predicted ss and sa.\n"; 
print SHELL "${script_dir}predict_ss_sa.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_ss.sh ${server_dir}predict_seq_sa.sh $script_dir \$1 \$2 \n";
close SHELL; 

#generate script of predicting sa one single sequence from scratch(generate alignments too) 
#use homology sa predictor
open(SHELL, ">${bin_dir}predict_acc.sh") || die "can't create predict_acc.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "#predict sa for a single sequence from scratch.\n";
#print SHELL "if [ \$# -ne 3 ]\n";
#print SHELL "then\n\techo \"need 3 parameters:seq_file(in fasta format), output_file, threshold(0:0%,1:5%,...,5:25%,...19:95%).\" \n\texit 1\nfi\n";
print SHELL "if [ \$# -ne 2 ]\n";
print SHELL "then\n\techo \"need 2 parameters:seq_file(in fasta format), output_file\" \n\texit 1\nfi\n";
print SHELL "#output: a file with predicted ss and sa.\n"; 
#print SHELL "${script_dir}predict_acc.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_sa.sh $script_dir \$1 \$2 \$3\n";
print SHELL "${script_dir}predict_acc.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_sa.sh $script_dir \$1 \$2 \n";
close SHELL; 

#generate script of predicting sa one single sequence from scratch(generate alignments too) 
#use homology sa predictor
open(SHELL, ">${bin_dir}predict_acc_multi.sh") || die "can't create predict_acc_multi.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "#predict sa for a single sequence from scratch.\n";
print SHELL "if [ \$# -ne 3 ]\n";
print SHELL "then\n\techo \"need 3 parameters:seq_file(in fasta format), output_file, threshold(0:0%,1:5%,...,5:25%,...19:95%).\" \n\texit 1\nfi\n";
print SHELL "#output: a file with predicted ss and sa.\n"; 
print SHELL "${script_dir}predict_acc_multi.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_sa.sh $script_dir \$1 \$2 \$3\n";
close SHELL; 

#generate script of predicting sa one single sequence from scratch(generate alignments too) 
#use ab-initio sa predictor
open(SHELL, ">${bin_dir}predict_acc_ab.sh") || die "can't create predict_acc_ab.sh.\n";
print SHELL "#!/bin/sh\n";
print SHELL "#predict sa for a single sequence from scratch.\n";
#print SHELL "if [ \$# -ne 3 ]\n";
#print SHELL "then\n\techo \"need 3 parameters:seq_file(in fasta format), output_file, threshold(0:0%,1:5%,...,5:25%,...19:95%).\" \n\texit 1\nfi\n";
print SHELL "if [ \$# -ne 2 ]\n";
print SHELL "then\n\techo \"need 2 parameters:seq_file(in fasta format), output_file\" \n\texit 1\nfi\n";
print SHELL "#output: a file with predicted ss and sa.\n"; 
print SHELL "${script_dir}predict_acc_ab.pl $blast_dir $big_db $nr_db ${server_dir}predict_seq_sa.sh $script_dir \$1 \$2 \n";
close SHELL; 

`chmod 755 ${bin_dir}*.sh`; 
`chmod 755 ${server_dir}*`;
`chmod 755 ${script_dir}*`;
