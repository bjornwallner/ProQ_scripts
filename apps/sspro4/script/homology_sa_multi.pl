#!/usr/bin/perl -w
###################################################
# homology_sa.pl
#  - makes a ACC prediction for a protein using
#    ACCpro and homology
#  Useage: homology.pl install_dir seq_file(name, sequence(in compact format)) align_file  output_file
#
# Author: Michael Sweredoski
# Date: 1/8/04
#
# Modified by Jianlin Cheng, 5/1/2004
# Integration with PSpro, change inteface, add error checking
# Output(chaged): name, seq, predicted_ss
#
# Modified by Mike Sweredoski: adapted from homology_ss
###################################################

if (@ARGV != 5)
{
	die "need 5 params: install dir, seq file(name, sequence in compact format), alignment file, output file acc_threshold\n"; 
}

$install_dir = shift @ARGV; 
$seq_file = shift @ARGV;
$align_file = shift @ARGV;
$output_file = shift @ARGV;
$acc_threshold = shift @ARGV;


if (! -d $install_dir)
{
	die "can't find install directory.\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
	$install_dir .= "/"; 
}
if (! -f $seq_file)
{
	die "input file doesn't exist.\n";
}
if (! -f $align_file)
{
	warn "alignment file doesn't exist.\n"; 
}
$temp_file = "${seq_file}.homo" . $$;
$dir = "${install_dir}script";

open(SEQ_FILE, "$seq_file") || die " can't open seq file.\n"; 
$seq_name = <SEQ_FILE>;
$seq_seq = <SEQ_FILE>; 
close SEQ_FILE; 

open(OUTPUT,">$output_file");

$accpro_pred = `cat $seq_file | $dir/getACCproPred_multi.pl $align_file $temp_file $install_dir $acc_threshold`;
$acc_homology_pred = `cat $seq_file | $dir/getACCAligns_multi.pl $temp_file $install_dir $acc_threshold`;

chomp($accpro_pred);
chomp($acc_homology_pred);

@ACC_pred = split "", $accpro_pred;
@ACC_HOM_pred = split "", $acc_homology_pred;
print OUTPUT $seq_name;
print OUTPUT $seq_seq; 
for($i = 0; $i <= $#ACC_pred; $i++) {
    if($ACC_HOM_pred[$i] ne "n") {
	print OUTPUT $ACC_HOM_pred[$i];
    } else {
	print OUTPUT $ACC_pred[$i];
    }
}
print OUTPUT "\n";
close OUTPUT; 
