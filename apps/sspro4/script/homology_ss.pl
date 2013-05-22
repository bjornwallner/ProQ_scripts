#!/usr/bin/perl -w
###################################################
# homology.pl
#  - makes a SS prediction for a protein using
#    SSpro and homology
#  Useage: homology.pl install_dir seq_file(name, sequence(in compact format)) align_file  output_file
#
# Author: Michael Sweredoski
# Date: 1/8/04
#
# Modified by Jianlin Cheng, 5/1/2004
# Integration with PSpro, change inteface, add error checking
# Output(chaged): name, seq, predicted_ss
#
###################################################

if (@ARGV != 4)
{
	die "need 4 params: install dir, seq file(name, sequence in compact format), alignment file, output file\n"; 
}

$install_dir = shift @ARGV; 
$seq_file = shift @ARGV;
$align_file = shift @ARGV;
$output_file = shift @ARGV;


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

$sspro_pred = `cat $seq_file | $dir/getSSproPred.pl $align_file $temp_file $install_dir`;
$ss_homology_pred = `cat $seq_file | $dir/getSSAligns.pl $temp_file $install_dir`;

chomp($sspro_pred);
chomp($ss_homology_pred);

@SS_pred = split "", $sspro_pred;
@SS_HOM_pred = split "", $ss_homology_pred;
print OUTPUT $seq_name;
print OUTPUT $seq_seq; 
for($i = 0; $i <= $#SS_pred; $i++) {
    if($SS_HOM_pred[$i] ne "-") {
	print OUTPUT $SS_HOM_pred[$i];
    } else {
	print OUTPUT $SS_pred[$i];
    }
}
print OUTPUT "\n";
close OUTPUT; 
