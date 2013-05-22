#!/usr/bin/perl -w
####################################
# getACCproPred
#  - reads an amino acid sequence from STD_IN 
#    and writes the SSpro prediction to STD_OUT
#  Useage: getSSproPred.pl $align_file $temp_file $dir < seq_file > output_preds
#
#  Author: Michael Sweredoski
#  Date: 1/7/04
# Modified by Jianlin Cheng: integrate with PSpro  ,use new sspro predictor
# Modified by Mike Sweredoski: adapted from getSSproPred
####################################

$align_file = shift @ARGV;
$temp_file = shift @ARGV;
$install_dir = shift @ARGV;
$acc_threshold = shift @ARGV;

<>;
$aa = <>;

#get align file name and dir
$pos = rindex($align_file, "/");
if ($pos >= 0)
{
	$ali_name = substr($align_file, $pos+1, length($align_file)-1 - $pos);
	$align_dir = substr($align_file, 0, $pos+1); 
}
else
{
	$ali_name = $align_file;
	$align_dir = "./"; 
}
# get SSpro prediction
open(AA,">$temp_file");
print AA "$ali_name\n$aa";
close(AA);
`$install_dir/bin/predict_seq_sa_multi.sh $temp_file $align_dir $temp_file.out $acc_threshold`;

# parse SSpro prediction
open(ACCpro,"<$temp_file.out");
<ACCpro>;
<ACCpro>;
$ACCpro_predict = <ACCpro>;
print "$ACCpro_predict";
`rm -fr $temp_file $temp_file.out`;

    
    
    
    
