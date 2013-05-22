#!/usr/bin/perl
###################################################
# getACCAligns
#  - finds alignments to AA sequence and outputs
#    predicted SS, using - for non-matches
#
#  Useage: getACCAligns.pl $temp_file $dir < seq_file
#
# Author: Michael Sweredoski
# Date: 1/7/04
#
# Modified by Jianlin cheng: change interface for integration with PSpro
# change name of temporary file, make it robust in multi-process environment. 
#
# Modified by Michael Sweredoski: adapted from getSSAligns
###################################################

$temp_file = shift @ARGV;
$install_dir = shift @ARGV;
$dir = $install_dir . "script"; 
$acc_threshold = shift @ARGV;
$acc_threshold *= 5; 

$threshold = .95;
#$threshold = .50;

$query_file = $temp_file;
$NCBI_PATH = "${install_dir}blast2.2.8";
$data_dir = "${install_dir}data/pdb_large";
#$data_dir = "${install_dir}data/pdb_small";
$query_header = <>;
$aa = <>;
open(TEMP_QUERY,">$query_file") or die "Couldn't open $temp_file for writing\n";
#add > to convert to fasta format
print TEMP_QUERY ">$query_header$aa";
close(TEMP_QUERY);

`$NCBI_PATH/blastall -i $query_file -d $data_dir/dataset -p blastp -F F -g F -o ${temp_file}foob.txt`;

`$dir/blast2index.pl ${temp_file}foob.txt > ${temp_file}fooi.bai`;

$aligns_file = "$data_dir/acc_data";
$index_file = "${temp_file}fooi.bai";

sub getACC {
    $name = shift @_;
    open(ALIGNS,"<$aligns_file") or die "Couldn't open aligns file $aligns_file\n";
    while($line = <ALIGNS>) {
	chomp($line);
	if($line eq $name) {
	    $acc = <ALIGNS>;
	    chomp($acc);
	    my @temp = split /\s+/,$acc;
	    $acc = "";
	    for(my $i = 0; $i <= $#temp; $i++) {
		$acc .= ($temp[$i] > $acc_threshold?"e":"-");
	    }
	    close(ALIGNS);
	    return $acc;
	} else {
	    <ALIGNS>;
	}
    }
    close(ALIGNS);
    return "";
}

open(INDEX,$index_file) or die "Couldn't open index file\n";

$line = <INDEX>;
($query_name,$query_len) = split " ",$line;

@aligns = ();

while($line = <INDEX>) {
    ($aligned,$gaps,$eval,$ident,$posi,$subj,$subj_len,$subj_start,$subj_end,$query_start,$query_end) = split "\t",$line;
#    if($subj ne $query_name) {
	if($subj_ACC = getACC($subj)) {
	    $subj_ACC =~ s/\s+//g;	    
	    $subj_section = substr($subj_ACC,$subj_start-1,$aligned+$gaps);

	    $pred_score = 0.3230 + .0260*log($aligned) + .1078*log($ident);

	    if($pred_score > 1) {
		$pred_score = 1;
	    } elsif($pred_score < 0) {
		$pred_score = 0;
	    }

	    @subj_CHARS = split "",$subj_section;

	    push @aligns, {
		p_score => $pred_score,
		s_section => [@subj_CHARS],
		q_start => $query_start-1,
		q_len => $aligned+$gaps
	    };
	}
#    }
}

@sorted_i = sort { $aligns[$b]{"p_score"} <=> $aligns[$a]{"p_score"} } (0..$#aligns);
@sorted_aligns = ();
for($i = 0; $i <= $#aligns; $i++) {
    $sorted_aligns[$sorted_i[$i]] = $aligns[$i];
} 

for($i = 0; $i < $query_len; $i++) {
    $found = 0;
    for($j = 0; $j <= $#sorted_aligns; $j++) {
	if($sorted_aligns[$j]{"q_start"} <= $i && $sorted_aligns[$j]{"q_len"}+$sorted_aligns[$j]{"q_start"} > $i && $sorted_aligns[$j]{"p_score"} > $threshold) {
	    if($found == 0) {
		print $sorted_aligns[$j]{"s_section"}[$i-$sorted_aligns[$j]{"q_start"}];
		$found = 1;
	    }
	}
    }
    if($found == 0) {
	print "n";
    }
}

print "\n";
`rm -f ${temp_file}foob.txt ${temp_file}fooi.bai $query_file`;




