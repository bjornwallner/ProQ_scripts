#!/usr/bin/perl
###################################################
# getSSAligns
#  - finds alignments to AA sequence and outputs
#    predicted SS, using - for non-matches
#
#  Useage: getSSAligns.pl $temp_file $dir < seq_file
#
# Author: Michael Sweredoski
# Date: 1/7/04
#
# Modified by Jianlin cheng: change interface for integration with PSpro
# change name of temporary file, make it robust in multi-process environment. 
###################################################

$temp_file = shift @ARGV;
$install_dir = shift @ARGV;
$dir = $install_dir . "script"; 

$threshold = .82;
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

$aligns_file = "$data_dir/secondary_structure";
$index_file = "${temp_file}fooi.bai";

sub getSS {
    $name = shift @_;
    open(ALIGNS,"<$aligns_file") or die "Couldn't open aligns file $aligns_file\n";
    while($line = <ALIGNS>) {
	chomp($line);
	if($line eq $name) {
	    $ss = <ALIGNS>;
	    chomp($ss);
	    close(ALIGNS);
	    return $ss;
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
   # if($subj ne $query_name) {
	if($subj_SS = getSS($subj)) {
	    $subj_SS =~ s/\s+//g;	    
	    $subj_section = substr($subj_SS,$subj_start-1,$aligned+$gaps);
	    $pred_score = -1.108 + .15*log($aligned) + .359*log($ident);
	    if($pred_score > 1) {
		$pred_score = 1;
	    } elsif($pred_score < 0) {
		$pred_score = 0;
	    }

	    $subj_section =~ tr/GI/H/;
	    $subj_section =~ tr/B/E/;
	    $subj_section =~ tr/TS./C/;
	    
	    @subj_CHARS = split "",$subj_section;

	    push @aligns, {
		p_score => $pred_score,
		s_section => [@subj_CHARS],
		q_start => $query_start-1,
		q_len => $aligned+$gaps
	    };
	}
   # }
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
	print "-";
    }
}

print "\n";
`rm -f ${temp_file}foob.txt ${temp_file}fooi.bai $query_file`;




