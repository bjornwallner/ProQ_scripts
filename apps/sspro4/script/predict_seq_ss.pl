#!/usr/bin/perl -w

#predict the secondary structure for a single sequence
#April 5, 2004, Jianlin Cheng

#input sequence file: 2-compacted line(no space separator): name, sequence. 
#parameters:
#1. ss_predictor_path
#2. sequence_file
#3. alignment path (notice******: alignment_file = alignment_path + seq_name, sequence name should not contain "." or " ", this could cause problem, we may fix it by changing sspro code later)
#4. output file


if (@ARGV != 4)
{
	die "need four parameters: ss_predictor, input_sequence, align_dir, output_file\n";
}

$ss_predictor = $ARGV[0];
$seq_file = $ARGV[1];
$align_dir = $ARGV[2];
$output_file = $ARGV[3]; 
open(OUTPUT, ">$output_file") || die "can't open the output file.\n"; 

if (! -f $ss_predictor)
{
	die "can't find secondary structure predictor.\n"; 
}

if (! -d $align_dir)
{
	print "the alignment directory doesn't exists.\n"; 
}

if ( substr($align_dir, length($align_dir) - 1, 1) ne "/" )
{
	$align_dir .= "/"; 
}

if (! -f $seq_file)
{
	die "can't find the sequence file.\n"; 
}
open(SEQ_FILE, "$seq_file") || die "can't open sequence file.\n";
@content = <SEQ_FILE>;
close(SEQ_FILE);

$name = shift @content;
chomp $name; 
$name =~ s/\s//g;

#check if the alignment file exists
$alignment_file = $align_dir . $name; 
if (! -f $alignment_file)
{
	print "Warning: alignment file: $alignment_file doesn't exist!\n"; 
}
#the alignment file name contains space or dot, it will cause trouble. 
if ($name =~ /[\. ]/)
{
	die "sequence name shouldn't contain dot or space.\n";
}

$sequence = shift @content;
chomp $sequence; 
$sequence =~ s/\s//g;

#make up a pseudo ss string required by SSpro
$sstructure = ""; 
for ($ss_i = 0; $ss_i < length($sequence); $ss_i++)
{
	$sstructure .= "H"; 
}

open(TMPFILE, ">$seq_file.tmp") || die "can't create temporary file.\n"; 
print TMPFILE "1 20 3\n$name\n$sequence\n$sstructure"; 
close TMPFILE; 

#predict the secondary structure with 3-line format
system("$ss_predictor $seq_file.tmp $align_dir 0 > $seq_file.sspro"); 
open(RES, "$seq_file.sspro") || die "can't open prediction result file.\n"; 
$org_ss = <RES>;
chomp $org_ss;
$pre_ss = <RES>;
chomp $pre_ss; 
close(RES);
$seq_length = length($org_ss);
if ( $seq_length != length($sequence) || $seq_length != length($pre_ss) )	
{
	die "sequence length doesn't match with output from SSpro.\n"; 
}

#remove the temporary files
`rm $seq_file.sspro $seq_file.tmp`;

#output the results
print OUTPUT "$name\n"; 
print OUTPUT "$sequence\n"; 
print OUTPUT "$pre_ss\n"; 
close OUTPUT; 
