#!/usr/bin/perl -w
#predict SA for a single sequence
#input parameters:
#1. sa predictor
#2. input sequence file(two lines: name(shouldn't include . or space in name), sequence(compact format, no space))
#3. alignment directory 
#Notice: invariance:  alignment_file = alignment_dir + name
#4. output file

#April 7, 2004, Jianlin Cheng

if (@ARGV != 4)
{
	die "need four parameters: SA_predictor, input_sequence, align_dir, output file\n";
}

$acc_predictor = $ARGV[0]; 
$seq_file = $ARGV[1];
$align_dir = $ARGV[2];
$output_file = $ARGV[3]; 
open(OUTPUT, ">$output_file") || die "can't open the output file.\n"; 
if (! -f $acc_predictor)
{
	die "can't find solvent acc predictor.\n"; 
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
#the alignment file name contains space or dot, it will cause trouble for ACCpro. 
if ($name =~ /[\. ]/)
{
	die "sequence name shouldn't contain dot or space.\n";
}

$sequence = shift @content;
chomp $sequence; 
$sequence =~ s/\s//g;

#create a temporary file and make prediction
open(TEMP, ">$seq_file.tmp") || die "can't create temporary file.\n";
print TEMP "1 20 3\n"; #create a title line required by ACCpro
print TEMP "$name\n"; 
print TEMP "$sequence"; 
close(TEMP);

#predict solvent accessibility with format 2 of ACCpro
$command = "$acc_predictor $seq_file.tmp $align_dir 2 > $seq_file.accpro"; 
system($command); 
open(RES, "$seq_file.accpro") || die "can't open SA prediction result file.\n"; 
$pre_acc = <RES>;
chomp $pre_acc; 
close(RES);

#remove the temporary files
`rm $seq_file.tmp $seq_file.accpro`;

$seq_length = length($pre_acc);
if ( $seq_length != length($sequence) )	
{
	print "pre: $pre_acc\n";
	print "seq: $sequence\n";
      die "sequence length doesn't match predicted solvent acc length.\n"; 
}

#convert buried flag "b" to "-"
$pre_acc =~ s/b/-/g; 
print OUTPUT "$name\n"; 
print OUTPUT "$sequence\n"; 
print OUTPUT "$pre_acc\n";
close OUTPUT; 

