#!/usr/bin/perl -w

#predict SS, SA for a single sequence, SS, SA is a combination of neural network and homology

#June 23, 2004, Jianlin Cheng

#input 
	# 1: blast_dir
	# 2: big_db
	# 3: nr_db
	# 4: ss_predictor
	# 5: sa_predictor
	# 6: script_dir
	# 7: sequence file(seq_file): one sequence in fasta format, 
	# 8: output file

#Output: pxml, casp, eva format contact map for both 12 and 8 A. 
#Notice: the white space or . in seq_file are replaced by _. 

if (@ARGV != 8)
{
	die "Usage: blast_dir big_db nr_db ss_predictor sa_predictor script_dir seq_file(fasta) output_file\n";
}

$blast_dir = shift @ARGV;
$big_db = shift @ARGV;
$nr_db = shift @ARGV;
$ss_predictor = shift @ARGV; 
$sa_predictor = shift @ARGV; 
$script_dir = shift @ARGV; 

################verify if all the things are there#########################
if (! -d $blast_dir)
{
	die "can't find blast directory.\n";
}
if ( substr($blast_dir, length($blast_dir) - 1, 1) ne "/" )
{
	$blast_dir .= "/"; 
}

if (! -d $script_dir)
{
	die "can't find perl script directory.\n";
}

if ( substr($script_dir, length($script_dir) - 1, 1) ne "/" )
{
	$script_dir .= "/"; 
}

if (! -f "${blast_dir}blastpgp")
{
	die "can't find blastpgp.\n"; 
}

if (! -f "${big_db}.pal") {
if (! -f "${big_db}.phr" || ! -f "${big_db}.pin" || ! -f "${big_db}.psq")
{
	die "can't find the big coil database.\n"; 	
} 
}

if (! -f "${nr_db}.pal") {
if (! -f "${nr_db}.phr" || ! -f "${nr_db}.pin" || ! -f "${nr_db}.psq" || ! -f "${nr_db}.pnd"
    || ! -f "${nr_db}.pni" || ! -f "${nr_db}.psd" || ! -f "${nr_db}.psi" )
{
	die "can't find the non-redundant database.\n"; 	
} 
}

if (! -f $ss_predictor)
{
	die "can't find secondary structure predictor.\n"; 
}

if (! -f $sa_predictor)
{
	die "can't find solvent accessibility predictor.\n"; 
}
#############################End of Verification#######################################



$seq_file = shift @ARGV;
$output_file = shift @ARGV; 
$output_dir = "./"; 

if (! -f $seq_file)
{
	die "can't find file: $seq_file.\n"; 
}

if (! -d $output_dir)
{
	die "the output directory doesn't exists.\n"; 
}

if ( substr($output_dir, length($output_dir) - 1, 1) ne "/" )
{
	$output_dir .= "/"; 
}

#extract sequence file name
$slash_pos = rindex($seq_file, "/");
if ($slash_pos != -1)
{
	$seq_filename = substr($seq_file, $slash_pos + 1, length($seq_file) - $slash_pos - 1); 
}
else
{
	$seq_filename = $seq_file; 
}
if (length($seq_filename) <= 0)
{
	die "sequence file name shouldn't be less or equal 0.\n"; 
}

#non-char and . is not allowed for ouput file name 
$seq_filename =~ s/\s/_/g; 
$seq_filename =~ s/\./_/g;  

#output prefix is used as the prefix name of output files
$output_prefix = $output_dir . $seq_filename; 
#set alignment file name
$seq_filename .= "alg";
$output_prefix_alg = $output_dir . $seq_filename; 

open(SEQ_FILE, "$seq_file") || die "can't open sequence file.\n";
@content = <SEQ_FILE>;
close(SEQ_FILE);
$name = shift @content;
chomp $name; 
$sequence = shift @content; 

#remove unseen dos format (cause the program fail)
$name =~ s/\s//g;
$sequence =~ s/\s//g;

#check the sequence format
if (substr($name, 0, 1) ne ">") 
{
	die "sequence file: $seq_file is not in fasta format.\n"; 
}
$name = substr($name, 1, length($name) - 1); 
$target_name = $name; 
if (length($target_name) == 0)
{
	$target_name = "unknown"; 
}
if (length($sequence) < 1)
{
	die "seqeunce is empty. \n"; 
}
$script_dir =~ /(.+)script\/$/;
$install_dir = $1;

#Generate Alignments for sequence
print "Generate alignment...\n";
#the output alignment file name is: $output_prefix
`${script_dir}generate_flatblast.pl $blast_dir $script_dir $big_db $nr_db $seq_file $output_prefix_alg >/dev/null`;  
print "done\n";

#Predict Secondary Stx
print "Predict secondary structure...";
#create input file for sspro
open(TMP, ">$output_prefix.tmp") || die "can't create temporary file for sspro.\n"; 
#a little bit ugly here: for sspro: alignment_file = align_dir + seq_filename
print TMP "$seq_filename\n";
print TMP "$sequence\n"; 
close(TMP); 
`${script_dir}homology_ss.pl $install_dir $output_prefix.tmp $output_prefix_alg  $output_prefix.sspro`; 
print "done\n";

#Predict Solvent Accessibility
print "Predict solvent accessibility...";
`${script_dir}homology_sa.pl $install_dir $output_prefix.tmp $output_prefix_alg  $output_prefix.accpro`;
print "done\n"; 

open(TMP, ">$output_file") || die "can't create temporary file for conpro.\n"; 
open(SSPRO, "$output_prefix.sspro") || die "can't open the ss result file.\n";
open(ACCPRO, "$output_prefix.accpro") || die "can't open the sa result file.\n";
@sspro = <SSPRO>;
@accpro = <ACCPRO>;
$name = shift @sspro;
$seq = shift @sspro;
$sstx = shift @sspro;
shift @accpro;
shift @accpro;
$acc = shift @accpro;
$acc =~ s/-/b/g; 
print TMP "$seq_filename\n$seq$sstx$acc";
close TMP;

#rename the alignment file
`mv $output_prefix_alg ${output_file}align`;
#remove pssm file
#`rm $output_prefix_alg.pssm`;
#remove the intermediate files
`rm $output_prefix.accpro $output_prefix.sspro`; 
`rm $output_prefix.tmp`; 
