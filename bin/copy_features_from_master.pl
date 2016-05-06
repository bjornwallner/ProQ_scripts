#!/usr/bin/perl -w
use Cwd 'abs_path';
use File::Basename;

my $INSTALL_DIR=dirname(abs_path($0));
if(scalar(@ARGV)==0) {
    
    print "\nThis script will copy sequence based features calculated for the full length sequence\n";
    print "USAGE:\n";
    print "\tcopy_features_from_master.pl <pdbfile/silent> <basename-template> <subunits optional>\n\n";
    exit;
}
my $subunits="";;

my $pdb=$ARGV[0];
my $basename=$ARGV[1];
$subunits=$ARGV[2] if(defined($ARGV[2]));
my $fasta="$pdb.fasta";
my $seq=`head $pdb|grep ^SEQUENCE|awk '{print \$2}'`;
if(length($seq)==0) { # assume file is a pdb...
    print STDERR "Assuming file is a pdb\n";
    $seq=`$INSTALL_DIR/aa321CA.pl $pdb`;
}

open(OUT,">$fasta");
print OUT ">$pdb\n$seq\n";
close(OUT);

`$INSTALL_DIR/acc_subset.pl $basename.acc $basename.fasta $fasta $pdb.acc $subunits`;
if(-e "$basename.mpSA") {
#    print "$INSTALL_DIR/mpSA_subset.pl $basename.mpSA $fasta $pdb.mpSA\n";
    #next;
    `$INSTALL_DIR/mpSA_subset.pl $basename.mpSA $fasta $pdb.mpSA`;
}
if(-e "$basename.zpred") {
    `$INSTALL_DIR/zpred_subset.pl $basename.zpred $fasta $pdb.zpred`;
}
if(-e "$basename.topcons.fa") {
    `$INSTALL_DIR/topcons_subset.pl $basename.topcons.fa $basename.fasta $fasta $pdb.topcons.fa`;
    my $spanfile="$pdb.span";
#    my ($seq2,$topo2)=parse_topcons("$pdb.topcons.fa");
    my $topo=`grep -v '>' $pdb.topcons.fa`;
    chomp($topo);
    my $span=Mspan($topo);
    open(SPAN,">$spanfile");
    print SPAN "Prediction from TOPCONS\n";
    print SPAN $span;
    close(SPAN);
}
`$INSTALL_DIR/profile_subset.pl $basename.psi $fasta $pdb.psi $subunits`;
`$INSTALL_DIR/profile_subset.pl $basename.mtx $fasta $pdb.mtx $subunits`;
`$INSTALL_DIR/ss2_subset.pl $basename.ss2 $fasta $pdb.ss2 $subunits`;


sub parse_topcons
{
    my $file=shift;
    my $topcons="";
    my $seq="";
    my $octopus="";
    my $get_prediction=0;
    my %pred=();
    $pred{'TOPCONS'}="";
    open(FILE,$file);
    while(<FILE>)
    {
	chomp;
	
	if($get_sequence) {
	    if(length($_)==0) {
		$get_sequence=0;
	    } else {
		$seq.=$_;
	    }
	}
	if($get_prediction) {
	    if(length($_)==0) {
		$get_prediction=0;
	    } else {
		$pred{$method}.=$_;
	    }
	}
	if(/(.+)\spredicted topology:/) {
	    $method=$1;
	    $get_prediction=1;
	}
	$get_sequence=1 if(/Sequence:/);


    }
    close(FILE);
    $seq=~s/\s+//g;
 #   print $seq."\n";
    foreach $method(keys(%pred)) {
	$pred{$method}=~s/\s+//g;
#	print "$method\n$pred{$method}\n";
    }
    $pred{'TOPCONS'}=$pred{'OCTOPUS'} if(length($pred{'TOPCONS'})==0);
    if(length($topcons)==0) {
	$topcons=$seq;
	$topcons=~s/./i/g;
    }
   # "FILE $file\n";
    #print "SEQ     $seq\n";
    #print "TOPCONS $topcons\n";

    return($seq,$pred{'TOPCONS'});
}


sub Mspan
{
    my $str=shift;
    my $span="";
    my @pred=split(//,$str);
    my $start=1;
    my $label=$pred[0];
    my $i=1;
    my $counter=0;
    for(;$i<scalar @pred;$i++)
    {
	if($pred[$i-1] ne $pred[$i])
	{
	    my $end=$i;
	#    print "$label $start $end\n";
	   
	    if($label eq "M")
	    {
		$span.=sprintf("%4d  %4d  %4d  %4d\n",$start,$end,$start,$end);
		$counter++;
	    }
	     $start=$i+1;
	    $label=$pred[$i];

	}
    }
    $end=$i;
#    print "$label $start $end\n";
    if($label eq "M")
    {
	$span.=sprintf("%d %d %d %d\n",$start,$end,$start,$end);
	$counter++;
    }
    my $len=length($str);
    $span="$counter $len\nantiparallel\nn2c\n$span";
    return $span;
}
