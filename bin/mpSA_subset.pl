#!/usr/bin/perl -w 

# extract subset of mpsa by given the target sequence and model sequence
use strict;
use Cwd 'abs_path';
use File::Basename;
my $rundir = dirname(abs_path($0));
require "$rundir/bjornlib.pl";

my $progname = basename($0);

my $usage="
Usage: $progname mpSAfile model_fastafile model_mpSAfile

Description:
    Create mpSA file for the model based on the mpSA file for the target
    sequence

Created 2014-05-05, updated 2014-05-05, Nanjiang Shu
";

my $numArgs = $#ARGV +1 ;

if ($numArgs < 3){
    print "$usage\n";
    exit;
}

my $mpsafile = $ARGV[0];
my $model_seqfile = $ARGV[1];
my $outfile = $ARGV[2];

my $seq=`grep -v '>' $model_seqfile`;
$seq=~s/\n//;

my ($profile_seq, $mpsa)=read_in_mpsa_all($mpsafile);
my @mpsa = @{$mpsa};

my ($a1,$a2) = align($profile_seq, $seq);

#print "P: $a1\nS: $a2\n";
my @a1=split(//,$a1);
my @a2=split(//,$a2);

open(OUT,">$outfile");
my $resno=0;
foreach my $line(@mpsa) {
    if($line=~/\d+/) {
        if($a1[$resno] eq $a2[$resno]) {
            print OUT $line;
        }
        $resno++;
    } else {
        #print OUT;
    }
}
close(OUT);

sub read_in_mpsa_all {#{{{
    my $file=shift;
    my $seq="";
    my @lines = "";
    open(FILE,"$file");# or die "Cannot open $file.\n";
    while(<FILE>)
    {
        if(/\d+/) {
            chomp;
            my @temp=split(/\s+/);
            $seq .= $temp[0];
            push(@lines,"$_\n");
        }
    }
    close(FILE);
    return ($seq,[@lines]);

}#}}}
