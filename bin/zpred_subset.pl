#!/usr/bin/perl -w 

# extract subset of zpred by given the target sequence and model sequence
use strict;
use Cwd 'abs_path';
use File::Basename;
my $rundir = dirname(abs_path($0));
require "$rundir/needle.pl";

my $progname = basename($0);

my $usage="
Usage: $progname zpredfile model_fastafile model_zpredfile

Description:
    Create zpred file for the model based on the zpred file for the target
    sequence

Created 2014-05-05, updated 2014-05-05, Nanjiang Shu
";

my $numArgs = $#ARGV +1 ;

if ($numArgs < 3){
    print "$usage\n";
    exit;
}

my $zpredfile = $ARGV[0];
my $model_seqfile = $ARGV[1];
my $outfile = $ARGV[2];

my $seq=`grep -v '>' $model_seqfile`;
$seq=~s/\n//;

my ($profile_seq, $zpred)=read_in_zpred_all($zpredfile);
my @zpred = @{$zpred};

my ($a1,$a2) = align($profile_seq, $seq);

#print "P: $a1\nS: $a2\n";
my @a1=split(//,$a1);
my @a2=split(//,$a2);

open(OUT,">$outfile");
my $resno=0;
foreach my $line(@zpred) {
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

sub read_in_zpred_all {#{{{
    my $file=shift;
    my $seq="";
    my @lines = "";
    open(FILE,"$file");# or die "Cannot open $file.\n";
    while(<FILE>)
    {
        if(/\d+/) {
            chomp;
            my @temp=split(/\s+/);
            $seq .= $temp[-1];
            push(@lines,"$_\n");
        }
    }
    close(FILE);
    return ($seq,[@lines]);

}#}}}
