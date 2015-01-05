#!/usr/bin/perl -w 

# extract subset of zpred by given the target sequence and model sequence
use strict;
use Cwd 'abs_path';
use File::Basename;
my $rundir = dirname(abs_path($0));
require "$rundir/needle.pl";

my $progname = basename($0);

my $usage="
Usage: $progname topfile seqfile model_seqfile model_topfile

Description:
    Create topology file (FASTA format) and the span file 
    for the model based on the topology file for the target sequence

Created 2014-05-05, updated 2014-05-05, Nanjiang Shu
";

my $numArgs = $#ARGV +1 ;

if ($numArgs < 4){
    print "$usage\n";
    exit;
}

my $topfile = $ARGV[0];
my $seqfile = $ARGV[1];
my $model_seqfile = $ARGV[2];
my $outfile = $ARGV[3];

my ($rootname_model_seqfile, $path_model_seqfile, $suffix_model_seqfile) = fileparse($model_seqfile, '\.[^\.]*');
my $basename_outfile = basename($outfile);

#`$rundir/acc_subset.pl $topfile $seqfile $model_seqfile $outfile`;
# 1. create subset for topcons.fa file
my $seq = ReadSeq($seqfile);
my $model_seq = ReadSeq($model_seqfile);
my $topo = ReadSeq($topfile);
my @arr_topo = split(//,$topo);

my ($a1,$a2)=align($seq,$model_seq);
my @a1=split(//,$a1);
my @a2=split(//,$a2);
open(OUT,">$outfile");
print OUT ">$basename_outfile\n";
my $resno=0;
foreach my $tp(@arr_topo) {
    if($a1[$resno] eq $a2[$resno]) {
        print OUT $tp;
    }
    $resno++;
}
print OUT "\n";
close(OUT);

#2. create subset for span file
my $model_topfile = $outfile;
my $model_spanfile = "$path_model_seqfile/$rootname_model_seqfile.span";
Create_subset_spanfile($model_topfile, $model_spanfile);

# extract span file
sub Create_subset_spanfile{#{{{
    my $topfile = shift;
    my $spanfile = shift;
    my $topo=ReadSeq($topfile);
    my $span = Mspan($topo);
    open(SPAN,">$spanfile");
    print SPAN "Prediction from TOPCONS on $ARGV[0]\n";
    print SPAN $span;
    close(SPAN);
}#}}}
sub Mspan {#{{{
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
    my $end=$i;
#    print "$label $start $end\n";
    if($label eq "M")
    {
        $span.=sprintf("%d %d %d %d\n",$start,$end,$start,$end);
        $counter++;
    }
    my $len=length($str);
    $span="$counter $len\nantiparallel\nn2c\n$span";
    return $span;
}#}}}
sub ReadSeq{#{{{
    my $file=shift;
    my $seq=`grep -v '>' $file`;
    $seq=~s/\n//g;
    return $seq;
}#}}}
