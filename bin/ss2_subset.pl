#!/usr/bin/perl -w 
use Cwd 'abs_path';
use File::Basename;
my $rundir = dirname(abs_path($0));
require "$rundir/needle.pl";

$subunits=1;
$ss2=$ARGV[0];
$seqfile=$ARGV[1];
$outfile=$ARGV[2];
$subunits=$ARGV[3] if(defined($ARGV[3]));
$seq=`grep -v '>' $seqfile`;
$seq=~s/\n//;

($profile_seq,$ss)=read_in_psipred_all($ss2);
@ss=@{$ss};

my $temp_seq="";
my @temp_ss=();
for(my $i=0;$i<$subunits;$i++) {
    $temp_seq.=$profile_seq;
    @temp_ss=(@temp_ss,@ss);
}
$profile_seq=$temp_seq;
@ss=@temp_ss;
#$ss_pred;
#print $seq."\n";
#print $profile_seq."\n";
#print $ss."\n";
($a1,$a2)=align($profile_seq,$seq);

#print "P: $a1\nS: $a2\n";
@a1=split(//,$a1);
@a2=split(//,$a2);

#$outfile=$seqfile;
#$outfile=~s/\.seq/\.ss2/g;
#print $outfile."\n";
open(OUT,">$outfile");
print OUT "# PSIPRED VFORMAT (PSIPRED V2.5 by David Jones) ss2_subset.pl reformat\n\n";
$resno=0;
#open(FILE,$ss2);
#while(<FILE>) 
foreach $line(@ss)
{
    if($line=~/\d+\s\w\s\w\s/)
    {
	if($a1[$resno] eq $a2[$resno]) {
	    
	    print OUT $line;
	}
	$resno++;
    } else {
	print OUT;
	
    }
    
}
close(OUT);


sub read_in_psipred_all
{
    my $file=shift;
    #my $target_seq=shift;
    my $seq="";
    my $ss_pred="";
    open(FILE,"$file");# or die "Cannot open $file.\n";
    while(<FILE>)
    {
	if(/^Pred:\s[HCE]+/)
	{
	    my @temp=split(/\s+/);
	    $ss_pred.=$temp[1] if(scalar @temp >=2);
	}
	elsif(/^  AA:/)
	{
	    my @temp=split(/\s+/);
	    $seq.=$temp[2] if(defined($temp[2]));
	}
    }
    close(FILE);
    #print length($seq)," ",length($ss_pred),"\n";
    #print $file,"\n";
    if(length($seq)==0 && length($ss_pred)==0)
    {
	#print $file,"\n";
	open(FILE,"$file");# or die "Cannot open $file.\n";
	while(<FILE>)
	{
	    if(/\d+\s\w\s\w\s+/) {
		chomp;
		@temp=split(/\s+/);
		$ss_pred.=$temp[3];
		$seq.=$temp[2];
		push(@lines,"$_\n");
	    }
	}
    }
    close(FILE);
    #print "$seq\n$ss_pred\n";
    #exit;
    #generate_top_code($ss_pred,10,110);
#    return ($seq,$ss_pred);
    return ($seq,[@lines]);

}
