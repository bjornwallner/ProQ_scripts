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
$seq=~s/\n//g;
#check for oligomer

#print "'".$seq."'"."\n";
$mtxfile=0;
$ext=".psi";
$count=0;
$header="";
$tail="";
$profile_seq="";
@profile=();

open(FILE,"$ss2");
while(<FILE>) {
    if($mtxfile) {
	if($count==1) {
	    chomp;
	    $profile_seq=$_;
	} elsif($count>=2 &&
		$count<=13) {
	    $header.=$_;
	} else {
	    
	    push(@profile,$_);
	}
	


    } 
    else
    {
	if($count>=1 && $count <= 2) {

	    $header.=$_;
	} 
	elsif(/\d+\s([A-Z])\s/)
	{
	    $profile_seq.=$1;
	    push(@profile,$_);

	}
	elsif($count>1) {
	    $tail.=$_;
	}
	
    }
    if($count==0 && /^\d+$/) {
	$mtxfile=1;
	$ext=".mtx";

    }
    $count++;
}

$len_seq=length($seq);
$len_profile=length($profile_seq);
#$subunits=ceil($len_seq/$len_profile);
print "$len_seq $len_profile $subunits subunits\n";
#exit;

my $temp_seq="";
my @temp_profile=();
for(my $i=0;$i<$subunits;$i++) {
    $temp_seq.=$profile_seq;
    @temp_profile=(@temp_profile,@profile);
}

@profile=@temp_profile;
$profile_seq=$temp_seq;
#print $mtxfile."\n";
#print $header;
#print $profile_seq."\n";
#$s=scalar(@profile);
#print $s."\n";
($a1,$a2)=align($profile_seq,$seq);

print "P: $a1\nS: $a2\n";
#exit;
@a1=split(//,$a1);
@a2=split(//,$a2);

#$outfile=$seqfile;
#$outfile=~s/\.seq/$ext/g;
#print $outfile."\n";

if($mtxfile) {
    $len=length($seq);
    $header="$len\n$seq\n$header";

}
else 
{
    $header="\n$header";
}
open(OUT,">$outfile");
print OUT $header;
$resno=0;

foreach $profile(@profile) {
    if($a1[$resno] eq $a2[$resno]) {
	    
	print OUT $profile;
    }
    $resno++;


}
print OUT $tail;
#print OUT "\n";


