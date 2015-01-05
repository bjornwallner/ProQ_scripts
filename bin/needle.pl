#!/usr/bin/perl -w 
use File::Temp qw/ tempfile /;

sub align   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    my $needle_linux="/usr/bin/needle";
    my $needle_mac="/opt/local/bin/needle";
    my $osname = $^O;
    my $input1=$seq1;
    my $input2=$seq2;
    $seq1=~s/-//g;
    $seq2=~s/-//g;
    #$seq1=remove_dashes($seq1);
    #$seq2=remove_dashes($seq2);
    $seq1=~s/\n//g;
    $seq2=~s/\n//g;
    $seq1=~s/\s+//g;
    $seq2=~s/\s+//g;
    
    my ($fh1,$file1)=tempfile("/tmp/seq.XXXXXXXX");
    my ($fh2,$file2)=tempfile("/tmp/seq.XXXXXXXX");
    my ($fh3,$file3)=tempfile("/tmp/ali.XXXXXXXX");
    close($fh3);
    print $fh1 ">seq1\n$seq1\n";
    close($fh1);
    print $fh2 ">seq2\n$seq2\n";
    close($fh2);

    if($osname eq "linux" && -e $needle_linux) {
	$needle=$needle_linux;
    }
    if($osname eq "darwin" && -e $needle_mac) {
	$needle=$needle_mac;
    }

    #print "needle -aseq $file1 -bseq $file2 -gapopen 10 -gapextend 0.5 -outfile $file3\n";
    `needle -aseq $file1 -bseq $file2 -gapopen 1 -gapextend 0.5 -outfile $file3 > /dev/null 2>&1`;
    #print $file3."\n";
    ($ali_return1,$ali_return2)=parse_needle_output($file3);
    `rm $file1 $file2 $file3`;
   # if(length($seq1)==1 || length($seq2)==1)
   # {
   #	my $len1=length($seq1);
   #	my $len2=length($seq2);
   #	my $len3=length($input1);
   #	my $len4=length($input2);
   #	if($len1==1)
   #	{
   #	   $dashes=length($input2)-length($input1);
   #	   $ali_return1=$input1._dashes($dashes);
   #	   $ali_return2=$input2;
   #	}
   #	if($len2==1)
   #	{
   #	   $dashes=length($input1)-length($input2);
   #	   $ali_return1=$input1;
   #	   $ali_return2=$input2._dashes($dashes);
   #	}
   #	return ($ali_return1,$ali_return2)
   # }
   #
   #
   # else
    #print "$ali_return1\n$ali_return2\n";
   
    return ($ali_return1,$ali_return2);
    
}

sub parse_needle_output
{
    my $file=shift;
    my $seq1="";
    my $seq2="";
    my $header1="";
    my $header2="";
    
    my $isFirst=1;
    open(FILE,$file);
    while(<FILE>){
	next if (/^#/);
        if(/^(.{13})(.{6}\d)\ ([\w\-]+)/){
	  #  print "header:$1, seqnro:$2, seq:$3|\n";
	    my $header = $1;
	    my $seq = $3;
	    if ($isFirst){
		$seq1.=$seq;
		$header1 = $header;
		$isFirst=0;
	    }
	    else {
		$seq2.=$seq;
		$header2 = $header;
		$isFirst=1;
	    }
        }
	
    }
    close(FILE);
    if(length($seq1) == 0) {

	print STDERR "needle from the EMBOSS package needs to be installed.\n";
    }
    return($seq1,$seq2);
}
1;

