#!/usr/bin/perl -w 

my ($seq,$topo)=parse_topcons($ARGV[0]);
my $span=Mspan($topo);

print "Prediction from TOPCONS on $ARGV[0]\n";
print $span;


sub parse_topcons
{
    my $file=shift;
    my $topcons="";
    my $seq="";
    my $octopus="";
    my $get_prediction=0;
    open(FILE,$file);
    while(<FILE>)
    {
	chomp;
	if($get_prediction) {
	    if(/^\s+TOPCONS\s+([iMo\s]+)/)
	    {
	#	print;
		#    print "\n";
		$topcons.=$1;
		
	    }
	    if(/^\s+OCTOPUS\s+([iMo\s]+)/)
	    {
		#    print;
		#    print "\n";
		$octopus.=$1;
		
	    }
	    if(/^\s+Seq\.\s+([A-Z\s]+)/) {
		
		#   print;
		#   print "\n";
		$seq.=$1;
	    }
	
	}  
	$get_prediction=1 if(/Sequence and predicted topologies:/);
    }
    $seq=~s/\s+//g;
    $topcons=~s/\s+//g;
    $octopus=~s/\s+//g;
    $topcons=$octopus if(length($topcons)==0);
   # print "FILE $file\n";
    #print "SEQ     $seq\n";
    #print "TOPCONS $topcons\n";
    close(FILE);
    return($seq,$topcons);
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



