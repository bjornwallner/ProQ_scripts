#!/usr/bin/perl -w

my ($seq,$topo)=parse_topcons($ARGV[0]);
print ">$ARGV[0]\n$topo\n";

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
    if(length($topcons)==0) {
	$topcons=$seq;
	$topcons=~s/./i/g;
    }
   # "FILE $file\n";
    #print "SEQ     $seq\n";
    #print "TOPCONS $topcons\n";
    close(FILE);
    return($seq,$topcons);
}


