#!/usr/bin/perl -w
use strict;

# THIS IS A MODIFIED SCRIPT DERIVED FROM HAKAN VIKLUND
# IT TAKES BLAST-M6-OUTPUT AS INPUT AND CONVERT IT TO A MSA
# IN FASTA-FORMAT

my $round   = 0;
my @no_hits;
    
my $msafile = $ARGV[0]; #$indir.$file;    
my $default = ' ;' x 60;

#for one blast_msa_in_file

my @seqs = ();
my $query_start = 0;
my $query_end = 0;
my $query = "";
my $fasta_length;
my $msa_length;
my $cur_block = 0;

my $fastafile = $ARGV[1]; 
if(-e $fastafile) {
    open(FASTAFILE, $fastafile) || die "could not open $fastafile\n";
}
else {
    print STDERR "could not open dir \n";
    next; 
}

<FASTAFILE>;
while(<FASTAFILE>) {
    my $line = $_;
    chomp($line);	
    last if ( $line =~ /^\>/ );
    $query .= $line;
}
$fasta_length = length $query;
#print "*$query*\n";
close(FASTAFILE);


if (!open(MSAFILE,$msafile)) {
    print STDERR "Could not open $msafile\n";
    next;
}

#check if file is ok
my $OK = -1;
while(<MSAFILE>) {
    if($_ =~ 'Query=') {
	$OK = 1;
    }
    if($_ =~ ' letters\)') {
	@_ = split;
	$msa_length = substr $_[0],1;
	$msa_length =~ s/,//g;
	last;
    }
}
if($OK < 0) {
    print STDERR "Error: msa_file $msafile is corrupt\n";
    next;
}

if($msa_length ne $fasta_length) {
    print STDERR "$msa_length $fasta_length\n";
    
    print STDERR "Error: non equal lengths of msa and fasta sequence\n";
    print STDERR "fasta_length = $fasta_length, msa_length = $msa_length \n";
    print STDERR "$fastafile\n";
    print STDERR "$msafile\n";
    print STDERR " $query\n";
    #exit;
    #next;
}

#find beginning of msa info
#print "$round\n";
while(<MSAFILE>) {
    chomp;
    if($round > 0) {
	if($_ =~ "Results from round $round" || $_ =~ 'CONVERGED!') {
	    #print "found results from round...\n";
	    print STDERR "$_\n";
	    last;
	}
	if($_ =~ 'No hits found') {
	    print STDERR "No hits found for $msafile\n";
	    $OK = -1;
	    push( @no_hits, "AA");
#		open(OUT,">".$outdir.$msaname.".mod");
#		$query =~ s/(.)/$+;/g;
#		printf OUT "<${query}>\n";
	    last;
	}
    }
    else {
	if($_ =~ 'No hits found') {
	    print STDERR "No hits found for $msafile\n";
	    $OK = -1;
	    push( @no_hits, "AA" );
#		open(OUT,">".$outdir.$msaname.".mod");
#		$query =~ s/(.)/$+;/g;
#		printf OUT "<${query}>\n";
	    last;
	}
	elsif($_ =~ 'Sequences producing significant alignments:') {
	    last;
	}
    }
}

if($OK < 0) {
    print STDERR "Error: strange format in msafile, could not find start of alignment info\n";
    next;
}


#find first sequence
while(<MSAFILE>) {
    chomp;
    #print "$_";
    if($_ =~ 'QUERY') {
	my $id;
	my $begin;
	my $end;
	my $seq;
	
	($id, $begin, $seq, $end) = split /\s+/, $_;
	$query_start = $begin;
	$query_end = $end;
	
	@seqs = (@seqs, $seq);
	last;
    }
}

my $seq_index = 1; 

while(<MSAFILE>) {
    if($_ eq "\n") {
	$cur_block++;
    }
    elsif($_ =~ 'Database') {
	last;
    }
    else {
	chomp;
	
	my $id;
	my $begin;
	my $end;
	my $seq;
	my @list = split /\s+/, $_;
	if($#list == 3) {
	    ($id, $begin, $seq, $end) = @list;
	}
	elsif($#list == 1 && substr($_,0,6) ne 'QUERY ') {
	    ($id, $seq) = @list;
	}
	elsif($#list == 2 && substr($_,0,6) ne 'QUERY ') {		
	    ($id, $begin, $seq) = @list;
	    if(length $seq < 50) {
		($id, $seq, $end) = @list;
	    }
	    if(length $seq < 50) {
		print STDERR "strange row in msafile: $_\n";
		exit;
	    }
	}
	else {
	    print STDERR "strange row in msafile: $_\n";
	    exit;		
	}
	
	if(substr($_,0,6) eq 'QUERY ') {
	    $query_end = $end;
	    $seqs[0] .= $seq;
	    $seq_index = 1;
	}
	else {
	    if($cur_block == 0) {
		@seqs = (@seqs, $seq);
	    }
	    else {
		if($seq_index > $#seqs) {
		    print STDERR "Error: msafile does not have the same number of aligned sequences in all blocks\n";
		    print STDERR "No output produced for AA\n";
		    exit;
		}
		#print "whole $_\n";
		#print "id $id\n";
		#print "begin $begin\n";
		    #print "end $end\n";
		#print "seq $seq\n";
		$seqs[$seq_index] .= $seq;
		$seq_index++;
	    }
	}
    }
}


#    print STDERR "$outfile\n";
#    open(OUT,">"."$outfile")
#	or die "could not open $outfile\n";
        
my $pos = $#seqs;
my $pream  = "";
my $restam = "";
if($query_start > 1) {
    
    printf ">AA\n";
    
    my $pre = substr($query,0,$query_start-1);
    #$pre =~ s/(.)/$+;/g;
    my @query_pre_list = split //, $pre;
    for(my $k = 0; $k <= $#query_pre_list; $k++) {
#	    printf OUT "$query_pre_list[$k]".";";
	printf "$query_pre_list[$k]";
    }
    my @query_seq_list = split //, $seqs[0];
    for(my $k = 0; $k <= $#query_seq_list; $k++) {
#	    printf OUT "$query_seq_list[$k]".";";
	    printf  "$query_seq_list[$k]";
    }
#	$pream = " ;" x ($query_start-1);
    $pream = "-" x ($query_start-1);
}
else {
    my @query_seq_list = split //, $seqs[0];
#	printf OUT "<";
    printf  ">AA\n";
    for(my $k = 0; $k <= $#query_seq_list; $k++) {
#	    printf OUT "$query_seq_list[$k]".";";
	printf "$query_seq_list[$k]";
    }
}
if($query_end < length($query)) {
    my $rest = substr($query,$query_end);
    #$rest =~ s/(.)/$+;/g;
    my @query_end_list = split //, $rest;
    for(my $k = 0; $k <= $#query_end_list; $k++) {
#	    printf OUT "$query_end_list[$k]".";";
	printf "$query_end_list[$k]";
    }
    $restam = "-" x ( length ($query) - $query_end);
}
#    printf OUT ">\n";
    printf  "\n";

for(my $j = 1; $j <= $pos; $j++) {
    printf ">$j\n";
    my @seq_list = split //, $seqs[$j];
    for(my $k = 0; $k <= $#seq_list; $k++) {
	if($seq_list[$k] ne '-') {
	    last;
	}
	else {
#		$seq_list[$k] = ' ';
	    $seq_list[$k] = '-';
	}
    }
    for(my $k = $#seq_list; $k >= 0; $k--) {
	    if($seq_list[$k] ne '-') {
		last;   #next?
	    }
	    else {
#		$seq_list[$k] = ' ';
		$seq_list[$k] = '-';
	    }
    }
#	printf OUT "<${pream}";
    printf "${pream}";
    for(my $k = 0; $k <= $#seq_list; $k++) {
#	    printf OUT "$seq_list[$k]".";";
	printf "$seq_list[$k]";
    }
#	printf OUT ">\n";
    
    printf "${restam}";
    printf "\n";
}

    #print STDERR "$counter Done modding $msaname\n";


#exit;

#close MSA_INDIR;


foreach my $line (@no_hits ){
    print "$line\n";
}
