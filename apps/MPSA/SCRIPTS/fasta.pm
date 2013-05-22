#!/usr/bin/perl

use warnings;
use strict;

#=====================================
sub get_first_seq_from_fasta_seq_file{
    my ( $file ) = @_;
    open( ALI_FILE , "$file") || die "cannot open $file \n";
    my $found = 0 ;
    my $Seq   = '';
    foreach my $line ( <ALI_FILE> ) {
	if ($line =~ /^>/ ) {
	    if ($found eq 1) {
		last;
	    }
	    $found = 1;
	    $Seq ='';
	    next;
	}
	if ($found eq 1) {
	    if ($line =~ /([\w\-\*\.]+)/ ) {	
		$Seq .= $1;
	    }
	}
    }
    $Seq =~s /\s+//g;
    close ALI_FILE;
    return $Seq;
}
#=====================================
sub get_seq_from_fasta_seq_file_with_known_ID{
    my ( $file, $Id ) = @_;
    open( ALI_FILE , "$file") || die "cannot open $file \n";
    my $found = 0 ;
    my $Seq   = '';
    foreach my $line ( <ALI_FILE> ) {
	if ($line =~ /^>/ ) {
	    $found = 0;
	}
	if ($line =~ /^>$Id[\s+\n]/ ) {
#	    print "$line \n";
	    if ($found eq 1 ){
		#and $Id eq $1) {
#		print $line;
		last;
	    }
	    $found = 1;
	    $Seq ='';
	    next;
	}
	if ($found eq 1) {
	    if ($line =~ /([\w\-\*\.]+)/ ) {	
		$Seq .= $1;
	    }
	}
    }
    close ALI_FILE;
    return $Seq;
}
#=====================================
sub get_seq_from_fasta_seq_file_with_known_ID2{
    my ( $file, $Id ) = @_;
    open( ALI_FILE , "$file") || die "cannot open $file \n";
    my $found = 0 ;
    my $Seq   = '';
    foreach my $line ( <ALI_FILE> ) {
	if ($line =~ /^>/ ) {
	    $found = 0;
	}
	if ($line =~ /^>*$Id*/ ) {
#	    print "$line \n";
	    if ($found eq 1 ){
		#and $Id eq $1) {
#		print $line;
		last;
	    }
	    $found = 1;
	    $Seq ='';
	    next;
	}
	if ($found eq 1) {
	    if ($line =~ /([\w\-\*\.]+)/ ) {	
		$Seq .= $1;
	    }
	}
    }
    close ALI_FILE;
    return $Seq;
}
#=====================================
sub get_all_descriptions_from_fasta_seq_file{
    my ( $file ) = @_;
    my @descriptions;
    open( ALI_FILE , "$file") || die "cannot open $file \n";
    foreach my $line ( <ALI_FILE> ) {
#	if ($line =~ /^>([\w\-]+)\s.+/ ) {
	if ($line =~ /^>([\w\-\/]+)/ ) {
#	    print "$1\n";
	    push( @descriptions, $1 );
	}
    }
    close ALI_FILE;
    return @descriptions;
}
#=====================================
sub print_fasta{
    my ($file, $descr, $seq ) = @_;
    open( FILE, ">$file" ) || die "cannot open file1 \n";
    print FILE ">$descr\n";  
    my $temp_line ='';
    my $line_pos = 0;
    for( my $i=0; $i< length( $seq); $i++){
	my $chr = substr( $seq, $i, 1 );
	$temp_line .= "$chr";
	if ( $line_pos eq 60 and $i ne length( $seq)-1 ){	
	    print FILE "$temp_line\n";
	    $line_pos=0;
	    $temp_line ='';
	}
	if ( $i eq length( $seq)-1 ){	
	    print FILE "$temp_line\n";
	    $line_pos=0;
	}
	else{
	    $line_pos++;
	}
    }
    close FILE;
}
#=====================================

1;
