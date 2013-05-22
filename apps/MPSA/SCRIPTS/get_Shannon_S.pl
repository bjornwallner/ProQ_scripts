#!/usr/bin/perl -w
use strict;

#use lib '/local/www/services/ProQM/MPSA/SCRIPTS/'; #/afs/pdc.kth.se/home/k/kriil/vol_07/OPM/PERL_SCRIPTS/';
#use lib '/home/simsuede/mpRAP/perllib/';
#use fasta;
#use OPM_module;

my $msa_file         = $ARGV[0];
my $pdb_AA_file      = $ARGV[1];

if ( ! -e $pdb_AA_file ) {
    print STDERR "ignoring beacuse pdb_AA_seq is missing. Probably polyALA-structure\n";
    next;
}

my $MSA_AA_seq      = get_first_seq_from_fasta_seq_file(        $msa_file );      
my $pdb_AA_seq      = get_first_seq_from_fasta_seq_file(       $pdb_AA_file );    
if ( length $MSA_AA_seq ne length $pdb_AA_seq ) { # 1pp9T.fa    # PDB_AA_SEQ AND ORIGINAL_AA_SEQ MAY 	
    print STDERR "sequences are not identical\n";               # BE NON-IDENTICAL AT COMMON SITES
    print STDERR "1) $msa_file\n";                                     
    print STDERR"3) $MSA_AA_seq\n";
    print STDERR"4) $pdb_AA_seq\n";
    exit;
}

# GET AA COUNTS FOR EACH COLUMN IN MSA
my @list_of_counters = get_list_of_aa_counts( $msa_file );
if ($list_of_counters[0] eq "ERROR" ){
    print STDERR "$msa_file contains seqs of different lengths skipping this sequence\n"; 
    exit;  #1jdmA.fa  are probably fixed now when creating psi-blast-msa
}    

#===================================================================
# GET CONSERVATION SCORE AND INDELFREQUENCIES FOR EACH COLUMN IN MSA
my @pdb_Shannon_indel;
my @pdb_Shannon_noindel;
my @pdb_indel;
my $gapped_AA_alphabet   = "ACDEFGHIKLMNPQRSTVXYW-";
my $ungapped_AA_alphabet = "ACDEFGHIKLMNPQRSTVXYW";            
for ( my $i=0;$i< length $MSA_AA_seq ; $i++){ 
    my $chr = substr( $MSA_AA_seq, $i, 1);
    if ( $chr eq "-") {
	print STDERR "WARNING: FIRST SEQUENCE IN PDB_MSA CONTAINS GAP\n";
	print STDERR "$MSA_AA_seq\n";
	exit;
    }
    if ( $chr eq "X" ){
	push( @pdb_Shannon_indel   , "X" );
	push( @pdb_Shannon_noindel , "X" );
	push( @pdb_indel           , "X" );
	next;
    }	
    my %AA_counter  = %{ $list_of_counters[ $i ] };
    my $col_with_gaps   = get_str_from_chr_hash( \%AA_counter,   \$gapped_AA_alphabet );
    my $col_without_gaps= get_str_from_chr_hash( \%AA_counter, \$ungapped_AA_alphabet );	
    my $msa_size               = length $col_with_gaps;
    my $number_of_ungap_in_col = length $col_without_gaps;
    my $gaps_in_col = $msa_size - $number_of_ungap_in_col;
    
    # in 4.0 and 2.0 msa there may be cols with only gaps. They are annotated as "X"; !!!!!!!!
    # Shannon entropy
    my $MSA_Shannon_indel    = get_Shannon_Entropy2( \%AA_counter , \$gapped_AA_alphabet  , \$msa_size );
    my $MSA_Shannon_noindel  = get_Shannon_Entropy2( \%AA_counter , \$ungapped_AA_alphabet, \$number_of_ungap_in_col );
    my $MSA_indel_freq       = $gaps_in_col / $msa_size;	
    push( @pdb_Shannon_indel   , $MSA_Shannon_indel  );
    push( @pdb_Shannon_noindel , $MSA_Shannon_noindel);
    push( @pdb_indel           , $MSA_indel_freq     );
}

#============================================================================    
# CALCULATING MEAN VALUES FOR ALL PDB-SITES  (WHAT ABOUT "X"s?)
my $number_of_aa_s      = 0;
my $indel_freq_sum      = 0;
my $Shannon_indel_sum   = 0;    
my $Shannon_noindel_sum = 0;    
for ( my $i=0;$i< length $MSA_AA_seq; $i++){
    my $chr = substr( $MSA_AA_seq, $i,1);
    next if ($chr eq "X" ); #may be if only gaps in msa column
    
    if ( $pdb_indel[$i] ne "X" ){
	$indel_freq_sum      +=           $pdb_indel[$i];
    }
    if ( $pdb_Shannon_indel[$i] ne "X" ){
	$Shannon_indel_sum   +=   $pdb_Shannon_indel[$i];
 	}
    if ( $pdb_Shannon_noindel[$i] ne "X" ){
	$Shannon_noindel_sum += $pdb_Shannon_noindel[$i];
    }
    $number_of_aa_s++;
}
if ( $number_of_aa_s == 0 ){
    print STDERR "zero residues in seq \n";
    exit;
}
my $mean_indel_freq      =      $indel_freq_sum / $number_of_aa_s;  
my $mean_Shannon_indel   =   $Shannon_indel_sum / $number_of_aa_s;
my $mean_Shannon_noindel = $Shannon_noindel_sum / $number_of_aa_s;
#====================================================================================
# CALCULATE STANDARD DEVIATION
my $indel_std = 0;
foreach my $item ( @pdb_indel ) {
    next if ( $item eq "X" );
    $indel_std += ( $item - $mean_indel_freq ) * ( $item - $mean_indel_freq ); 
}
$indel_std = $indel_std / ($number_of_aa_s - 1 );
$indel_std = sqrt( $indel_std);

# calculating standard deviation of sample
my $Shannon_indel_std = 0;
foreach my $item ( @pdb_Shannon_indel ) {
    next if ( $item eq "X" );
    $Shannon_indel_std += ( $item - $mean_Shannon_indel ) * ( $item - $mean_Shannon_indel ); 
}
$Shannon_indel_std = $Shannon_indel_std / ($number_of_aa_s - 1 );
$Shannon_indel_std = sqrt( $Shannon_indel_std);
my $Shannon_noindel_std = 0;
foreach my $item ( @pdb_Shannon_noindel ) {
    next if ( $item eq "X" );
    $Shannon_noindel_std += ( $item - $mean_Shannon_noindel ) * ( $item - $mean_Shannon_noindel ); 
}
$Shannon_noindel_std = $Shannon_noindel_std / ($number_of_aa_s - 1 );
$Shannon_noindel_std = sqrt( $Shannon_noindel_std);
#======================================================================
# CALCULATE NORMALIZED VALUES FOR INDEL FREQUENCIES AND SHANNON ENTROPY
my @pdb_norm_indel;
my @pdb_norm_Shannon_indel;
my @pdb_norm_Shannon_noindel; 

if ( scalar @pdb_indel ne length $MSA_AA_seq ) {
    print STDERR "length error4\n";
    exit;
}

foreach my $item ( @pdb_indel ) {
    if ( $item eq "X" ){
	push( @pdb_norm_indel, $item );
    }
    else{
	my $new_item = ($item - $mean_indel_freq ) / ($indel_std + 0.0001) ;	    
	push( @pdb_norm_indel, $new_item );
    }
}    
foreach my $item ( @pdb_Shannon_indel ) {
    if ( $item eq "X" ){
	push( @pdb_norm_Shannon_indel, $item );
    }
    else{
	my $new_item = ($item - $mean_Shannon_indel ) / ($Shannon_indel_std + 0.0001);	    
	push( @pdb_norm_Shannon_indel, $new_item );
    }
}
foreach my $item ( @pdb_Shannon_noindel ) {
    if ( $item eq "X" ){
 	    push( @pdb_norm_Shannon_noindel, $item );
    }
    else{
	my $new_item = ($item - $mean_Shannon_noindel ) / ($Shannon_noindel_std + 0.0001);	    
	push( @pdb_norm_Shannon_noindel, $new_item );
    }
}
#=========================================================================
# ALIGN VECTOR OF INDEL AND SHANNON VALUES

my $pdb_norm_Shannon_noindel_str = join( "_", @pdb_norm_Shannon_noindel);
my $pdb_norm_Shannon_indel_str   = join( "_", @pdb_norm_Shannon_indel  );        
my $pdb_abs_Shannon_noindel_str = join( "_", @pdb_Shannon_noindel);
#=========================================================================
# PRINT OUTFILE     
print  ">NORM_Shannon_noindel\n";
print  "$pdb_norm_Shannon_noindel_str\n";
print  ">ABS_Shannon_noindel\n";
print  "$pdb_abs_Shannon_noindel_str\n";
print  ">AA\n";
print  "$MSA_AA_seq\n";



#=========================================================================    
sub get_list_of_aa_counts{
    my ( $msa_infile ) = @_;
    # INIT ONE AA-COUNTER
    my $aa_alphabet = "ACDEFGHIKLMNPQRSTVXYW-";    # obs X- counted
    my %chr_counter;
    for (my $i=0;$i<length $aa_alphabet; $i++){
	my $chr = substr( $aa_alphabet, $i, 1);
	$chr_counter{"$chr"} = 0;
    }
    # INIT ONE AA-COUNTER PER SITE (COLUMN IN MSA)
    my @list_of_counters;
    my $first_AA_seq = get_first_seq_from_fasta_seq_file( $msa_infile );    
    
    #print "A:$first_AA_seq\n";

    for (my $i=0;$i<length $first_AA_seq; $i++){
	push( @list_of_counters, { %chr_counter } );
    }
    #print "@list_of_counters\n";
    # SKIP/EXIT LARGE MSA-FILES
    my $filesize = -s "$msa_infile";
    #print "$filesize\n";
 # NOCHECK FOR THIS 16GB MEMORY SHOULD MANAGE.... /BW    if ( $filesize > 1000000) {
 # NOCHECK FOR THIS 16GB MEMORY SHOULD MANAGE.... /BW 	print STDERR "large_filesize $filesize $msa_infile\n";	
 # NOCHECK FOR THIS 16GB MEMORY SHOULD MANAGE.... /BW 	exit;
 # NOCHECK FOR THIS 16GB MEMORY SHOULD MANAGE.... /BW 	#next;
 # NOCHECK FOR THIS 16GB MEMORY SHOULD MANAGE.... /BW    }    
    #==========================================    
    open( MSA_INFILE , "$msa_infile") || die "cannot open infile \n";
    # SKIP FIRST SEQUENCE
    <MSA_INFILE>;
    <MSA_INFILE>;
    # LOOP ALL OTHER SEQUENCES IN MSA-FILE    
    my $curr_descr= "";
    my $curr_seq  = $first_AA_seq;
    while (<MSA_INFILE>){
	my $line = $_;
	#print "L";
	if ($line =~ /^>([\w\-\/]+)/ ) {
	    $curr_descr = $1;
	    if ( length $curr_seq ne length $first_AA_seq ) {
		print STDERR "WARNING1: length error in $msa_infile!!!!! \n";
		print "$curr_seq\n";
		print "$first_AA_seq\n",
		return ("ERROR","");
	    }
	    # ADD RESIDUE-COUNT FOR EACH SITE (COLUMN) OF CURRENT SEQUENCE 
	    #print "$curr_seq\n";
	    #print "@list_of_counters\n";
	    @list_of_counters = count_chr_in_seq( \@list_of_counters, \$curr_seq );
	    $curr_seq   = "";
	}
	else{	   
	    if ($line =~ /([\w\-\*\.]+)/ ) {	
		$curr_seq .= $1;
		#print "$curr_seq\n";
	    }
	    else {
		print "error\n";
		print "$line";
		exit;
	    }
	}
    }
    @list_of_counters = count_chr_in_seq( \@list_of_counters, \$curr_seq );
    close MSA_INFILE;
    return @list_of_counters;    
}
#=========================================================================    
sub get_Shannon_Entropy2 {
    my ($c, $a, $l ) = @_;
    my %type_counts = %{$c};
    my $alpha       = ${$a};
    my $Shannon     = 0;
    my $col_len     = ${$l};
    for (my $i=0 ; $i<length  $alpha; $i++) {
	my $key = substr($alpha,$i,1);
	my $counts = $type_counts{$key};
	if ( $col_len eq 0 ){
	    return "X";
	}
	my $freq = $counts / $col_len;
	next if $freq eq 0;
	$Shannon = $Shannon + $freq * ( log ( $freq )/ log 2 ) ;
    }
    return -$Shannon;
}
#=========================================================================    
sub get_str_from_chr_hash{
    my ( $h, $a ) =@_;
    my %aa_hash = %{$h};
    my $aa_alphabet = ${$a};
    my $temp_str = "";
    for ( my $i=0; $i< length $aa_alphabet; $i++){
	my $chr = substr( $aa_alphabet, $i, 1);
	for ( my $j=0; $j< $aa_hash{$chr}; $j++){
	    $temp_str .= $chr;
	}
    }
    return $temp_str;
}    
#=========================================================================    
sub count_chr_in_seq{
    my ( $l_c, $s )          = @_;
    my @old_list_of_counters = @{$l_c};
    my $seq                  = ${$s};
    my @new_list_of_counters;
    for (my $i=0;$i<length $seq; $i++){
#	print "*$seq*\n";
	my %temp_counter = %{ $old_list_of_counters[$i] };
	my $chr = substr( $seq, $i, 1);
	$temp_counter{$chr}++;
	push( @new_list_of_counters, { %temp_counter } );
    }
    return @new_list_of_counters;
}
#=========================================================================    



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
