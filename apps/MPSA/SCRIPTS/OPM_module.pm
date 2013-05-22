use strict;
use warnings;

sub drag_seq_1_to_3_using_12_and_23_alignments{
    my ( $ali1, $ali2, $ali3) = @_;
    my $ungapped1 = remove_gaps( $ali1 );
    my $ungapped2 = remove_gaps( $ali2 );
    my $ungapped3 = remove_gaps( $ali3 );
    if ( $ali1 =~ /\-/  or $ali2 =~ /\-/ ){      #no gaps in ali1 and ali2???  but x's is allowed #no!! palign change X to A
	print STDERR "seqs contains gaps\n";
	print STDERR "$ali1\n";
	print STDERR "$ali2\n";
	#exit;
    }    
    #=================================================================
    #    
    # WANT TO ALIGN ALI1 WITH THE USE OF ALIGNMENT BETWEEN ALI2 AND ALI3    Ex     #1j4nB.fa
    # 
    # ali1: CCSTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCSSCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHTTSSCHHHHHHHHHHHHHHHHHHHHHHHHHTTTCTTCCTTCCCCCTTCCTTGGGHHHHHHHHHHHHHHHHHTCSSCCCCCSCHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHTTCCTTTTHHHHHHHHHHHHHHHHHHTTTSCCSSCHHHHHGGGTC
    # ali2: MASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS
    # ali3: ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK
    #
    # ALIGN ali2 and ali3. Check no gaps in new_ali3
    #
    # new_ali2: MASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS----------------------
    # new_ali3: -ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK
    #
    # gap found in new_ali3. This site will be removed in ali1, ali2, new_ali1, new_ali2;
    #
    # new_ali2: ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK
    # new_ali3:ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS----------------------
    # ali1: CSTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCSSCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHTTSSCHHHHHHHHHHHHHHHHHHHHHHHHHTTTCTTCCTTCCCCCTTCCTTGGGHHHHHHHHHHHHHHHHHTCSSCCCCCSCHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHTTCCTTTTHHHHHHHHHHHHHHHHHHTTTSCCSSCHHHHHGGGTC
    # ali2: ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS
    
    
#    print "palign...\n";
    my ($new_ali3,$sub_CIS_sep_seq,$new_ali2) = make_palign_global_sequence_alignment( $ungapped3 ,$ungapped2 ,"S");   
#    print "ALI3:$new_ali3\n";
#    print "ALI2:$new_ali2\n\n";
    
    my $gaps = count( '-', $new_ali3 );
#    print "gaps: $gaps\n";
    my $seqid_ = get_SeqId_of_aligned( $new_ali3, $new_ali2 );
#    print "seqid: $seqid_\n";
    if ( $gaps ne 0 ){  #see example above..
	if ( get_SeqId_of_aligned( $new_ali3, $new_ali2 ) > 0.80 ){  #1ijdA seqid=0.88	    
	    #residue must be removed from the 2 input seqs as well ...	    
#	    print STDERR "WARNING ANNOTATED SEQUENCES CONTAINS INSERTION: THIS WILL BE REMOVED \n\n";
#	    print STDERR "$new_ali3\n";
#	    print STDERR "$new_ali2\n";
#	    print STDERR "changed to ...\n";
	    # remove inserted sites ....
	    my $new3 = "";
	    my $new2 = "";
	    my $new_ungapped1 = "";  #input seqs..
	    my $new_ungapped2 = "";
	    my $ungapped_pos  = 0;
	    for( my $i=0;$i< length $new_ali3; $i++){
		my $chr3 = substr( $new_ali3,$i,1);
		my $chr2 = substr( $new_ali2,$i,1);
		my $ung1 = substr( $ali1, $ungapped_pos, 1);
		my $ung2 = substr( $ali2, $ungapped_pos, 1);
		if ( $chr3 eq "-"){
		    $ungapped_pos++;
		}
		else{
		    $new3 .= $chr3;
		    $new2 .= $chr2;
		    if ( $chr2 ne '-'){
			$new_ungapped1 .= $ung1;
			$new_ungapped2 .= $ung2;
			$ungapped_pos++;
		    }
		}
	    }
	    $new_ali3 = $new3;
	    $new_ali2 = $new2;
	    $ali1 = $new_ungapped1;
	    $ali2 = $new_ungapped2;
#	    print STDERR "$new_ali3\n";
#	    print STDERR "$new_ali2\n";
#	    print STDERR "ungapped\n";
#	    print STDERR "$ali1\n";
#	    print STDERR "$ali2\n\n";
	}
	else{
	    print STDERR "WARNING ANNOTATED SEQUENCES CONTAINS INSERTION \n";
	    print STDERR "$new_ali3\n";
	    print STDERR "$new_ali2\n";
	    exit;
	}
    }
    my $aligned_ali1_3     = align_annot( $new_ali2, $ali2, $ali1  );
    my $final_aligned_ali1 = align_annot( $ali3    , $ungapped3, $aligned_ali1_3 );
    
    #                                         may contain gaps in annotation-seq : should not be so!!!
    if ( remove_gaps( $final_aligned_ali1) ne $ali1 or length $final_aligned_ali1 ne length $ali3){
	print STDERR "error seq changed filename\n";
	print STDERR"$ali2\n";
	print STDERR"$ali1\n";
	print STDERR"$final_aligned_ali1\n";
	print STDERR"$ali2\n";
	print STDERR"$ali3\n";
	print STDERR"$aligned_ali1_3\n";
	exit;
    }
    return $final_aligned_ali1;
}

sub drag_real_value_seq_1_to_3_using_12_and_23_alignments{
    my ( $a1, $a2, $a3) = @_;
    my @ali1 = @{$a1};
    my $ali2 = ${$a2};
    my $ali3 = ${$a3};
    my @ungapped1;

#    print STDERR "hoho\n";

    foreach my $item ( @ali1 ){
	if ( $item eq "X" ){
	    next;
	}
	elsif ( $item eq '-' ){
	    print STDERR "error: gap in aligned Z\n";
	}
	else{
	    push ( @ungapped1,  $item );
	}
    }
    
    my $ungapped2 = remove_gaps( $ali2 );
    my $ungapped3 = remove_gaps( $ali3 );
    if ( $ali2 =~ /\-/ ){      #no gaps in ali1 and ali2???  but x's is allowed #no!! palign change X to A
	print STDERR "seqs contains gaps\n";
	print STDERR "$ali2\n";
	exit;
    }    
    #=================================================================
    #    
    # WANT TO ALIGN ALI1 WITH THE USE OF ALIGNMENT BETWEEN ALI2 AND ALI3    Ex     #1j4nB.fa
    # 
    # ali1: CCSTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCSSCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHTTSSCHHHHHHHHHHHHHHHHHHHHHHHHHTTTCTTCCTTCCCCCTTCCTTGGGHHHHHHHHHHHHHHHHHTCSSCCCCCSCHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHTTCCTTTTHHHHHHHHHHHHHHHHHHTTTSCCSSCHHHHHGGGTC
    # ali2: MASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS
    # ali3: ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK
    #
    # ALIGN ali2 and ali3. Check no gaps in new_ali3
    #
    # new_ali2: MASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS----------------------
    # new_ali3: -ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK
    #
    # gap found in new_ali3. This site will be removed in ali1, ali2, new_ali1, new_ali2;
    #
    # new_ali2: ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK
    # new_ali3:ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS----------------------
    # ali1: CSTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCSSCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHTTSSCHHHHHHHHHHHHHHHHHHHHHHHHHTTTCTTCCTTCCCCCTTCCTTGGGHHHHHHHHHHHHHHHHHTCSSCCCCCSCHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHTTCCTTTTHHHHHHHHHHHHHHHHHHTTTSCCSSCHHHHHGGGTC
    # ali2: ASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTS
   
    
#    print "Hej6\n";
    
    my ($new_ali3,$sub_CIS_sep_seq,$new_ali2) = make_palign_global_sequence_alignment( $ungapped3 ,$ungapped2 ,"S");   

#    print "Hej7\n";

    my $gaps = count( '-', $new_ali3 );
    if ( $gaps ne 0 ){  #see example above..
	if ( get_SeqId_of_aligned( $new_ali3, $new_ali2 ) > 0.80 ){	    
	    #residue must be removed from the 2 input seqs as well ...	    
#	    print STDERR "WARNING ANNOTATED SEQUENCES CONTAINS INSERTION: THIS WILL BE REMOVED \n\n";
#	    print STDERR "$new_ali3\n";
#	    print STDERR "$new_ali2\n";
#	    print STDERR "changed to ...\n";
	    # remove inserted sites ....
	    my $new3 = "";
	    my $new2 = "";
	    my @new_ungapped1     ;    #input seqs..
	    my $new_ungapped2 = "";
	    my $ungapped_pos  = 0;
	    for( my $i=0;$i< length $new_ali3; $i++){
		my $chr3 = substr( $new_ali3,$i,1);
		my $chr2 = substr( $new_ali2,$i,1);		
		my $ung1 =         $ali1[ $ungapped_pos];
		my $ung2 = substr( $ali2, $ungapped_pos, 1);
		if ( $chr3 eq "-"){
		    $ungapped_pos++;
		}
		else{
		    $new3 .= $chr3;
		    $new2 .= $chr2;
		    if ( $chr2 ne '-'){
			push ( @new_ungapped1, $ung1 );
			$new_ungapped2 .= $ung2;
			$ungapped_pos++;
		    }
		}
	    }
	    $new_ali3 = $new3;
	    $new_ali2 = $new2;
	    @ali1 = @new_ungapped1;	    
	    $ali2 = $new_ungapped2;
#	    print STDERR "$new_ali3\n";
#	    print STDERR "$new_ali2\n";
#	    print STDERR "ungapped\n";
#	    print STDERR "@ali1\n";
#	    print STDERR "$ali2\n\n";
	}
	else{
#	    print STDERR "WARNING ANNOTATED SEQUENCES CONTAINS INSERTION \n";
	    print STDERR "$new_ali3\n";
	    print STDERR "$new_ali2\n";
	    print STDERR "SEQID: ", get_SeqId_of_aligned( $new_ali3, $new_ali2 ), "\n";
	    exit;
	}
    }
    
    my $stringed_ali1 = join( "_", @ali1 );
    my $str_aligned_ali1_3     = align_real_value_annot( $new_ali2, $ali2     , $stringed_ali1  );  
    $str_aligned_ali1_3 =~ s/\s+/\_/g;
    my $str_final_aligned_ali1 = align_real_value_annot( $ali3    , $ungapped3, $str_aligned_ali1_3 );
    my @final_aligned_ali1 = split( " ", $str_final_aligned_ali1 );
 
    #                                         may contain gaps in annotation-seq : should not be so!!!
    if ( scalar @final_aligned_ali1 ne length $ali3){
	print STDERR "error seq changed\n";
	print STDERR"$ali2\n";
	print STDERR"@ali1\n";
	print STDERR"$ali2\n";
	print STDERR"$ali3\n";
	print STDERR"$str_aligned_ali1_3\n";
	exit;
    }
    return @final_aligned_ali1;
}

#-------------------------------------------------------------------------------------------
sub make_palign_global_sequence_alignment{
    my ($subali_old1, $subali_old2,$type ) = @_;
    my $subali_old_ungapped1  = remove_gaps( $subali_old1 ) ;
    my $subali_old_ungapped2  = remove_gaps( $subali_old2 ) ;
    
    my $temp1 = $subali_old_ungapped1;
    my $temp2 = $subali_old_ungapped2;    

    my $timer_start = time(); 

    my $prot1_file = "ali1.$timer_start.seq";
    my $prot2_file = "ali2.$timer_start.seq";
    my $temp_descr1 = "prot1";
    my $temp_descr2 = "prot2";
    print_fasta2( $prot1_file,$temp_descr1, $temp1); 
    print_fasta2( $prot2_file,$temp_descr2, $temp2); 
    
#    if (length $subali_old_ungapped1 > 1000) {
#	$temp1 = substr( $subali_old_ungapped1, 0, 1000);
#	$temp1.= "\n";
#	$temp1.= substr( $subali_old_ungapped1, 1000);
#    }
#    if (length $subali_old_ungapped2 > 1000) {
#	$temp2 = substr( $subali_old_ungapped2, 0, 1000);
#	$temp2.= "\n";
#	$temp2.= substr( $subali_old_ungapped2, 1000);
#    }
#    `echo -e ">ali1\n$temp1" > ali1.seq`;
#    `echo -e ">ali2\n$temp2" > ali2.seq`;
#    `module add palign`;
#    my @ali_res = `palign -global -ali -seq ali1.seq - -prof ali2.seq -`;
    
    my @ali_res = `/afs/pdc.kth.se/home/k/kriil/vol_07/palign/linux/palign -global -ali -endsfree  $prot1_file $prot2_file`;
    #always use this alternative??
        
    `rm $prot1_file $prot2_file`;    
    
    #iff palign crashes
    my $len_data = scalar ( @ali_res );
    if ( not( @ali_res ) or  $len_data <= 1 )  {
	print "palign has crashed\n";
	die "Error1 in make_palign_global_sequence_alignment"; 
    }  

    if (  $len_data <= 3 )  {
	if ( $ali_res[0] =~/ERROR> Sequence too long !/ ) {
	    print "palign has crashed\n";
	    print "*@ali_res*\n";
	    die "Error4 in make_palign_global_sequence_alignment"; 
	}
    }  

#    `rm temp.file`;
    
    my $found       =  0 ;
    my $subali_new1 = '' ;
    my $subali_new2 = '' ;
    my $first       = 'X';
    my $last        =  0 ;
    my $aligned     = 'F';
    my $old_pos1    = -1 ;
    my $old_pos2    = -1 ;
    my $OK_ali      = 'T';
    
    foreach my $line ( @ali_res ) {
	#while (   <ALI_RES_FILE>  ) {
	#$line = $_;


#	print "$line";




	if ( $line =~ /^SCORE> / ) {
	    $found = 1;
	    next;
	}
	if ( $found eq 1 ) {
	    my $AA1  = substr( $line, 3, 1 );
	    my $AA2  = substr( $line, 21,1 );
	    my $pos1 = substr( $line, 9, 4 );
	    my $pos2 = substr( $line, 28,4 );
	    $pos1    =~s/\s+//g;
	    $pos2    =~s/\s+//g;

	    if ( $pos1  eq ''){
		$OK_ali = 'F';
	    }
	    
	    if ( $pos1 eq 0 or $pos2 eq 0 ) {
		$aligned = 'F';
	    }
	    else {
		$aligned = 'T';
	    }
	    if ( $aligned eq 'T' and ($old_pos1 eq -1 or $old_pos2 eq -1 )) {
		$first = 'T';
	    }
	    if ( $pos2 ne 0 ) {
		$last = $pos2;
	    }
	    if ( $aligned eq 'T' and $first eq 'T' and $pos2 ne 1 ) {  
		for (my $j = 0; $j < $pos2-1 ; $j++) {
		    $subali_new1.= '-';
		    $subali_new2.= substr( $subali_old_ungapped2, $j , 1 );
		}
	    }			
	    if ( $aligned eq 'T' and $first eq 'F'and $old_pos2+1 ne $pos2 ) {  #if not prev.pos2+1 eq pos2
		for (my $j = $old_pos2; $j < $pos2-1 ; $j++) {
		    $subali_new1.= '-';
		    $subali_new2.= substr( $subali_old_ungapped2, $j , 1 );
		}
	    }			
	    if ( $first eq 'T') {
		$first = 'F';
	    }
	    $subali_new1.= $AA1;
	    $subali_new2.= $AA2;
	    $old_pos1 = $pos1;
	    if ( $pos2 ne 0 and $old_pos2 ne 0) {
		$old_pos2 = $pos2;
	    }
	}
    }
    
    my $sub_sep_seq = '';
    if ( $OK_ali eq 'F') {
	print "ali is not OK. unalign... \n";
	($subali_new1,$sub_sep_seq,$subali_new2) = unalign_sequence( $subali_old_ungapped1,$subali_old_ungapped2,$type);  #??
	print "$subali_new1 \n";
	print "$sub_sep_seq \n";
	print "$subali_new2 \n";
	exit;
    }

    if ( count('-',$subali_old2) eq length($subali_old2) ) { #only gaps (should not be!)
	print "Only gaps in one seq \n";
	exit;
    }
    my $len2 = length ( remove_gaps($subali_old2) );

    if ( $last-1 ne $len2-1 ) {                    #len2
	for (my $j = $last; $j < $len2 ; $j++) {
	    $subali_new1.= '-';
	    $subali_new2.= substr( $subali_old_ungapped2, $j , 1 ); 
	}
    }
        
    for (my $j = 0; $j < length ($subali_new1) ; $j++) {
	$sub_sep_seq .= $type;  
    }
    $subali_new1  =~s/\./\-/g;		#change'.' to '-'
    $subali_new2  =~s/\./\-/g; 
    if ( ($subali_old_ungapped1 ne remove_gaps($subali_new1) ) or ($subali_old_ungapped2 ne remove_gaps($subali_new2)) ){ 

	if ($subali_old_ungapped1 ne remove_gaps($subali_new1) ){
#	    print STDERR "1) \n";
	    my $seq1 = $subali_old_ungapped1 ;
	    my $seq2 = remove_gaps($subali_new1) ;
#	    print STDERR length($seq1)," ",length($seq2), " \n";
	    #for ( my $i=0; $i<length($seq1);$i++){
		#my $chr1 = substr( $seq1, $i, 1);
		#my $chr2 = substr( $seq2, $i, 1);
		#if ($chr1 ne $chr2) {
		#    print STDERR "$i: $chr1 $chr2 \n";
		#}
	    #}
	    if ( length( $subali_old_ungapped1) ne length( remove_gaps($subali_new1)) ){		
#		print STDERR "lengths have changed in seq1. \n";
##		print        "lengths have changed in seq1. \n";
#		print STDERR "ERROR3 in make_palign_global_sequence_alignment \n";		
		die "ERROR MAKE PALIGN GLOBAL SEQUENCE ALIGNMENT - ERROR 3\n";
	    }
	}
	if ($subali_old_ungapped2 ne remove_gaps($subali_new2) ){
#	    print STDERR "2) \n";
	    my $seq1 = $subali_old_ungapped2 ;
	    my $seq2 = remove_gaps($subali_new2) ;
	    my $new_subali_new2 ;
#	    print STDERR length($seq1)," ",length($seq2), " \n";	    
	    my $i = -1;
	    for ( my $j=0; $j<length($subali_new2);$j++){		
		my $chr2 = substr( $subali_new2, $j, 1);
		if ( $chr2 ne '-'){
		    $i++;
		}
		my $chr1 = substr( $seq1, $i, 1);
		if ($chr1 ne $chr2 and $chr1 eq "X" and $chr2 eq "A" ){
#		    print STDERR "Palign has changed X to A at $i .. now change back\n";
		    $chr2 = "X";
		}
		elsif ($chr1 ne $chr2 and $chr2 ne '-') {		
#		    print STDERR "$i: $chr1 $chr2 \n";
		}
		$new_subali_new2 .= $chr2;
	    }	   
	    $subali_new2 = $new_subali_new2;
	    if ( length( $subali_old_ungapped2) ne length( remove_gaps($subali_new2) ) ){		
#		print "lengths have changed in seq2. \n";
		print STDERR "lengths have changed in seq2. \n";
		print STDERR "ERROR4 in make_palign_global_sequence_alignment \n";				
		exit;
	    }
	}	
#	print        "Palign warning: Sequences have been changed, but lengths remain same \n";
	print STDERR "Palign warning: Sequences have been changed, but lengths remain same \n";
#	print STDERR length( $subali_old_ungapped1 )," ",length(remove_gaps( $subali_new1 )),"\n";
#	print STDERR length( $subali_old_ungapped2 )," ",length(remove_gaps( $subali_new2 )),"\n\n";
#	print "ERROR in make_palign_global_sequence_alignment \n";
#	exit;
    }
    return ( $subali_new1, $sub_sep_seq, $subali_new2 );
}
#===============================================
sub unalign_sequence{                    #remove sep_Seq
    my ( $subali1, $subali2, $type ) = @_; 
    my $subali_old_ungapped1  = remove_gaps( $subali1 ) ;
    my $subali_old_ungapped2  = remove_gaps( $subali2 ) ;
    my $len1 = length( $subali_old_ungapped1 );
    my $len2 = length( $subali_old_ungapped2 );
    my $new1 = '';
    my $new2 = '';
    my $new_sep_seq = '';
    for( my $i= 0; $i < $len1; $i++ ) {
	my $chr1 = substr( $subali_old_ungapped1, $i , 1 );
	$new1        .= $chr1;
	$new2        .= '-';
	$new_sep_seq .= $type;
    }
    for( my $i= 0; $i < $len2; $i++ ) {
	my $chr2 = substr( $subali_old_ungapped2, $i , 1 );
	$new1 .= '-';
	$new2 .= $chr2;
	$new_sep_seq .= $type;
    }        
    return ($new1,$new_sep_seq,$new2); 
}
#==============================================================================================================
sub print_fasta2{
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
#==============================================================================================================
sub count{
    my ($chr,$seq) = @_;
    my $counter = 0;
    my $len = length ($seq);
    if ($len eq 0 ) {
	print STDERR "no input in count \n";
	return 0;
    }
    for(my $i= 0; $i<$len; $i++) {
	my $seq_chr = substr($seq,$i,1);
	if ($chr eq $seq_chr) {
	    $counter++;
	}
    }
    return $counter;
}

sub remove_gaps{
    my ($seq) = @_;    
    my $new_seq;
    my $len = length ($seq );
    for (my $i=0; $i< $len; $i++) {
	my $chr = substr($seq,$i,1);
#	next if ($chr eq '-' or $chr eq '.');
	next if ($chr eq '-' );
	$new_seq .= $chr;
    }
    return $new_seq;
}

sub get_SeqId_of_aligned{
    my ($ali_seq1,$ali_seq2) = @_;
    my $ID          = 0;
    my $nonID       = 0;
    my $non_aligned = 0;
    my $chr1        = '';
    my $chr2        = '';
    my $seqID       = 0.0;
    my $ambigs      = 0;
    foreach (my $pos = 0; $pos < length $ali_seq1;++$pos) {
	$chr1 = substr($ali_seq1,$pos,1);
	$chr2 = substr($ali_seq2,$pos,1);
	if ($chr1 eq '.' or $chr1 eq '-' or $chr2 eq '.' or $chr2 eq '-') { 
	    $non_aligned += 1 ;
	}
	elsif ($chr1 eq 'X' or $chr2 eq 'X') { 
	    $ambigs++;
	}
	elsif ($chr1 eq $chr2) {
	    $ID += 1;
	}
	else {
	    $nonID += 1;
	}
    }
    if ($nonID + $ID == 0){
	$seqID = 0;
    }
    else {
	$seqID = $ID / ($nonID + $ID) ;
    }
    return $seqID;
}
#----------------------------------------------------------------------------------

sub align_annot{ 
    my ($ali_seq,$seq,$annot) = @_;
    chomp($seq);
    chomp($annot);
    chomp($ali_seq);
    my $len = length($ali_seq);
    my $len2 = length($seq);
    my $len3 = length($annot);
    
#    if ( not($len2 eq $len3)){ # and $len2 ne 0 and $len ne 0 ) ) {
#	print STDERR "input error in align annot \n";
#	return '';
#    }
    my $ali_pos   = 0;
    my $seq_pos   = 0;
    my $annot_pos = 0;    
    my $ali_chr;
    my $seq_chr;
    my $annot_chr;
    my $new_ali;
    my $new_seq;
    my $new_annot;
    
    #warn if wrong sequence length..!!
    
    while ($ali_pos < $len) {
	$ali_chr   = substr($ali_seq,$ali_pos  ,1);
	$seq_chr   = substr($seq    ,$seq_pos  ,1);
	$annot_chr = substr($annot  ,$annot_pos,1);
	while  ( (uc($ali_chr) ne uc($seq_chr) ) and ($ali_pos < $len) ) {
	    $new_ali   .= $ali_chr;
	    $new_seq   .= '-';
	    $new_annot .= '-';
	    $ali_pos   +=1;    
	    $ali_chr   = substr($ali_seq,$ali_pos  ,1);
	}
	$new_ali   .= $ali_chr;
	$new_seq   .= $seq_chr;
	$new_annot .= $annot_chr;	
	$ali_pos   +=1;
	$seq_pos   +=1;
	$annot_pos +=1;
    }
    return $new_annot;
}

sub align_real_value_annot{ 
    my ($ali_seq,$seq,$annot) = @_;
    chomp($seq);
    chomp($annot);
    chomp($ali_seq);
    my $len          = length($ali_seq);
    my $len2         = length($seq);    
    my @annotations  = split( /\_/, $annot);    
#    my @annotations  = split( /\s+/, $annot);
    my $len3         = scalar( @annotations );
    my $ungapped_ali = remove_gaps( $ali_seq );
    my $len4         = length ($ungapped_ali );
    
    if ( $len2 ne $len3 or $len3 ne $len4 ) {
	print "length error in align real value annot \n";
	print "$annot\n";
	print "$len2 $len3 $len4 \n";
	print "*$annot*\n";
	print "$ali_seq\n";
	print "$seq\n";
	print "$ungapped_ali\n";
	print "@annotations\n";       

	die "LENGTH ERROR IN ANNOTATION SEQUENCE\n";
    }

    my $ali_pos   = 0;
    my $seq_pos   = 0;
#    my $annot_pos = 0;    
    my $ali_chr;
    my $seq_chr;
    my $annot_val;
    my $new_ali;
    my $new_seq;
    my $new_annot;
        
    chomp( $seq ) ;
    chomp( $ali_seq );

    while ($ali_pos < $len) {
	$ali_chr   = substr($ali_seq,$ali_pos  ,1);	
	if ( $seq_pos < $len2 ){  #if outside string
	    $seq_chr   = substr($seq    ,$seq_pos  ,1);
	    $annot_val = $annotations[$seq_pos];
	}
	else {
	    $seq_chr   = '-';
	    $annot_val = '-';
	}
	while  ( (uc($ali_chr) ne uc($seq_chr) ) and ($ali_pos < $len) ) {
	    $new_ali    .= $ali_chr;
	    $new_seq    .= '-';
	    #$new_annot .= '-_';
	    $new_annot  .= '- ';
	    $ali_pos   +=1;    
	    $ali_chr   = substr($ali_seq,$ali_pos  ,1);
	}
	$new_ali   .= $ali_chr;
	$new_seq   .= $seq_chr;
	#$new_annot .= "$annot_val"."_";	
	$new_annot .= "$annot_val"." ";	
	$ali_pos   +=1;
	$seq_pos   +=1;  #may plus outside string position

    }
    chop( $new_annot ); #remove last '_'
#    print "$ali_seq\n";
 #   print "$seq \n";
 #   print "$new_annot\n";
    return $new_annot;
}

sub DSSP2HSL{
    my ($old_seq) = @_;    
    my %transl;
    $transl{"H"} = "H"; # H = alpha helix
    $transl{"G"} = "H"; # G = 3-helix (3/10 helix)
    $transl{"I"} = "H"; # I = 5 helix (pi helix)
    $transl{"E"} = "S"; # E = extended strand, participates in beta ladder
    $transl{"B"} = "L"; # B = residue in isolated beta-bridge
    $transl{"C"} = "L"; # C = random coil. A blank in the DSSP secondary structure determination stands for loop or irregular.
    $transl{"T"} = "L"; # T = hydrogen bonded turn
    $transl{"S"} = "L"; # S = bend    
    $transl{"X"} = "X"; # S = bend    
    my $new_seq;
    foreach (my $pos = 0; $pos < length $old_seq;++$pos) {
	my $chr = substr($old_seq,$pos,1);
	if( defined( $transl{$chr} ) ){
	    $new_seq .= $transl{$chr};
	}
	else {
	    $new_seq .= 'X';
	    print STDERR "DSSP2HSL $chr -> X \n";
	    #exit
	}
    }
    return $new_seq;
}
#-------------------------------------------------------------------------------------------
sub get_seq_tract_info_hash{
    my ( $seq,$start_pos,$stop_pos) = @_;
    my %tract_info;
    $tract_info{"TYPE"  } = substr( $seq , $start_pos, 1 );
    $tract_info{"START" } = $start_pos;
    $tract_info{"STOP"  } = $stop_pos;
    $tract_info{"LENGTH"} = $stop_pos-$start_pos+1;
    return (%tract_info);
}
#-------------------------------------------------------------------------------------------
sub get_list_of_tract_info_hash{   
    my ( $seq ) = @_;    
    my $len  = length $seq;
    my $pos              = 0 ;    
    my $start_pos        = 0;
    my $stop_pos         = -1;
    my @list_of_tract_info;
    my $last_position = $len - 1;
    
    while ($pos < $last_position ) {	   
	my $chr = substr( $seq, $pos, 1);
	my $tract_type = $chr;
	my $temp_seq = '';
	$start_pos = $pos;
	while ($chr eq $tract_type) {                       
	    $temp_seq .= $chr;
	    last if ($pos eq $last_position);
	    $pos++;
	    $chr = substr( $seq, $pos, 1);
	}
	my $previous_char = substr( $seq,$pos-1, 1);
	my $current_char  = substr( $seq,$pos  , 1);
	if     ($pos eq $last_position and $previous_char eq $current_char ) {   #last position do not change
	    $stop_pos = $pos;
	}
	elsif ($pos eq $last_position and $previous_char ne $current_char) { #last position do change
	    $stop_pos = $pos-1;
	    my %tract_info = get_seq_tract_info_hash( $seq,$start_pos,$stop_pos);
	    push( @list_of_tract_info, {%tract_info} );	    
	    $tract_type = $chr;
	    $start_pos  = $pos;
	    $stop_pos   = $pos;
	}
	elsif ( $pos ne $last_position ) { #middle position character
	    $stop_pos = $pos-1;
	}
	my %tract_info = get_seq_tract_info_hash( $seq,$start_pos,$stop_pos);
	push( @list_of_tract_info, {%tract_info} );	   
    }          
    return ( @list_of_tract_info );
}

#------------------------------------------------------------------------------------
sub get_list_of_simple_tract_info{   #only used for separating helices/strand/loop-tracts
   my ( $seq ) = @_;    
   my $len  = length ( $seq );    
   my $pos  = 0;
   my $sum  = 0;
   my $stop = 0;
   my $tract_info_row = '';
   my @list_of_rows;

   if ($len eq 1) {
       return "$seq 0 0 1 1 1\n";
   }
   while ( $pos < $len-1 ) {
	my $chr = substr( $seq, $pos, 1);
	my $tract_type = $chr;
	my $tract_length = 0;
	my $start = $pos;
	my $temp_seq .= '';
	while ($chr eq $tract_type) {
	    $temp_seq .= $chr;
	    $tract_length++;
	    last if ($pos eq $len-1);
	    $pos++;
	    $chr = substr( $seq, $pos, 1);
	}
	$sum += $tract_length;
	my $prev = substr( $seq,$pos-1, 1);
	my $next = substr( $seq,$pos  , 1);
	if ($pos eq $len-1 and $prev eq $next ) {  #last position do not change
	    $stop = $pos;
	}
	elsif ($pos eq $len-1 and $prev ne $next ) {  #last position do  change
#	    print "analayze tract: last res changes \n";
	    $stop = $pos-1;
	    $tract_info_row = "$tract_type $start $stop $tract_length $sum $len \n";
	    push( @list_of_rows, $tract_info_row ); 
	    $tract_type = $chr;
	    $start = $pos;
	    $stop  = $pos;
	    $tract_length = 1;
	    $sum += $tract_length;
	}
	elsif ($pos ne $len-1 and $prev eq $next ) {  #not last position  character did not change
	    print "Error in analyze tracts.... pos ne len-1 and prev eq next \n";
	    print "$seq \n";
	    exit;
	    die "ANALYZE_TRACTS_ERROR";
	}
	elsif ( $pos ne $len-1 and $prev ne $next ) { #not last position  character did change
	    $stop = $pos-1;
	}

	$tract_info_row = "$tract_type $start $stop $tract_length $sum $len \n";
	push( @list_of_rows, $tract_info_row ); 
    }
    if ( $sum ne $len )            {
	print "ANALYZE_TRACTS_ERROR \n";  #put together HSL-sequence from list_of_rows = HSL-seq?
	exit;
    }
    return ( @list_of_rows );
}
#=========================


sub change_gap_to_X_in_list{
    my (@list) = @_;
    my @new_list;
    for my $item ( @list){
   	if ( $item eq '-' ) {
    	    $item = "X";
    	}
    	push( @new_list, $item );
    }
    return @new_list;
}


1;
