=Head NAME

Bio::Pdb - Object for handling pdb-files

=head1 VERSION

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Overview

=head1 AUTHOR

Björn Wallner
bjorn@sbc.su.se

=head1 COPYRIGHT

Copyright (c)  2001, Björn Wallner. All Rights reserved.
This module is free software. It may be used, redistributed
and/or modified under the same terms as Perl itself.

=cut

package Pdb;
use lib '/afs/pdc.kth.se/home/k/kriil/vol_07/OPM/PERL_SCRIPTS/';

#use OPM_module;
#use fasta;

use strict;


my %aa321=('ALA', 'A',
	'ARG', 'R',
	'ASN', 'N',
	'ASP', 'D',
	'CYS', 'C',
	'GLN', 'Q',
	'GLU', 'E',
	'GLY', 'G',
	'HIS', 'H',
	'ILE', 'I',
	'LEU', 'L',
	'LYS', 'K',
	'MET', 'M',
	'PHE', 'F',
	'PRO', 'P',
	'SER', 'S',
	'THR', 'T',
	'TRP', 'W',
	'TYR', 'Y',
	'VAL', 'V',
        'ASX', 'B',
        'GLX', 'Z',
        'XXX', 'A',
        'MSE', 'M',
        'FME', 'M',
        'PCA', 'E',
        '5HP', 'E',
        'SAC', 'S',
        'CCS', 'C');

my %standard_aa321=('ALA', 'A',
	'ARG', 'R',
	'ASN', 'N',
	'ASP', 'D',
	'CYS', 'C',
	'GLN', 'Q',
	'GLU', 'E',
	'GLY', 'G',
	'HIS', 'H',
	'ILE', 'I',
	'LEU', 'L',
	'LYS', 'K',
	'MET', 'M',
	'PHE', 'F',
	'PRO', 'P',
	'SER', 'S',
	'THR', 'T',
	'TRP', 'W',
	'TYR', 'Y',
	'VAL', 'V');



my %aa123=('A','ALA',
	   'R','ARG', 
	   'N','ASN', 
	   'D','ASP',
	   'C','CYS',
	   'Q','GLN', 
	   'E','GLU', 
	   'G','GLY', 
	   'H','HIS',
	   'I','ILE',
	   'L','LEU', 
	   'K','LYS', 
	   'M','MET', 
	   'F','PHE',
	   'P','PRO',
	   'S','SER', 
	   'T','THR', 
	   'W','TRP',
	   'Y','TYR', 
	   'V','VAL',
           'B','ASX',
           'Z','GLX');

 my %atomtypes=("ALA C",1,
 "ALA CA",1,
 "ALA CB",1,
 "ALA N",1,
 "ALA O",1,
 "ARG C",1,
 "ARG CA",1,
 "ARG CB",1,
 "ARG CD",1,
 "ARG CG",1,
 "ARG CZ",1,
 "ARG N",1,
 "ARG NE",1,
 "ARG NH1",1,
 "ARG NH2",1,
 "ARG O",1,
 "ASN C",1,
 "ASN CA",1,
 "ASN CB",1,
 "ASN CG",1,
 "ASN N",1,
 "ASN ND2",1,
 "ASN O",1,
 "ASN OD1",1,
 "ASP C",1,
 "ASP CA",1,
 "ASP CB",1,
 "ASP CG",1,
 "ASP N",1,
 "ASP O",1,
 "ASP OD1",1,
 "ASP OD2",1,
 "CYS C",1,
 "CYS CA",1,
 "CYS CB",1,
 "CYS N",1,
 "CYS O",1,
 "CYS SG",1,
 "GLN C",1,
 "GLN CA",1,
 "GLN CB",1,
 "GLN CD",1,
 "GLN CG",1,
 "GLN N",1,
 "GLN NE2",1,
 "GLN O",1,
 "GLN OE1",1,
 "GLU C",1,
 "GLU CA",1,
 "GLU CB",1,
 "GLU CD",1,
 "GLU CG",1,
 "GLU N",1,
 "GLU O",1,
 "GLU OE1",1,
 "GLU OE2",1,
 "GLY C",1,
 "GLY CA",1,
 "GLY N",1,
 "GLY O",1,
 "HIS C",1,
 "HIS CA",1,
 "HIS CB",1,
 "HIS CD2",1,
 "HIS CE1",1,
 "HIS CG",1,
 "HIS N",1,
 "HIS ND1",1,
 "HIS NE2",1,
 "HIS O",1,
 "ILE C",1,
 "ILE CA",1,
 "ILE CB",1,
 "ILE CD1",1,
 "ILE CG1",1,
 "ILE CG2",1,
 "ILE N",1,
 "ILE O",1,
 "LEU C",1,
 "LEU CA",1,
 "LEU CB",1,
 "LEU CD1",1,
 "LEU CD2",1,
 "LEU CG",1,
 "LEU N",1,
 "LEU O",1,
 "LYS C",1,
 "LYS CA",1,
 "LYS CB",1,
 "LYS CD",1,
 "LYS CE",1,
 "LYS CG",1,
 "LYS N",1,
 "LYS NZ",1,
 "LYS O",1,
 "MET C",1,
 "MET CA",1,
 "MET CB",1,
 "MET CE",1,
 "MET CG",1,
 "MET N",1,
 "MET O",1,
 "MET SD",1,
 "PHE C",1,
 "PHE CA",1,
 "PHE CB",1,
 "PHE CD1",1,
 "PHE CD2",1,
 "PHE CE1",1,
 "PHE CE2",1,
 "PHE CG",1,
 "PHE CZ",1,
 "PHE N",1,
 "PHE O",1,
 "PRO C",1,
 "PRO CA",1,
 "PRO CB",1,
 "PRO CD",1,
 "PRO CG",1,
 "PRO N",1,
 "PRO O",1,
 "SER C",1,
 "SER CA",1,
 "SER CB",1,
 "SER N",1,
 "SER O",1,
 "SER OG",1,
 "THR C",1,
 "THR CA",1,
 "THR CB",1,
 "THR CG2",1,
 "THR N",1,
 "THR O",1,
 "THR OG1",1,
 "TRP C",1,
 "TRP CA",1,
 "TRP CB",1,
 "TRP CD1",1,
 "TRP CD2",1,
 "TRP CE2",1,
 "TRP CE3",1,
 "TRP CG",1,
 "TRP CH2",1,
 "TRP CZ2",1,
 "TRP CZ3",1,
 "TRP N",1,
 "TRP NE1",1,
 "TRP O",1,
 "TYR C",1,
 "TYR CA",1,
 "TYR CB",1,
 "TYR CD1",1,
 "TYR CD2",1,
 "TYR CE1",1,
 "TYR CE2",1,
 "TYR CG",1,
 "TYR CZ",1,
 "TYR N",1,
 "TYR O",1,
 "TYR OH",1,
 "VAL C",1,
 "VAL CA",1,
 "VAL CB",1,
 "VAL CG1",1,
 "VAL CG2",1,
 "VAL N",1,
 "VAL O",1);





=head2 new

   Title    : new
   Usage    : $pdb    = Bio::Pdb->new( -file => 'pdbXXXX.ent',
                                       -chain => 'A',
                                       -pdbcode => 'XXXX');

   Function : Returns a new pdb object for
              the pdb file defined by the file string or from the pdbcode
   Returns  : a new Bio::Pdb object

=cut

sub new 
{
    my $FTPURL="ftp://ftp.ebi.ac.uk/pub/databases/pdb/all_entries/compressed_files/";
    my($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    my $self = {};
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
    my $nohydrogen=0;
    my $skip_distance=0;
    my $skip_altloc=0;
    my $skip_hetatoms=0;
    my $only_standard_aa=0;

    #print "muuh!\n";

    if(defined($param{-nohydrogen}))
    {
	$nohydrogen=1;
    }
    if(defined($param{-skip_dist}))
    {
	$skip_distance=$param{-skip_dist};
    }
    if(defined($param{-skip_altloc}))
    {
	$skip_altloc=$param{-skip_altloc};
    }
    if(defined($param{-skip_hetatoms}))
    {
	$skip_hetatoms=$param{-skip_hetatoms};
    }
    if(defined($param{-only_standard_amino_acids}))
    {
	$only_standard_aa=$param{-only_standard_amino_acids};
    }

    if(defined($param{-file}))
    {
	if(defined $param{-chain})
	{
	    #$self=_read_pdb($param{-file},uc($param{-chain}),$nohydrogen,$skip_distance,$skip_altloc,$skip_hetatoms,$only_standard_aa);
	    $self=_read_pdb($param{-file},($param{-chain}),$nohydrogen,$skip_distance,$skip_altloc,$skip_hetatoms,$only_standard_aa);
	}
	else
	{
	    $self=_read_pdb($param{-file}," ",$nohydrogen,$skip_distance,$skip_altloc,$skip_hetatoms,$only_standard_aa);
	}
    }
    elsif(defined($param{-pdbcode}))
    {
	my $file="pdb$param{-pdbcode}.ent";
	`ncftpget $FTPURL$file`;
	if(defined $param{-chain})
	{
	    #$self=_read_pdb($file,uc($param{-chain}),$nohydrogen,$skip_distance,$skip_altloc,$skip_hetatoms,$only_standard_aa);
	    $self=_read_pdb($file,($param{-chain}),$nohydrogen,$skip_distance,$skip_altloc,$skip_hetatoms,$only_standard_aa);
	}
	else
	{
	    $self=_read_pdb($file," ",$nohydrogen,$skip_distance,$skip_altloc,$skip_hetatoms,$only_standard_aa);
	}
	`rm $file`;
	
    }
    if($self)
    {
	bless ($self,$class);
	return $self;
    }
    else
    {
	return 0;
    }
}
=head2 _read_pdb

  Title   : _read_pdb
  Usage   : Internal function used by the new function.
  Returns : A a pdb object to be blessed.

=cut

sub _read_pdb 
{
    my ($file,$chain,$nohydrogen,$skip_distance,$skip_alt_loc,$skip_hetatoms,$only_standard_aa)=@_;
    #print "\ninput: $file $chain\n";
    my $no_seqres=1;   # Deals with case when the pdb SEQRES is missing
    my @list=split(/\//,$file);
    my $name=$list[length(@list)-1];
    my @pdbfield=();
    my @pdbatomno=();
    my @pdbatomtype=();
    my @pdbatom_alt_loc=();
    my @pdbres=();
    my @pdbresno=();
    my @pdbinsertion_code=();
    my @x=();
    my @y=();
    my @z=();
    my @resbegin_temp=();
    my @resCA_temp=();
    my @resbegin=();
    my @resCA=();
    my @bfactor=();
    my $seqres="";
    my $seqca="";
    my $SeqRes;
    #my $SeqCA;
    #my $SeqCA2;
    my $SeqRes_ali;
    my $SeqCA_ali;
    my $atomcount=0;
    my $resname="";
    my $old_resname="undef";
    my $old_alt_loc=" ";
    my $old_res="undef";
    my $last_added="";
    my @PDBHEADER=('^HEADER','^COMPND','^SOURCE','^AUTHOR','^TITLE',
		   ,'^REVDAT','REMARK', '^JRNL','^HET ', '^HETNAM','^MODRES',
		   ,'^HELIX','^SHEET','^LOOP','^TURN', 'MODEL');    # Specifies which parts of the pdb header to take.
    my @pdbhead=();
    my $seq_begin="";   # A sequence for the residue of first atom in a residue,
                        # for doing an alignment with CA_res to determine which residues to
                        # print out.
    #my $success=1;
    #open(FILE,"$file") or $success=0;
    #my $done=0;
    #if($success)
    #{
    my $chain_found=0;
    my $first_chain="undef";
#    if($chain eq " ") # if the specified chain is " ", it checks that, that chain exists and if it doesn't it takes the first chain.
#    {
	#print $file."\n";
	#my $out=`cat $file`;
	#print $out;
    #print "CHAIN: $chain\n";
    open (FILE,"$file")|| die "Cannot open $file.\n";
    while(<FILE>)     # checking if the specified chain exist otherwise try the first chain
    {
	if(/^ATOM/)
	{
	    my $current_chain=substr($_,21,1);
	    my $res=substr($_,17, 3);
	#    if($file=~/pdbblast/)
	#    {my $dir=`pwd`;
	#	print $dir;
	#    }
#	    print "$file $current_chain $res\n";
	    if($first_chain eq "undef" && defined($aa321{$res}))
	    {
		$first_chain=$current_chain;
		#print "FIRST CHAIN is \"$first_chain\"\n";
	    }

	   # print "hej\n";
	    if($chain eq "*" || $current_chain eq $chain && defined($aa321{$res}))
	    {
	#	print "hej\n";
		$chain_found=1;
		last;
	    }
	}
	
    }
#    if($first_chain eq "undef")
#    {
#	warn "The chain \"$chain\" does not exist in the pdb file: $file\nNo object created...Bla..\n";
#	return 0;
#    }
    if(not($chain_found))
    {
	warn "Specified chain: \"$chain\" not found \n";
	return 0;
	#warn "Specified chain: \"$chain\" not found taking the first chain $first_chain instead.\n";
	#$chain=$first_chain;
    }
	
#    }
    close(FILE);
    open(FILE,"$file") || die "Cannot open $file.\n";
    my $reached_TER=0;
    while(<FILE>)
    {
	chomp;
	#print "$_\n";
	foreach my $key (@PDBHEADER)
	{
	    if(/$key/)
	    {
		push(@pdbhead,$_);
		last;
	    }
	}
	
	if(/^SEQRES/)
	{
	    my $current_chain=substr($_,11,1);
	    if($current_chain eq $chain || $chain eq "*")
	    {
		$no_seqres=0;
		my $atoms=substr($_, 19, 52);
		my @list=split(/\s+/,$atoms);
		for(my $i=0;$i<=$#list;$i++)
		{
		    if($only_standard_aa)
		    {
			$seqres.=$aa321{$list[$i]} if(defined($standard_aa321{$list[$i]}));
		    }
		    else
		    {
			$seqres.=$aa321{$list[$i]} if(defined($aa321{$list[$i]}));
		    }
		}	
	    }
	}

	if(not($reached_TER) && ((/^ATOM/)))# || (/^HETATM/ && (/MSE/ || /PCA/ || /5HP/ || /SAC/ || /CCS/ || /FME/))))
	{
	    #print $_."\n";
	    my $current_chain=substr($_,21,1);
	    if($current_chain eq $chain || $chain eq "*" || $current_chain eq "-")
	    {
	#	print;
	#	print "\n";
		my $field = substr($_,0,6);
		$field=~s/ //g;
		my $atomno=substr($_, 6, 5);
		my $atomtype=substr($_, 12, 4);
		$atomtype=~s/ //g;
		my $alt_loc=substr($_,16,1);
		my $res=substr($_,17, 3);
		#next if(not(defined($standard_aa321{$res})));
		my $resno=substr($_, 22, 4);
		#A hack that fixes a problems with negtive resno from structal (minus sign is taken as the chain identifier)
		if($skip_alt_loc==1 && $current_chain eq "-") 
		{
		    $resno=~s/ //g;
		    $resno=sprintf("%4s",-$resno);
		}
		my $insertion_code=substr($_,26,1);
		$resname="$resno$insertion_code";
		$atomno=~s/ //g;
		$res=~s/ //g;
	#	print "1: $_";
	#	print "\n";
		$old_alt_loc=" " if($resname ne $old_resname);
		next if(not(defined($standard_aa321{$res})));
		next if($skip_hetatoms &&
			($res=~/UNK/ || $res=~/GLX/ || $res=~/ASX/ || $res=~/MSE/ || $res=~/PCA/ || $res=~/5HP/ || $res=~/SAC/ || $res=~/CCS/ || $res=~/FME/));
		next if($nohydrogen==1 &&
			$atomtype=~/HA/ ||
			$atomtype=~/HB/ ||
			$atomtype=~/HD/ ||
			$atomtype=~/HE/ ||
			$atomtype=~/HG/ ||
			$atomtype=~/HH/ ||
			$atomtype=~/HZ/ ||
			$atomtype=~/^H/);
		
		next if($atomtype=~/E21/ ||
			$atomtype=~/E22/ ||
			$atomtype=~/D21/ ||
			$atomtype=~/D22/);
		next if($skip_alt_loc==1 &&
			$resname eq $old_resname &&
			$alt_loc ne $old_alt_loc &&
			$alt_loc ne " " &&
			$old_alt_loc ne " ");
		
		#print "2: $_";
		#print "\n";
		if(defined($aa321{$res}))
		{
		    $resno=~s/ //g;
		    $pdbfield[$atomcount]=$field;
		    $pdbatomno[$atomcount]=$atomno;
		    $pdbatomtype[$atomcount]=$atomtype;
		    if($skip_alt_loc)
		    {
			$pdbatom_alt_loc[$atomcount]=" ";
		    }
		    else
		    {
			$pdbatom_alt_loc[$atomcount]=$alt_loc ;
		    }

		    $pdbres[$atomcount]=$res;
		    $pdbresno[$atomcount]=$resno;
		    $pdbinsertion_code[$atomcount]=$insertion_code;
		    $x[$atomcount]=substr($_,30,8);
		    $y[$atomcount]=substr($_,38,8);
		    #print "$res $chain $resno $y[$atomcount]\n";
		    $z[$atomcount]=substr($_,46,8);
		    $x[$atomcount]=~s/ //g;
		    $y[$atomcount]=~s/ //g;
		    $z[$atomcount]=~s/ //g;
		    if(length($_)>=66)
		    {
			$bfactor[$atomcount]=substr($_,60,6);
		    }
		    else
		    {
			$bfactor[$atomcount]=9.99;
		    }
		    #print $_."\n" if($atomtype eq "CA ");
		    if($x[$atomcount]=~/\d+\-\d+/ ||
		       $y[$atomcount]=~/\d+\-\d+/ ||
		       $z[$atomcount]=~/\d+\-\d+/)
		    {
			warn "Skipping\n";
			return 0;

		    }
		    
		    #THIS TRANSLATES ALL HETATM TO ATOM AND TO THE CORRESPONDING RESIDUE
		    if($field=~/HETATM/)
		    {
			$pdbres[$atomcount]=$aa123{$aa321{$res}};
			$pdbfield[$atomcount]="ATOM";
		    }
		    if($resname ne $old_resname || $atomcount == 0)
			#if($atomtype eq "N  " || $atomcount == 0)
		    {
			#print "$atomtype\n";
			push(@resbegin_temp,$atomcount);
			$seq_begin.=$aa321{$res}; #store the "FIRST ATOM IN RES" sequence for which coordinates exists.
			#print "$res $resname $resno $atomcount $atomtype\n";
		    }
		    if($resname eq $old_resname && $res ne $old_res)
		    {

			warn "[WARNING] Residue changed from $old_res to $res resname ($resname) is still unchanged, possible error in the pdb!\n";
		    }
		    if($atomtype eq "CA" && $last_added ne $resname)
			#if($atomtype eq "CA ")
		    {
			if(scalar @resCA_temp != 0)
			{
			    #print $_."\n";
			    my $last_index=$#resCA_temp;
			    #my $dist=distance($x[$atomcount],$y[$atomcount],$z[$atomcount],$x[$resCA_temp[$#resCA_temp]],$y[$resCA_temp[$#resCA_temp]],$z[$resCA_temp[$#resCA_temp]]);
			    my $dist=3.4;
			    if($skip_distance != 0 && $dist>$skip_distance)
			    {
				warn "$file Distance between CA $last_added $resname is $dist A.\n"
			    }
			    #if($dist>8 || $dist < 2)
			    #{
				#warn "Distance between CA $last_added $resname is $dist A.\n";

				#print "$file Distance between CA $last_added $resname is $dist A.\n" if($resname-$last_added==1);
				#print "($x[$atomcount],$y[$atomcount],$z[$atomcount],$x[$resCA_temp[$#resCA_temp]],$y[$resCA_temp[$#resCA_temp]],$z[$resCA_temp[$#resCA_temp]])\n";
				#print $_."\n";
			    #}
			}
			#print "$x[$atomcount],$y[$atomcount],$z[$atomcount],$x[$resCA_temp[$#resCA]],$y[$resCA_temp[$#resCA]],$z[$resCA_temp[$#resCA]]"."\n";
			$seqca.=$aa321{$res};  #store the CA sequence for which coordinates exists.
			push(@resCA_temp,$atomcount);
			#print " $resno $atomcount $atomtype".scalar @resCA_temp." ".scalar @resbegin_temp."\n";
			$last_added=$resname;
		    }
		    $atomcount++;
		    $old_res=$res;
		    $old_resname=$resname;
		    $old_alt_loc="$alt_loc";
		}
	    }
	}
	$reached_TER=1 if((/^TER/ || /^ENDMDL/) && $atomcount>0);   # In order to avoid reading HETATM after termination of residue coordinates.
    }
    
    if($atomcount == 0)
    {
	warn "The chain \"$chain\" does not exist in the pdb file: $file\nNo object created.\n";
	return 0;
    }
    else
    {
   			
	# Aligment between CA res and "first atom res", to filter out the
	# cases when CA is missing, when CA is missing but some other atom is present that atom
	# is put into the CA position.
	my $seqca2=$seqca;	

	#print "$seq_begin\n$seqca\n\n";
#	exit;

	my @data = `/afs/pdc.kth.se/home/k/kriil/vol_07/OPM/PERL_SCRIPTS/palignpair.pl $seq_begin $seqca`;	
	my $Seq_begin_ali2= $data[0];
	my $SeqCA_ali2    = $data[1];
	chomp( $Seq_begin_ali2 );
	chomp(  $SeqCA_ali2    );
	
	#my ($Seq_begin_ali2,$sep_seq,$SeqCA_ali2)=make_palign_global_sequence_alignment($seq_begin,$seqca ,"S");   

	###my($Seq_begin_ali2,$SeqCA_ali2)=_align($seq_begin,$seqca);
#	print "$Seq_begin_ali2\n\n$SeqCA_ali2\n\n";
	my $len1 = scalar @resbegin_temp;
	my $len2 = scalar @resCA_temp;
	if($len1 != $len2)
	{
	    print STDERR "CA is probably missing in some residue $file $len1 $len2\n";
	    print STDERR "Don't worry the program should be able to handle this\n";
	    print STDERR "Displaying alignment between sequence and CA-sequence\n";
	    print STDERR "$Seq_begin_ali2\n\n$SeqCA_ali2\n\n";
	    
	    #exit(1);
	}
	    
	my @ali=split(//,$SeqCA_ali2);  # Aligment between CA res and "first atom res", to filter out the
	                                # cases when CA is missing.
	my @ali_begin=split(//,$Seq_begin_ali2);
	my $CA_index=0;
	my $resbegin_index=0;
	for(my $i=0;$i<scalar @ali;$i++)
	{
	    if($ali[$i] eq "-" && $ali_begin[$i] ne "-")  #if CA not aligned with "first atom" in resdiue change CA pointer to $resbegin
	    {
		    #changed 2002-01-29 /BW
		$resbegin[$i]=$resbegin_temp[$resbegin_index];
		$resCA[$i]=$resbegin_temp[$resbegin_index];  #"-";
		$ali[$i]=$ali_begin[$i];  #OBSERVE this adds another residue to seqca which in fact is not CA but the first atom belonging to that res.
		$resbegin_index++;
		
		#$resbegin[$i]="-";
		#$resCA[$i]="-";
		#print "IF: $ali[$i] $resbegin[$i] $resCA[$i]\n";
	    }
	    else
	    {
		$resbegin[$i]=$resbegin_temp[$resbegin_index];
		$resCA[$i]=$resCA_temp[$CA_index];
		#print "ELSE: $ali[$i] $resbegin[$i] $resCA[$i]\n";
		#my $temp=$aa321{$pdbres[$resbegin_temp[$j]]};
		$CA_index++;
		$resbegin_index++;
	    }
	}
	
	$seqca=join('',@ali);
	#print $seqca."\n";	

#	$SeqCA = Bio::Seq->new( -moltype => "protein",
#				-seq => $seqca,
#				-id  => 'SEQCA',);	
#	$SeqCA2 = Bio::Seq->new( -moltype => "protein",
#				 -seq => $seqca2,
#				 -id  => 'SEQCA2',);
	if(length($seqres)>0)   #Produce Alignment
	{
#	    $SeqRes = Bio::Seq->new( -moltype => "protein",
#				     -seq => $seqres,
#				     -id  => 'SEQRES',);

	    my @data = `/afs/pdc.kth.se/home/k/kriil/vol_07/OPM/PERL_SCRIPTS/palignpair.pl $seqres $seqca`;	
	    my $SeqRes_ali= $data[0];
	    my $SeqCA_ali = $data[1];
	    chomp( $SeqRes_ali );
	    chomp( $SeqCA_ali  );

	    #my ($SeqRes_ali,$sep_seq,$SeqCA_ali)=make_palign_global_sequence_alignment($seqres,$seqca,"S");   
	    #($SeqRes_ali,$SeqCA_ali)=_align($seqres,$seqca);
#	    print "$SeqRes_ali\n$SeqCA_ali\n";
	}
	else
	{
	    $no_seqres=1;
	    #print "Hello\n";
#	    $SeqRes=0;
	    $SeqRes_ali=0;
	    $SeqCA_ali=$seqca;
	}

	
    }
    my $self={name => $name,
	      header => [@pdbhead],
	      chain => $chain,
#	      seqres => $SeqRes,
	      no_seqres => $no_seqres,

	      seqCA  => $seqca,

#	      seqCA => $SeqCA,
#	      seqCA2 => $SeqCA2,
	      seqres_ali => $SeqRes_ali,
	      seqCA_ali => $SeqCA_ali,
	      field => [@pdbfield],
	      atomno => [@pdbatomno],
	      atomtype => [@pdbatomtype],
	      atom_alt_loc=> [@pdbatom_alt_loc],
	      res => [@pdbres],
	      resno => [@pdbresno],
	      insertion_code => [@pdbinsertion_code],
	      x => [@x],
	      y => [@y],
	      z => [@z],
	      resbegin => [@resbegin],   #NO:pointer to begin of a residue if == "-" CA is missing for that residues
	      resCA => [@resCA],        #NO:pointer to CA of a residue if == "-" CA is missing for that residues
	      bfactor => [@bfactor]};
	      #table321 => {%aa321}};

    return $self;
}

=head2 get_coord

    Title   : get_coord
    Usage   : Prints out a pdb all atoms, where the SEQRES and
              CA residues match.
    Function: 
    Example : $pdb_obj->print_pdb(*FILEHANDLE, %lookup_hash)
    Returns : Void
    Args    : A filehandle  and a lookup_hash for the numbering of the residues (Optional);

=cut

sub get_coord
{
    my ($self,$rescount)=@_;
    my $old_resname="undef";
    # my @CAali = @{$self->{resCA}};
    my $counter=0;
    my $print_res=0;
    my $j=0;
    for(my $i=0;$i<scalar @{$self->{res}};$i++)
    { 
	my $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	if($old_resname ne $current_resname)
	{
	     if($self->{resCA}[$j] ne "-")  # Just print the residue for which a CA exists
	     {
		 $counter++;
		 $print_res=1 if($counter==$rescount);
		     
	     }
	     $j++;
	 }
	 if($print_res && $self->{atomtype}[$i] eq "CA")
	{
	    #printf ("%-6s %4d  %-3s%1s%3s %1s%5s   %8.3f%8.3f%8.3f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$self->{chain},$self->{resno}[$i].$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
	    return($self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
	}
	 $old_resname=$current_resname;
    }
    return 0 if(!$print_res);
}





sub translate
{
    my ($self,$x,$y,$z)=@_;
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i] ------ ";
	$self->{x}[$i]=sprintf("%8.3f",$self->{x}[$i]+$x);
	$self->{y}[$i]=sprintf("%8.3f",$self->{y}[$i]+$y);
	$self->{z}[$i]=sprintf("%8.3f",$self->{z}[$i]+$z);

	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]\n";
    }
    return $self;
}


sub shrink
{
    my ($self,$factor)=@_;
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i] ------ ";
	$self->{x}[$i]=sprintf("%8.3f",$self->{x}[$i]*$factor);
	$self->{y}[$i]=sprintf("%8.3f",$self->{y}[$i]*$factor);
	$self->{z}[$i]=sprintf("%8.3f",$self->{z}[$i]*$factor);	
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]\n";
    }
}



sub rotate_KI
{
    my ($self,$trans)=@_;
    my @trans=@{$trans};
#    my $last_ca=0;
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	if($self->{atomtype}[$i] eq "CA")
	{
	    #print "TEST1 $self->{x}[$i] $self->{z}[$i] $self->{y}[$i]\n" if($self->{atomtype}[$i] eq "CA");
	    #print "TEST1a $trans[0][0]*$self->{x}[$i]+$trans[0][1]*$self->{y}[$i]+$trans[0][2]*$self->{z}[$i]\n";
	    #print "TEST2a $trans[1][0]*$self->{x}[$i]+$trans[1][1]*$self->{y}[$i]+$trans[1][2]*$self->{z}[$i]\n";
	    #print "TEST3a $trans[2][0]*$self->{x}[$i]+$trans[2][1]*$self->{y}[$i]+$trans[2][2]*$self->{z}[$i]\n";
	}
#	print "$trans[0][0]\n";
	my $x=$trans[0][0]*$self->{x}[$i]+$trans[0][1]*$self->{y}[$i]+$trans[0][2]*$self->{z}[$i];
	my $y=$trans[1][0]*$self->{x}[$i]+$trans[1][1]*$self->{y}[$i]+$trans[1][2]*$self->{z}[$i];
	my $z=$trans[2][0]*$self->{x}[$i]+$trans[2][1]*$self->{y}[$i]+$trans[2][2]*$self->{z}[$i];
	
	$self->{x}[$i]=$x;
	$self->{y}[$i]=$y;
	$self->{z}[$i]=$z;
#	$last_ca=$i if($self->{atomtype}[$i] eq "CA");
#	print "TEST2 $self->{x}[$i] $self->{y}[$i] $self->{z}[$i]\n" if($self->{atomtype}[$i] eq "CA");
	#print "TEST2a $x $y $z\n" if($self->{atomtype}[$i] eq "CA"); 
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i] ------ ";
	#$self->{x}[$i]=sprintf("%8.3f",$self->{x}[$i]+$x);
	#$self->{y}[$i]=sprintf("%8.3f",$self->{y}[$i]+$y);
	#$self->{z}[$i]=sprintf("%8.3f",$self->{z}[$i]+$z);
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]\n";
    }
    #print "TEST2 $self->{x}[$last_ca] $self->{z}[$last_ca] $self->{y}[$last_ca]\n";


    return $self;  #

}

sub rotate
{
    my ($self,$trans)=@_;
    my @trans=@{$trans};
#    my $last_ca=0;
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	if($self->{atomtype}[$i] eq "CA")
	{
	    #print "TEST1 $self->{x}[$i] $self->{z}[$i] $self->{y}[$i]\n" if($self->{atomtype}[$i] eq "CA");
	    #print "TEST1a $trans[0][0]*$self->{x}[$i]+$trans[0][1]*$self->{y}[$i]+$trans[0][2]*$self->{z}[$i]\n";
	    #print "TEST2a $trans[1][0]*$self->{x}[$i]+$trans[1][1]*$self->{y}[$i]+$trans[1][2]*$self->{z}[$i]\n";
	    #print "TEST3a $trans[2][0]*$self->{x}[$i]+$trans[2][1]*$self->{y}[$i]+$trans[2][2]*$self->{z}[$i]\n";
	}
#	print "$trans[0][0]\n";
	my $x=$trans[0][0]*$self->{x}[$i]+$trans[0][1]*$self->{y}[$i]+$trans[0][2]*$self->{z}[$i];
	my $y=$trans[1][0]*$self->{x}[$i]+$trans[1][1]*$self->{y}[$i]+$trans[1][2]*$self->{z}[$i];
	my $z=$trans[2][0]*$self->{x}[$i]+$trans[2][1]*$self->{y}[$i]+$trans[2][2]*$self->{z}[$i];
	
	$self->{x}[$i]=$x;
	$self->{y}[$i]=$y;
	$self->{z}[$i]=$z;
#	$last_ca=$i if($self->{atomtype}[$i] eq "CA");
#	print "TEST2 $self->{x}[$i] $self->{y}[$i] $self->{z}[$i]\n" if($self->{atomtype}[$i] eq "CA");
	#print "TEST2a $x $y $z\n" if($self->{atomtype}[$i] eq "CA"); 
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i] ------ ";
	#$self->{x}[$i]=sprintf("%8.3f",$self->{x}[$i]+$x);
	#$self->{y}[$i]=sprintf("%8.3f",$self->{y}[$i]+$y);
	#$self->{z}[$i]=sprintf("%8.3f",$self->{z}[$i]+$z);
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]\n";
    }
    #print "TEST2 $self->{x}[$last_ca] $self->{z}[$last_ca] $self->{y}[$last_ca]\n";



}

sub center
{
    my ($self)=shift;
    my $xcen=0;
    my $ycen=0;
    my $zcen=0;
    my $natoms=0;
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	if($self->{atomtype}[$i] eq "CA")
	{
	    $xcen+=$self->{x}[$i];
	    $ycen+=$self->{y}[$i];
	    $zcen+=$self->{z}[$i];
	    $natoms++;
	}
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i] ------ ";
	#$self->{x}[$i]=sprintf("%8.3f",$self->{x}[$i]+$x);
	#$self->{y}[$i]=sprintf("%8.3f",$self->{y}[$i]+$y);
	#$self->{z}[$i]=sprintf("%8.3f",$self->{z}[$i]+$z);
	#print "$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]\n";
    }
    $xcen/=$natoms;
    $ycen/=$natoms;
    $zcen/=$natoms;
   # print "$xcen $ycen $zcen\n";
    my $last_ca;
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	$self->{x}[$i]-=$xcen;
	$self->{y}[$i]-=$ycen;
	$self->{z}[$i]-=$zcen;
	$last_ca=$i if($self->{atomtype}[$i] eq "CA");
	
    }
  #  print "TEST2 $self->{x}[$last_ca] $self->{z}[$last_ca] $self->{y}[$last_ca]\n";

}

sub radius
{
    my ($self)=shift;
    my @radius=();
    for(my $i=0;$i< scalar @{$self->{x}};$i++)
    {
	if($self->{atomtype}[$i] eq "CA")
	{
	    my $radius=distance(0,0,0,$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
	    push(@radius,$radius);
	  #  printf("%-5s %5s %8.5f\n", $self->{res}[$i], "$self->{resno}[$i]$self->{insertion_code}[$i]",$radius);
	    
	    #$xcen+=$self->{x}[$i];
	    #$ycen+=$self->{y}[$i];
	    #$zcen+=$self->{z}[$i];
	    #$natoms++;
	}
    }
    return @radius;
}


sub get_chain
{
    my $self=shift;
    return $self->{chain};

}

sub get_ali
{
    my $self=shift;
    my $coord;
    my @list=$self->{z};
}

sub number_of_res
{
    my $self=shift;
    return length($self->{seqCA});
}


sub CAres
{
    my $self=shift;
    return $self->{seqCA};
}


sub CAres2
{
    my $self=shift;
    return $self->{seqCA2};
}

sub get_atomno
{
    my $self=shift;
    return @{$self->{atomno}};
}
sub get_atomtype
{
    my $self=shift;
    return @{$self->{atomtype}};
}
sub get_resno
{
    my $self=shift;
    return @{$self->{resno}};
}
sub change_Bfactor
{
    my ($self,@B)=@_;
    my $B_index=-1;
    my $old_resname="";
    my $resname="";
    for(my $i=0;$i<scalar @{$self->{bfactor}};$i++)
    { 
	$resname=$self->{resno}[$i].$self->{insertion_code}[$i];
	$B_index++ if($old_resname ne $resname);
	$self->{bfactor}[$i]=$B[$B_index];
	$old_resname=$resname;
    }
    return $self;
}
sub add_dummy_atom
{
    my $self= shift;
    
}

#sub CAresname
#{
#    my $self=shift;
 #   my @resno=();
 #   if(length($self->{seqCA_ali}) > 1)
  #  {
#	my $len=length($self->{seqCA_ali});
#	my @ali=split('',$self->{seqCA_ali});

# Function which returns a vector with all 


sub get_seqres
{
    my $self=shift;
    if(not($self->{no_seqres}))
    {
	return $self->{seqres}->seq();
    }
    else
    {
	warn "WARNING: No SEQRES in pdb, the chain specified is probably inconsistent\n";
	return 0;
    }
}

sub get_seqres_CA_ali
{
    my $self=shift;
    return ($self->{seqres_ali},$self->{seqCA_ali});
    
}



=head2 CAresname

    Title   : CAresname
    Usage   : Helper function for fix_numbering

    Function: 
    Example : @array=$pdb_obj->CAresname()
    Returns : An array with all
    Args    : None besides the object of course.

=cut

sub CAresname
{
    my $self=shift;
    my @resno=();
    #if(scalar @{$self->{resCA}} > 1)
    #{
	#my $len=length($self->{resCA});
    my @ali=@{$self->{resCA}};
    for(my $i=0;$i<=$#ali;$i++)
    {
	if($ali[$i] ne "-")
	{
	    push(@resno,$self->{resno}[$self->{resCA}[$i]].$self->{insertion_code}[$self->{resCA}[$i]]);
	    #print $self->{res}[$self->{resCA}[$i]]." ".$self->{resno}[$self->{resCA}[$i]].$self->{insertion_code}[$self->{resCA}[$i]]."\n";
	}
    }
    #}
    #else
    #{
	#my $len=scalar(@{$self->{resCA}});
	#for(my $i=0;$i<$len;$i++)
	#{
	#    push(@resno, $self->{resno}[$self->{resCA}[$i]].$self->{insertion_code}[$self->{resCA}[$i]]);
	#    print "då\n";
	#}
    #}
    #foreach my $elem (@resno)
    #{
#	print $elem."\n";
#    }
    return @resno;
}

=head2 print_pdb

    Title   : print_pdb
    Usage   : Prints out a pdb all atoms, where the SEQRES and
              CA residues match.
    Function: 
    Example : $pdb_obj->print_pdb(*FILEHANDLE, %lookup_hash)
    Returns : Void
    Args    : A filehandle  and a lookup_hash for the numbering of the residues (Optional);

=cut





sub get_CA_coord{
    my ($self,$start)=@_;    
    my @CAali = @{$self->{resCA}}; #aligned to "first atom res"

#sub number_of_res
#{
#    my $self=shift;
#    return length($self->{seqCA}->seq());
#}




    my $len = number_of_res($self);
    my $counter = 0;
    my $res_count=0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    my @header=@{$self->{header}};
    
    
    if( $start > -2 ){
	$len=1;
	my $i=0;
	my $j=0;
	my $print_res=0;
	my $counter=-1;
	my $old_resname="";
	my $current_resname="";
	
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{   
	    $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($current_resname ne $old_resname)
	    {
		if($self->{resCA}[$j] eq "-")
		{
		    $print_res=0
		}
		else
		{
		    $print_res=1;
		}
		$counter++;
	        $j++;
	    }
	    
	    if(($counter==$start) && $print_res && $self->{atomtype}[$i] eq "CA" ){
		#printf STDERR "%8.3f%8.3f%8.3f\n", $self->{x}[$i],$self->{y}[$i],$self->{z}[$i] ;
		#printf STDERR ("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		return (  $self->{x}[$i],$self->{y}[$i],$self->{z}[$i] );
	    }
	    $old_resname=$current_resname;
	}
    }
    else
    {
	warn "get_CA_coordb: error\n";
	return 0;
	
    }
    return 0;
}





















sub print_pdb
{
#    print "IN PRINT\n";
    my ($self,$filehandle,%lookup_hash)=@_;
    my $use_lookup=scalar(keys(%lookup_hash));
    #my @ali = split(//,$self->{seqCA_ali});
    my @CAali = @{$self->{resCA}};
    		    
#print $self->{seqCA_ali};
    #print length($self->{seqCA_ali});
    #print "\n\n";
    #print $self->{seqCA}->seq();
    #print "\n\n";
    my $len = number_of_res($self);
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    my @header=@{$self->{header}};
    if($use_lookup)
    {
	print $filehandle "REMARK The numbering of the residue changed to match the\nREMARK correct structure: ".scalar(localtime)."./BW\n";
    }
    #foreach my $line(@{$self->{header}})
    #print $#ali;
    #print "\n";
    #foreach my $line(@header)
    #{
#	print $filehandle $line."\n";
#    }
    
    for (my $i=0;$i<=$#CAali;$i++)
    {
	if($CAali[$i] ne "-")
	{
	    if($counter%13==0 && $counter != 0)
	    {
		
		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		#printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		$rows++;
		$residues="";
	    }
	    my $current_index=$self->{resCA}[$i];
	    #print $current_index."\n";
	    #print $i."\n";
	    $residues.=$self->{res}[$current_index]." ";
	    #print $residues."\n";
	    $counter++;
	}
    }
    if(length($residues) != 0)
    {
	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	#printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	$residues="";
	#$counter++;
    }
    
    if($use_lookup)
    {   	
	my $i=0;
	my $j=0;
	my $old_resname="undef";
	my $print_res=0;
	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{  
	    my $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    #print "$current_resname\n";
	    if($old_resname ne $current_resname)
	    {
		if($self->{resCA}[$j] eq "-")  # Just print the residue for which a CA exists
		{
		    $print_res=0;  
		}
		else
		{
		    $print_res=1;
		}
		$j++;
	    }
	    
	    if($print_res)
	    { 
		#print "$current_resname $lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]}\n";
		#print "$self->{bfactor}[$i]\n";
		if(length($self->{atomtype}[$i])==4)
		{
		    printf $filehandle ("%-6s%5d %-4s%1s%3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]},$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		else
		{
		    #printf ("%-6s%5d  %-3s%1s%3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]},$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		    #exit;
		    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]},$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		#printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%5s   %8.3f%8.3f%8.3f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]},$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
		$last_added_index=$i;
	    }
	    $old_resname=$current_resname;
	}
	printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%5s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$lookup_hash{$self->{resno}[$last_added_index].$self->{insertion_code}[$last_added_index]});
	print $filehandle "END\n";
	close($filehandle);
    }
    else
    {
	my $i=0;
	my $j=0;
	my $old_resname="undef";
	my $print_res=0;
	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{  
	    #print "$self->{res}[$i] $self->{resno}[$i]$self->{insertion_code}[$i] $old_resname\n";
	    my $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($old_resname ne $current_resname)
	    {
		#print "In here $self->{resbegin}[$j]\n";

		if($self->{resCA}[$j] eq "-")
		{
		    $print_res=0;   # Changed to 1 from 0, which makes the whole if statement obsolete,
		}                   # I can't recall why I added it in the first place. In the case
		                    # CA is missing I think it shouldn't print the res.
		else
		{
		    $print_res=1;
		}
		$j++;
	    }
	    
	    if($print_res)
	    {
		if(length($self->{atomtype}[$i])==4)
		{
		    printf $filehandle ("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		else
		{
		    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		$last_added_index=$i;
	    }
		$old_resname=$current_resname;
	}
	printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$self->{resno}[$last_added_index],$self->{insertion_code}[$last_added_index]);
	print $filehandle "END\n";
	close($filehandle);
    }
}

=head2 print_sliced_pdb

    Title   : print_sliced_pdb
    Usage   : Prints out a pdb all atoms, where the SEQRES and
              CA residues match in the interval specified.
    Function: 
    Example : $pdb_obj->print_sliced_pdb(*FILEHANDLE,$start,$end)
    Returns : Void
    Args    : A filehandle  and start and end points

=cut

sub print_sliced_pdb
{
    my ($self,$filehandle,$start,$end)=@_;
    #my @ali = split(//,$self->{seqCA_ali});
    #print $self->{seqCA_ali};
    #print length($self->{seqCA_ali});
    #print "\n\n";
    #print $self->{seqCA}->seq();
    #print "\n\n";
    my @CAali = @{$self->{resCA}}; #aligned to "first atom res"
   # print @CAali;
    my $len = number_of_res($self);
    my $counter = 0;
    my $res_count=0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    my @header=@{$self->{header}};
    if(($end-$start)<=$len || $start<=$end)
    {
	$len=$end-$start+1;
	#print $len."\n";
	print $filehandle "REMARK Only the residues from $start to $end /BW\n";
	#foreach my $line(@{$self->{header}})
	#print $#ali;
	#print "\n";
	foreach my $line(@header)
	{
	    print $filehandle $line."\n";
	}
	for (my $i=0;$i<=$#CAali;$i++)
	{
	    #print "$CAali[$i] n\n";
	    if($CAali[$i] ne "-" && $counter>=$start-1 && $counter<=$end-1)
	    {
		
		if($res_count%13==0 && $res_count != 0)
		{
		    
		    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		    #printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		    $rows++;
		    $residues="";
		}
		my $current_index=$self->{resCA}[$i];
		#print $current_index."\n";
		#print $i."\n";
		$residues.=$self->{res}[$current_index]." ";
		#print $residues."\n";
		$res_count++;
	    }
	    $counter++ if($CAali[$i] ne "-");
	}
	if(length($residues) != 0)
	{
	    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    #printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    $residues="";
	    #$counter++;
	}
	my $i=0;
	my $j=0;
	my $print_res=0;
	my $counter=-1;
	my $old_resname="";
	my $current_resname="";

	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{   
	    $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($current_resname ne $old_resname)
	    {
		if($self->{resCA}[$j] eq "-")
		{
		    $print_res=0
		}
		else
		{
		    $print_res=1;
		}
		$counter++;
	        $j++;
	    }

	    if(($counter>=$start-1 && $counter<=$end-1) && $print_res)
	    {

		if(length($self->{atomtype}[$i])==4)
		{
		    printf $filehandle ("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		else
		{
		    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		$last_added_index=$i;
	    }
	    $old_resname=$current_resname;
	}
	printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$self->{resno}[$last_added_index],$self->{insertion_code}[$last_added_index]);
	print $filehandle "END\n";
	close($filehandle);
    }
    else
    {
	warn "print_sliced_pdb: The starting point ($start) is after the end point ($end)\n" if($start>$end);
	warn "print_sliced_pdb: Specified interval $start-$end is larger than the length\n" if($end-$start>$len);
	return 0;
	    
    }
}

sub get_sliced_pdb
{
    my ($self,$start,$end)=@_;
    #my @ali = split(//,$self->{seqCA_ali});
    #print $self->{seqCA_ali};
    #print length($self->{seqCA_ali});
    #print "\n\n";
    #print $self->{seqCA}->seq();
    #print "\n\n";
    my @CAali = @{$self->{resCA}}; #aligned to "first atom res"
   # print @CAali;
    my $len = number_of_res($self);
    my $counter = 0;
    my $res_count=0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    my @header=@{$self->{header}};
    my $pdb="";
    if(($end-$start)<=$len || $start<=$end)
    {
	$len=$end-$start+1;
	
	if(length($residues) != 0)
	{
	 #  printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    #printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    $residues="";
	    #$counter++;
	}
	my $i=0;
	my $j=0;
	my $print_res=0;
	my $counter=-1;
	my $old_resname="";
	my $current_resname="";

	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{   
	    $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($current_resname ne $old_resname)
	    {
		if($self->{resCA}[$j] eq "-")
		{
		    $print_res=0
		}
		else
		{
		    $print_res=1;
		}
		$counter++;
	        $j++;
	    }

	    if(($counter>=$start-1 && $counter<=$end-1) && $print_res)
	    {

		if(length($self->{atomtype}[$i])==4)
		{
		    $pdb.=sprintf("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		else
		{
		    $pdb.=sprintf("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i],1,$self->{bfactor}[$i]);
		}
		$last_added_index=$i;
	    }
	    $old_resname=$current_resname;
	}
#	printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$self->{resno}[$last_added_index],$self->{insertion_code}[$last_added_index]);
	#print $filehandle "END\n";
	#close($filehandle);
    }
    else
    {
	warn "print_sliced_pdb: The starting point ($start) is after the end point ($end)\n" if($start>$end);
	warn "print_sliced_pdb: Specified interval $start-$end is larger than the length\n" if($end-$start>$len);
	return 0;
	    
    }
    return $pdb;
}



=head2 print_ca_pdb

    Title   : print_ca_pdb
    Usage   : Prints out a pdb with the CA atoms only, where the SEQRES and
              CA atom residues match.
    Function: 
    Example : $pdb_obj->print_pdb(*filehandle)
    Returns : Void
    Args    : A filehandle 

=cut


sub print_ca_pdb
{
    my ($self,$filehandle)=@_;
#    my @ali = split(//,$self->{seqCA_ali});
    my @ali = split(//,$self->{seqCA}->seq());
    #print @ali;
    my $len = scalar @ali;
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    
   #SEQRES
    if($self->{seqres} != 0)
    {
        for (my $i=0;$i<=$#ali;$i++)
	{
	    if($ali[$i] ne "-")
	    {
		if($counter%13==0 && $counter != 0)
		{
		    
		    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		    $rows++;
		    $residues="";
		}
		my $current_index=$self->{resCA}[$i];
		$residues.=$self->{res}[$current_index]." ";
		$counter++;
		#$residues.=$aa123{$res}." ";
		#$counter++;
	    }
	}
	if(length($residues) != 0)
	{
	    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    $residues="";
	    #$counter++;
	}
    }
    for(my $i=0;$i<scalar @ali;$i++)
    {   
	if($ali[$i] ne "-")
	{
	    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f\n", $self->{field}[$self->{resCA}[$i]],$self->{atomno}[$self->{resCA}[$i]],$self->{atomtype}[$self->{resCA}[$i]],$self->{atom_alt_loc}[$self->{resCA}[$i]],$self->{res}[$self->{resCA}[$i]],$chain,$self->{resno}[$self->{resCA}[$i]],$self->{insertion_code}[$self->{resCA}[$i]],$self->{x}[$self->{resCA}[$i]],$self->{y}[$self->{resCA}[$i]],$self->{z}[$self->{resCA}[$i]]);

	    #printf ("%-6s %4d  %-3s%1s%3s %1s%5s    %7.3f %7.3f %7.3f\n", $self->{field}[$self->{resCA}[$i]],$self->{atomno}[$self->{resCA}[$i]],$self->{atomtype}[$self->{resCA}[$i]],$self->{atom_alt_loc}[$self->{resCA}[$i]],$self->{res}[$self->{resCA}[$i]],$chain,$self->{resno}[$self->{resCA}[$i]],$self->{insertion_code}[$self->{resCA}[$i]],$self->{x}[$self->{resCA}[$i]],$self->{y}[$self->{resCA}[$i]],$self->{z}[$self->{resCA}[$i]]);
	}
    }
}
#=head2 remove_residues
#
#    Title   : remove_residues
#    Usage   : remove_residues not specified by the intervall
#    Function: 
#    Example : $pdb_obj->remove_residues($start,$end)
#    Returns : Void
#    Args    : $start and $end
#
#=cut
#sub remove_residues
#{
#    my ($self,$start,$end)=@_;
#
#}





=head2 print_ca_model

    Title   : print_ca_model
    Usage   : Prints out a pdb with the CA atoms only, where the SEQRES and
              CA atom residues match the provided alignment
    Function: 
    Example : $pdb_obj->print_ca_model(*filehandle,$pdb_ali,$target_ali)
    Returns : Void
    Args    : A filehandle,alignment
    $pdb_ali sequence has to be taken from the correct pdbfile -> that sequence must match the exact CAres from the pdb

=cut


sub print_ca_model
{
    my ($self,$filehandle,$pdb_ali,$target_ali)=@_;
    
    my @index_to_CAres = @{$self->{resCA}};
    
    my @pdb_ali = split(//,$pdb_ali);
    my @target_ali = split(//,$target_ali);
    
#    print "1: $pdb_ali\n\n2: $target_ali\n\n";
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    #my $chain=$self->{chain};
    my $chain=" ";
#count length,
    my $len = 0;
    for (my $i=0;$i<=$#pdb_ali;$i++)
    {
	$len++ if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-');
    }


#   #SEQRES
#    for (my $i=0;$i<=$#pdb_ali;$i++)
#    {
#	if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-')
#	{
#	    if($counter%13==0 && $counter != 0)
#	    {
#		
#		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#		$rows++;
#		$residues="";
#	    }
#	    $residues.=$aa123{$target_ali[$i]}." ";
#	    $counter++;
#	}
#    }
#    
#    if(length($residues) != 0)
#    {
#	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#	$residues="";
#	#$counter++;
#    }
    
    my $pdb_index=0;
    my $target_counter=0;
    for(my $target_index=0;$target_index<=$#target_ali;$target_index++)
    {
	$target_counter++ if($target_ali[$target_index] ne '-');
	#print "1: $pdb_ali[$target_index] 2:$target_ali[$target_index]\n";
	if($pdb_ali[$target_index] ne '-' && $target_ali[$target_index] ne '-')
	{
	    $counter++;
	    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f\n", $self->{field}[$index_to_CAres[$pdb_index]],$self->{atomno}[$index_to_CAres[$pdb_index]],"CA",$self->{atom_alt_loc}[$index_to_CAres[$pdb_index]],$aa123{$target_ali[$target_index]},$chain,$target_counter," ",$self->{x}[$index_to_CAres[$pdb_index]],$self->{y}[$index_to_CAres[$pdb_index]],$self->{z}[$index_to_CAres[$pdb_index]]);
	#    printf  ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f\n", $self->{field}[$index_to_CAres[$pdb_index]],$self->{atomno}[$index_to_CAres[$pdb_index]],"CA",$self->{atom_alt_loc}[$index_to_CAres[$pdb_index]],$aa123{$target_ali[$target_index]},$chain,$target_counter," ",$self->{x}[$index_to_CAres[$pdb_index]],$self->{y}[$index_to_CAres[$pdb_index]],$self->{z}[$index_to_CAres[$pdb_index]]);
	    
	}
	$pdb_index++ if($pdb_ali[$target_index] ne '-');
    }
    if($counter > 0)
    {
	printf $filehandle "TER\n"
    }
    else
    {
	printf STDERR "No residues aligned!\n";
    }
    close($filehandle);
    #print "$pdb_ali\n$target_ali\n";
    
}


=head2 print_backbone

    Title   : print_backbone
    Usage   : Prints out a pdb with the backbone (C,N,O,CA) atoms only, where the SEQRES and
              backbone atoms residues match the provided alignment
    Function: 
    Example : $pdb_obj->print_ca_model(*filehandle,$pdb_ali,$target_ali)
    Returns : Void
    Args    : A filehandle,alignment
    $pdb_ali sequence has to be taken from the correct pdbfile -> that sequence must match the exact CAres from the pdb

=cut


sub print_backbone
{
    my ($self,$filehandle,$pdb_ali,$target_ali)=@_;
    print_backbone_model($self,$filehandle,$self->{seqCA}->seq(),$self->{seqCA}->seq());
}

=head2 print_backbone_model

    Title   : print_backbone_model
    Usage   : Prints out a pdb with the backbone (C,N,O,CA) atoms only, where the SEQRES and
              backbone atoms residues match the provided alignment
    Function: 
    Example : $pdb_obj->print_ca_model(*filehandle,$pdb_ali,$target_ali)
    Returns : Void
    Args    : A filehandle,alignment
    $pdb_ali sequence has to be taken from the correct pdbfile -> that sequence must match the exact CAres from the pdb

=cut



sub print_backbone_model
{
    my ($self,$filehandle,$pdb_ali,$target_ali)=@_;
    
    my @index_to_res = @{$self->{resbegin}};
    
    my @pdb_ali = split(//,$pdb_ali);
    my @target_ali = split(//,$target_ali);
    
   # print "$pdb_ali\n\n$target_ali\n\n";
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    #my $chain=$self->{chain};
    my $chain=" ";
#count length,
    my $len = 0;
    for (my $i=0;$i<=$#pdb_ali;$i++)
    {
	$len++ if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-');
    }


#   #SEQRES
#    for (my $i=0;$i<=$#pdb_ali;$i++)
#    {
#	if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-')
#	{
#	    if($counter%13==0 && $counter != 0)
#	    {
#		
#		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#		$rows++;
#		$residues="";
#	    }
#	    $residues.=$aa123{$target_ali[$i]}." ";
#	    $counter++;
#	}
#    }
#    
#    if(length($residues) != 0)
#    {
#	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#	$residues="";
#	#$counter++;
#    }
    
    my $pdb_index=0;
    my $target_counter=0;
    for(my $target_index=0;$target_index<=$#target_ali;$target_index++)
    {
	$target_counter++ if($target_ali[$target_index] ne '-');
#	print "$pdb_ali[$target_index] $target_ali[$target_index] $pdb_index $index_to_res[$pdb_index] $aa321{$self->{res}[$index_to_res[$pdb_index]]}\n";
	if($pdb_ali[$target_index] ne '-' && $target_ali[$target_index] ne '-')
	{
	    #Print all N,C,O and CA from file.
	    my $res=$aa123{$target_ali[$target_index]};
	    my $index=$index_to_res[$pdb_index];
	    my $current_resname="$self->{resno}[$index]$self->{insertion_code}[$index]";
	    my $old_resname=$current_resname;
	    while($current_resname eq $old_resname)
	    {
	#	print "$current_resname $self->{atomno}[$index] $self->{res}[$index] $self->{atomtype}[$index] $index $pdb_index\n";
		if($self->{atomtype}[$index] eq "C" ||
		   $self->{atomtype}[$index] eq "O" ||
		   $self->{atomtype}[$index] eq "N" ||
		   $self->{atomtype}[$index] eq "CA")
		{
		    #print "$self->{atomno}[$index] $self->{atomtype}[$index] $self->{res}[$index] $res ($target_ali[$target_index])\n";
		    #print "$self->{atomno}[$index] $self->{res}[$index]\n";
		    #print "$self->{field}[$index],$self->{atomno}[$index],$self->{atomtype}[$index],$self->{atom_alt_loc}[$index],$res,$chain,$target_counter, ,$self->{x}[$index],$self->{y}[$index],$self->{z}[$index]\n";
		    my $bfactor=9.99;
		    #print "$self->{bfactor}[$index]\n";
		    $bfactor=$self->{bfactor}[$index] if($self->{bfactor}[$index]=~/\d/);
		    #print $bfactor."\n";

		    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$index],$self->{atomno}[$index],$self->{atomtype}[$index],$self->{atom_alt_loc}[$index],$res,$chain,$target_counter," ",$self->{x}[$index],$self->{y}[$index],$self->{z}[$index],1,$bfactor);
		}
		$index++;
		if(defined($self->{resno}[$index]))
		{
		    $old_resname="$self->{resno}[$index]$self->{insertion_code}[$index]";
		}
		else
		{
		    $old_resname="";
		}
		    
	    }
	    
	}
	$pdb_index++ if($pdb_ali[$target_index] ne '-');
    }
    printf $filehandle "TER\n";
    close($filehandle);
}

sub print_all_identical_atoms_model
{
    my ($self,$filehandle,$pdb_ali,$target_ali)=@_;
    
    my @index_to_res = @{$self->{resbegin}};
    
    my @pdb_ali = split(//,$pdb_ali);
    my @target_ali = split(//,$target_ali);
    
   # print "$pdb_ali\n\n$target_ali\n\n";
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    #my $chain=$self->{chain};
    my $chain=" ";
#count length,
    my $len = 0;
    for (my $i=0;$i<=$#pdb_ali;$i++)
    {
	$len++ if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-');
    }


#   #SEQRES
#    for (my $i=0;$i<=$#pdb_ali;$i++)
#    {
#	if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-')
#	{
#	    if($counter%13==0 && $counter != 0)
#	    {
#		
#		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#		$rows++;
#		$residues="";
#	    }
#	    $residues.=$aa123{$target_ali[$i]}." ";
#	    $counter++;
#	}
#    }
#    
#    if(length($residues) != 0)
#    {
#	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#	$residues="";
#	#$counter++;
#    }
    
    my $pdb_index=0;
    my $target_counter=0;
    for(my $target_index=0;$target_index<=$#target_ali;$target_index++)
    {
	$target_counter++ if($target_ali[$target_index] ne '-');
#	print "$pdb_ali[$target_index] $target_ali[$target_index] $pdb_index $index_to_res[$pdb_index] $aa321{$self->{res}[$index_to_res[$pdb_index]]}\n";
	if($pdb_ali[$target_index] ne '-' && $target_ali[$target_index] ne '-')
	{
	    #Print all N,C,O and CA from file.
	    my $res=$aa123{$target_ali[$target_index]};
	    my $index=$index_to_res[$pdb_index];
	    my $current_resname="$self->{resno}[$index]$self->{insertion_code}[$index]";
	    my $old_resname=$current_resname;
	    while($current_resname eq $old_resname)
	    {
	#	print "$current_resname $self->{atomno}[$index] $self->{res}[$index] $self->{atomtype}[$index] $index $pdb_index\n";
		if($pdb_ali[$target_index] eq $target_ali[$target_index] ||
		   defined($atomtypes{"$res $self->{atomtype}[$index]"}))
#		   $self->{atomtype}[$index] eq "C" ||
#		   $self->{atomtype}[$index] eq "O" ||
#		   $self->{atomtype}[$index] eq "N" ||
#		   $self->{atomtype}[$index] eq "CA")
		{
		    #print "$self->{atomno}[$index] $self->{atomtype}[$index] $self->{res}[$index] $res ($target_ali[$target_index])\n";
		    #print "$self->{atomno}[$index] $self->{res}[$index]\n";
		    #print "$self->{field}[$index],$self->{atomno}[$index],$self->{atomtype}[$index],$self->{atom_alt_loc}[$index],$res,$chain,$target_counter, ,$self->{x}[$index],$self->{y}[$index],$self->{z}[$index]\n";
		    my $bfactor=9.99;
		    #print "$self->{bfactor}[$index]\n";
		    $bfactor=$self->{bfactor}[$index] if($self->{bfactor}[$index]=~/\d/);
		    #print $bfactor."\n";

		    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$index],$self->{atomno}[$index],$self->{atomtype}[$index],$self->{atom_alt_loc}[$index],$res,$chain,$target_counter," ",$self->{x}[$index],$self->{y}[$index],$self->{z}[$index],1,$bfactor);
		}
		$index++;
		if(defined($self->{resno}[$index]))
		{
		    $old_resname="$self->{resno}[$index]$self->{insertion_code}[$index]";
		}
		else
		{
		    $old_resname="";
		}
		    
	    }
	    
	}
	$pdb_index++ if($pdb_ali[$target_index] ne '-');
    }
    printf $filehandle "TER\n";
    close($filehandle);
}


=head2 print_backbone_model2

    Title   : print_backbon_model
    Usage   : Prints out a pdb with the backbone (C,N,O,CA) atoms only if not the residue is the same 
              in which case it prints out the whole residue, where the SEQRES and
              backbone atoms residues match the provided alignment
    Function: 
    Example : $pdb_obj->print_ca_model(*filehandle,$pdb_ali,$target_ali)
    Returns : Void
    Args    : A filehandle,alignment
    $pdb_ali sequence has to be taken from the correct pdbfile -> that sequence must match the exact CAres from the pdb

=cut


sub print_backbone_model2
{
    my ($self,$filehandle,$pdb_ali,$target_ali)=@_;
    
    my @index_to_res = @{$self->{resbegin}};
    
    my @pdb_ali = split(//,$pdb_ali);
    my @target_ali = split(//,$target_ali);
    
    print "$pdb_ali\n\n$target_ali\n\n";
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    #my $chain=$self->{chain};
    my $chain=" ";
#count length,
    my $len = 0;
    for (my $i=0;$i<=$#pdb_ali;$i++)
    {
	$len++ if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-');
    }


#   #SEQRES
#    for (my $i=0;$i<=$#pdb_ali;$i++)
#    {
#	if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-')
#	{
#	    if($counter%13==0 && $counter != 0)
#	    {
#		
#		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#		$rows++;
#		$residues="";
#	    }
#	    $residues.=$aa123{$target_ali[$i]}." ";
#	    $counter++;
#	}
#    }
#    
#    if(length($residues) != 0)
#    {
#	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#	$residues="";
#	#$counter++;
#    }
    
    my $pdb_index=0;
    my $target_counter=0;
    for(my $target_index=0;$target_index<=$#target_ali;$target_index++)
    {
	$target_counter++ if($target_ali[$target_index] ne '-');
#	print "$pdb_ali[$target_index] $target_ali[$target_index] $pdb_index $index_to_res[$pdb_index] $aa321{$self->{res}[$index_to_res[$pdb_index]]}\n";
	if($pdb_ali[$target_index] ne '-' && $target_ali[$target_index] ne '-')
	{
	    #Print all N,C,O and CA from file.
	    my $res=$aa123{$target_ali[$target_index]};
	    my $index=$index_to_res[$pdb_index];
	    my $current_resname="$self->{resno}[$index]$self->{insertion_code}[$index]";
	    my $old_resname=$current_resname;
	    while($current_resname eq $old_resname)
	    {
	#	print "$current_resname $self->{atomno}[$index] $self->{res}[$index] $self->{atomtype}[$index] $index $pdb_index\n";
		if($self->{atomtype}[$index] eq "C" ||
		   $self->{atomtype}[$index] eq "O" ||
		   $self->{atomtype}[$index] eq "N" ||
		   $self->{atomtype}[$index] eq "CA" ||
		   $pdb_ali[$target_index] eq $target_ali[$target_index])
		{
		    #print "$self->{atomno}[$index] $self->{atomtype}[$index] $self->{res}[$index] $res ($target_ali[$target_index])\n";
		    #print "$self->{atomno}[$index] $self->{res}[$index]\n";
		    #print "$self->{field}[$index],$self->{atomno}[$index],$self->{atomtype}[$index],$self->{atom_alt_loc}[$index],$res,$chain,$target_counter, ,$self->{x}[$index],$self->{y}[$index],$self->{z}[$index]\n";
		    my $bfactor=9.99;
		    $bfactor=$self->{bfactor}[$index] if($self->{bfactor}[$index]=~/\d/);
		    printf $filehandle ("%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", $self->{field}[$index],$self->{atomno}[$index],$self->{atomtype}[$index],$self->{atom_alt_loc}[$index],$res,$chain,$target_counter," ",$self->{x}[$index],$self->{y}[$index],$self->{z}[$index],1,$bfactor);
		}
		$index++;
		if(defined($self->{resno}[$index]))
		{
		    $old_resname="$self->{resno}[$index]$self->{insertion_code}[$index]";
		}
		else
		{
		    $old_resname="";
		}
		    
	    }
	    
	}
	$pdb_index++ if($pdb_ali[$target_index] ne '-');
    }
    printf $filehandle "TER\n";
    close($filehandle);
}



=head2 fix_numbering

    Title   : fix_numbering
    Usage   : prints out a pdb file with the residue numbers (residue names) provided by the 
              correct pdb file

    Function: 
    Example : $pdb_obj->fix_numbering($correct_pdb_obj,*FILE);
    Returns : Void
    Args    : A correct pdb file, a filehandle

=cut

sub fix_numbering
{
    my ($model_obj,$correct_obj,$filehandle)=@_;
    my $model_seq=$model_obj->CAres()->seq();
    my $correct_seq=$correct_obj->CAres()->seq();
    my @resname_correct=$correct_obj->CAresname();
    my @resname_model=$model_obj->CAresname();
    
    #foreach my $name(@resname_model)
    #{
#	print $name."\n";
#    }

    my %lookup_hash=();
    #print scalar @resname_correct;
    #print "\n";
    #print length $correct_seq;
    #print "\n";
	


    my($model_seq2,$correct_seq2)=_align($model_seq,$correct_seq);
#    print "$model_seq2\n\n$correct_seq2\n\n";

    
    my @model_seq2=split('',$model_seq2);
    my @correct_seq2=split('',$correct_seq2);
    my $len=length($model_seq2);
    my $model_count=0;
    my $correct_count=0;
    my $missing_res_count=0;
    my $last_resnum=0;
    #my $special_counter=-999;
    my @temp=split(/[A-Z]/,$resname_model[$correct_count]);
    my $current_resnum=$temp[0];
    for(my $i=0;$i<$len;$i++)
    {
	if($correct_seq2[$i] ne "-" && $model_seq2[$i] ne "-")
	{
	    $lookup_hash{$resname_model[$model_count]}=$resname_correct[$correct_count];
	    my @temp=split(/[A-Z]/,$resname_correct[$correct_count]);
	    #print " $resname_model[$model_count] $resname_correct[$correct_count]\n";
	    $current_resnum=$temp[0];
	    $last_resnum=$current_resnum;
	    $model_count++;
	    $correct_count++;
	    #print "$resname_correct[$correct_count] $current_resnum\n";
	    #   print "$correct_seq2[$i] $model_seq2[$i] $resname_correct[$correct_count-1]\n";
	}
	elsif($correct_seq2[$i] ne "-")  # $model_seq2[$i] is not eq "-" because of first if.
	{
	  #  print $correct_seq2[$i]."\n";
	    
	    my @temp=split(/[A-Z]/,$resname_correct[$correct_count]);
	    $current_resnum=$temp[0];
	    $last_resnum=$current_resnum;
	    $correct_count++;
		
	}
	elsif($model_seq2[$i] ne "-")  # $correct_seq2[$i] is not eq "-" because of first if.
	{
	    $last_resnum++;
	    $missing_res_count=$last_resnum;
	    print "$missing_res_count"."X $resname_model[$model_count]\n";
	    $lookup_hash{$resname_model[$model_count]}="$missing_res_count"."X";
	    $model_count++;
	}
	
    }
    #foreach my $key (sort keys(%lookup_hash))
    #{
#	print "$key $lookup_hash{$key}\n";
#    }

    $model_obj->print_pdb($filehandle,%lookup_hash);

}

# sub _align   # Takes two strings removes all dashes and returns the alignment. 
# {
#     my ($seq1,$seq2)=@_;
#     $seq1=_remove_dashes($seq1);
#     $seq2=_remove_dashes($seq2);
#     #print $seq1."\n";
#     #print $seq2."\n\n";
#     my ($ali_return1,$ali_return2);
# #    my $factory=new Bio::Tools::pSW('-matrix' => '/afs/pdc.kth.se/home/k/kriil/foo3/trunk/Bio/Ext/Align/blosum62.bla','-gap' => 1,'-ext' => 0);
#     my $factory=new Bio::Tools::pSW('-matrix' => '/afs/pdc.kth.se/home/b/bjornw/bjorn/bioperl/bioperl-live/examples/blosum62.bla','-gap' => 1,'-ext' => 0);
#     my $seq_obj1=Bio::Seq->new(-moltype => "protein", -seq => $seq1, -id => "seq1");
#     my $seq_obj2=Bio::Seq->new(-moltype => "protein", -seq => $seq2, -id => "seq2");

#     #print "S:$seq1\n";
#     #xxx
# #    if ( length $seq1 == 1 and ( $seq1 eq "A" or $seq1 eq "G" or $seq1 eq "C" or $seq1 eq "T" )){
# 	#return ($seq1, $seq2 );
# #	print "tut\n";
# #	return ( $seq_obj1, $seq_obj2 );
# 	#$seq_obj1=Bio::Seq->new( -seq => $seq1, -id => "seq1");
# 	#$seq_obj2=Bio::Seq->new( -seq => $seq2, -id => "seq2");
# 	#$seq_obj1=Bio::Seq->new(-moltype => "dna", -seq => $seq1, -id => "seq1");
# 	#$seq_obj2=Bio::Seq->new(-moltype => "dna", -seq => $seq2, -id => "seq2");
# #    }
#     #print "aligning...\n";
    
#     my $aln = $factory->pairwise_alignment($seq_obj1,$seq_obj2);

#     #print "fixing...\n";

#   #  $factory->align_and_show($seq_obj1,$seq_obj2,*STDOUT);
#     ($ali_return1,$ali_return2)=_fix_alignment($aln,$seq_obj1,$seq_obj2);
#     #print "L:$ali_return1\n";
#     #print "L:$ali_return2\n";
#     return ($ali_return1,$ali_return2);
# }


sub distance
{
    my ($x1,$y1,$z1,$x2,$y2,$z2)=@_;
    #print "x1=$x1 y1=$y1 z1=$z1 x2=$x2 y2=$y2 z2=$z2\n";
    #print "\n";
    return sqrt(($x1-$x2)*($x1-$x2) + ($y1-$y2)*($y1-$y2) + ($z1-$z2)*($z1-$z2));
}

sub _remove_dashes
{
    my $temp=shift;
    my @list=split(/\-+/,$temp);
    my $new_str=join('',@list);
    return $new_str;

}



sub _fix_alignment
{
    my($aln,$seq1,$seq2)=@_;
    
    # Parse the alignment
    my ($seq,@start, @end,@ali_seqs);
    my $i=0;

    foreach $seq ($aln->each_seq())
    {
	$start[$i]=$seq->start();
	$end[$i]=$seq->end();
	$ali_seqs[$i]=$seq->seq();
	$i++;
    }
    #print $seq1->seq()."\n\n";
    #print $seq2->seq()."\n\n";
    #print $start[0]."\n".$start[1]."\n";
# Reformat alignment so that it contain all resides and -.
    
    my $ali_seq1="";
    my $ali_seq2="";
    
# Fix the begining
    if($start[0]!=1)
    {
	$ali_seq1.=$seq1->subseq(1,$start[0]-1);
	$ali_seq2.=_dashes($start[0]-1);
    }
    if($start[1]!=1)
    {
	$ali_seq1.=_dashes($start[1]-1);
	$ali_seq2.=$seq2->subseq(1,$start[1]-1);
    }
    # Add the alignment
    
    $ali_seq1.=$ali_seqs[0];
    $ali_seq2.=$ali_seqs[1];
    
# Fix alignment end;
    if($end[0]<$seq1->length())
    {
	my $len=$seq1->length();
	$ali_seq1.=$seq1->subseq($end[0]+1,$len);
	$ali_seq2.=_dashes($len-$end[0]);
    }
    
    if($end[1]<$seq2->length())
    {
	my $len=$seq2->length();
	$ali_seq1.=_dashes($len-$end[1]);
	$ali_seq2.=$seq2->subseq($end[1]+1,$len);
    }
    #print $ali_seq1."\n\n".$ali_seq2."\n\n";
    return ($ali_seq1,$ali_seq2);
}

sub _dashes
{
    my $number=shift;
    my $str="";

    for(my $i=0;$i<$number;$i++)
    {
	$str.="-";
    }
    return $str;

}
1;



