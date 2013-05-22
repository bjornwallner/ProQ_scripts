#! /usr/bin/perl -w
use strict;
use Math::Trig;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Env qw(LD_LIBRARY_PATH);
#use vars qw();
#my $FASTAfile=$ARGV[0];
my $INSTALL_DIR=dirname(abs_path($0));
#print $INSTALL_DIR."\n";
#exit;



if(scalar(@ARGV) % 2 !=0 || scalar(@ARGV) == 0){
    print "\n|||-----------------------------------------------------------|||\n";
    print "|||  ZPRED - Predicting the distance to the membrane center   |||\n";
    print "|||-----------------------------------------------------------|||\n";
#     print "Usage: zpred -mode [fasta/profile/profile_topology] -sequence [fastafile] -profile [optional .psi file] -topology [topology file] -prediction [tmhmm/modhmm] -out [name of output file]\n\n";
    print "Usage: zpred -mode [fasta/profile/profile_topology] -sequence [fastafile] -profile [optional .psi file] -topology [topology file] -prediction [tmhmm/modhmm] -out [name of output file] -tmpdir [name of dir to put temporary file]\n\n";
    #maybe want to add plotting?
    exit(1);
}

my %FLAGS=@ARGV;
my ($mode,$fastafile,$profilefile,$outputfile,$topologyfile,$topology,$prediction_method);

# Base directories for auxiliary files
# my $zpreddir = "/afs/pdc.kth.se/home/a/aronh/0vol/embrace/zpred/zpred-src";
# my $tmhmmdir = "/afs/pdc.kth.se/home/a/aronh/prodiv_tmhmm_v0.92b";
# my $blastdir = "/afs/pdc.kth.se/home/a/aronh/1vol/blast-2.2.16";

#my $zpreddir = "/afs/pdc.kth.se/home/k/kriil/vol_03/Programs/zpred/zpred-src";
my $zpreddir = "$INSTALL_DIR/"; #afs/pdc.kth.se/home/k/kriil/vol_03/Programs/zpred/zpred-src";
my $tmhmmdir = "/afs/pdc.kth.se/home/k/kriil/vol_03/Programs/topology_predictors/prodiv_tmhmm";
my $blastdir = "/afs/pdc.kth.se/home/a/arnee/MODULES/seqtools/blast-2.2.18_x86_64";
my $modhmmdir = "/afs/pdc.kth.se/home/k/kriil/vol_03/Programs/modhmm";

if(defined($FLAGS{-mode})){
    $mode=$FLAGS{-mode};
    if(($mode ne "fasta") and ($mode ne "profile") and ($mode ne "profile_topology")){
	print "Please specify fasta, profile or profile_topology\n";
	exit(0);
    }
}

if($mode eq "fasta"){
    if(defined($FLAGS{-sequence})){
	$fastafile=$FLAGS{-sequence};
	if(!-e $fastafile){
	    print "Specified fastafile: \"$fastafile\" does not exist\n";
	    exit(0);
	}
    }
}
if($mode eq "profile"){
    if(defined($FLAGS{-sequence})){
	$fastafile=$FLAGS{-sequence};
    }
    elsif(defined($FLAGS{-profile})){
	$profilefile=$FLAGS{-profile};
    }
    else{
	print "Please provide fastasequence or profile (.psi) file\n";
	exit(0);
    }
}
if($mode eq "profile_topology"){
    if(defined($FLAGS{-sequence})){
	$fastafile=$FLAGS{-sequence};
    }
    elsif(defined($FLAGS{-profile})){
	$profilefile=$FLAGS{-profile};
    }
    else{
	print "Please provide fastasequence or profile (.psi) file\n";
	exit(0);
    }
}

my $tmpdir;
if ( defined($FLAGS{-tmpdir}) ) {
    $tmpdir = $FLAGS{-tmpdir};
}
else {
    print "Please provide a directory for temporary files\n";
    exit(0);
}

my $stamp=$$;
if(defined($FLAGS{-out})){
    $outputfile=$FLAGS{-out};
}
else{
    $outputfile=$stamp.".temp";
}


my $current_dir=getcwd;
chdir "$tmpdir" or `mkdir $tmpdir`;
chdir $current_dir;
#check to see whether erikgr exists, creates it if it doesn't...


my $winsize=9; #DO NOT CHANGE WINSIZE!
my (@NET1,@NET2,@NET3,@NET4,@NET5,@MEAN_NET);
my @SEQUENCE;

###for profiles and tmhmm
##############################
if($mode eq "profile_topology"){
    my $profile;

    if(defined($FLAGS{-sequence})) {
    #if(-e $fastafile){
	run_PSIBLAST($fastafile);
	
	my @fields=split(/\//,$fastafile);

	###This part is not good!!!;
	if(scalar(@fields) == 1){
	    $profile="$tmpdir/".$fastafile.".psi";
	}
	else{
	    $profile="$tmpdir/".$fields[scalar(@fields)-1].".psi";
	}
    }
    elsif(-e $profilefile){
	$profile=$profilefile;
    }
    else{
	print "no sequence or profile provided...\n";
	exit(0);
    }
    my %PSIblasthash;
    
    my $net1file="$zpreddir/data/net1_tmhmm_profile";
    my $net2file="$zpreddir/data/net2_tmhmm_profile";
    my $net3file="$zpreddir/data/net3_tmhmm_profile";
    my $net4file="$zpreddir/data/net4_tmhmm_profile";
    my $net5file="$zpreddir/data/net5_tmhmm_profile";
    
    
    my @PSIBLASTarray;
    for(my $j=0;$j<$winsize;$j++){
	push @PSIBLASTarray,"0;-1;-1;-1;-2;-1;-1;-1;-1;-1;-1;-1;-1;-1;-2;0;0;-2;-1;-1;";
    }
    
    
    %PSIblasthash=read_psiblast_profile($profile);
    foreach my $hashkey (sort numerically keys %PSIblasthash){
#    print $hashkey,"\t";
#    print $PSIblasthash{$hashkey}{residue},"\t";
	push @PSIBLASTarray, $PSIblasthash{$hashkey}{A}.";".$PSIblasthash{$hashkey}{R}.";".$PSIblasthash{$hashkey}{N}.";".$PSIblasthash{$hashkey}{D}.";".$PSIblasthash{$hashkey}{C}.";".$PSIblasthash{$hashkey}{Q}.";".$PSIblasthash{$hashkey}{E}.";".$PSIblasthash{$hashkey}{G}.";".$PSIblasthash{$hashkey}{H}.";".$PSIblasthash{$hashkey}{I}.";".$PSIblasthash{$hashkey}{L}.";".$PSIblasthash{$hashkey}{K}.";".$PSIblasthash{$hashkey}{M}.";".$PSIblasthash{$hashkey}{F}.";".$PSIblasthash{$hashkey}{P}.";".$PSIblasthash{$hashkey}{S}.";".$PSIblasthash{$hashkey}{T}.";".$PSIblasthash{$hashkey}{W}.";".$PSIblasthash{$hashkey}{Y}.";".$PSIblasthash{$hashkey}{V}.";";
	push @SEQUENCE, $PSIblasthash{$hashkey}{residue};
    }
    
    for(my $j=0;$j<$winsize;$j++){
	push @PSIBLASTarray,"0;-1;-1;-1;-2;-1;-1;-1;-1;-1;-1;-1;-1;-1;-2;0;0;-2;-1;-1;";
    }

######READ OR DO THE TOPOLOGY...
    my @TEMPtopo;
    
    if(defined($FLAGS{-topology})){
	##read the topology check whether it is prediction or topology
	$topologyfile=$FLAGS{-topology};
	if(!-e $topologyfile){
	    print "Specified topologyfile: \"$topologyfile\" does not exist\n";
	    exit(0);
	}
	else{
	    $topology=read_topology_file($topologyfile);
	    @TEMPtopo=convert_topology($topology);
	}
    }
    else{
	##do a topology (use TMHMM at the moment..)
	if(defined($FLAGS{-prediction})){
	    $prediction_method=$FLAGS{-prediction};
	}
	else{
	    print "choose known topology or prediction\n";
	    exit(0);
	}
	
	if($prediction_method eq "tmhmm"){
	    #TMHMM not possible at the moment.
	}
	elsif($prediction_method eq "modhmm"){
	    run_MODHMM($fastafile);
	    my $sname = `basename $fastafile`; chomp $sname;
	    my $modhmmfile=$tmpdir."/".$sname.".modhmm";
	    $topology=parse_PRODIV_results($modhmmfile);
	    @TEMPtopo=convert_topology($topology);
	}
	else{
	    print "Please specify prediction method tmhmm or modhmm, not: ",$prediction_method,"\n";
	    exit(0);
	}
    }
    
    
    
######END TOPOLOGY STUFF...

    my @temparray;
    my $tempinput;

    if(scalar(@TEMPtopo) != scalar(@PSIBLASTarray)){
	print "length of topology and length of sequence/profile are not the same\n";
	exit(0);
    }
    for(my $i=0;$i<scalar(@PSIBLASTarray)-18;$i++){
	$tempinput="";
	
	for(my $j=0;$j<19;$j++){
	    $tempinput.=$TEMPtopo[$j+$i];
	}
	for(my $j=0;$j<19;$j++){
	    @temparray=split(/\;/,$PSIBLASTarray[$j+$i]);
	    foreach my $element (@temparray){
		$tempinput.=logistic($element)." ";
	    }
	}
	
#    print $tempinput,"\n";
	push @NET1, NN_fwd($tempinput,$net1file);
	push @NET2, NN_fwd($tempinput,$net2file);
	push @NET3, NN_fwd($tempinput,$net3file);
	push @NET4, NN_fwd($tempinput,$net4file);
	push @NET5, NN_fwd($tempinput,$net5file);
    }
}

###for profiles...
##############################
if($mode eq "profile"){
    my $profile;
    if(-e $fastafile){
#	my $PSIBLASTscript="/afs/pdc.kth.se/home/e/erikgr/yrojekt/WWW/run_PSIBLAST_ARGV.pl";
	run_PSIBLAST($fastafile);
#	print @PSIBLASTtext;

	my @fields=split(/\//,$fastafile);
	if(scalar(@fields) == 1){
	    $profile="$tmpdir/".$fastafile.".psi";
	}
	else{
	    $profile="$tmpdir/".$fields[scalar(@fields)-1].".psi";
	}
    }
    elsif(-e $profilefile){
	$profile=$profilefile;
    }
    else{
	print "no sequence or profile provided...\n";
	exit(0);
    }
    my %PSIblasthash;

    my $net1file="$zpreddir/data/net1_profile";
    my $net2file="$zpreddir/data/net2_profile";
    my $net3file="$zpreddir/data/net3_profile";
    my $net4file="$zpreddir/data/net4_profile";
    my $net5file="$zpreddir/data/net5_profile";

    
    my @PSIBLASTarray;
    for(my $j=0;$j<$winsize;$j++){
	push @PSIBLASTarray,"0;-1;-1;-1;-2;-1;-1;-1;-1;-1;-1;-1;-1;-1;-2;0;0;-2;-1;-1;";
    }


    %PSIblasthash=read_psiblast_profile($profile);
    foreach my $hashkey (sort numerically keys %PSIblasthash){
#    print $hashkey,"\t";
#    print $PSIblasthash{$hashkey}{residue},"\t";
	push @PSIBLASTarray, $PSIblasthash{$hashkey}{A}.";".$PSIblasthash{$hashkey}{R}.";".$PSIblasthash{$hashkey}{N}.";".$PSIblasthash{$hashkey}{D}.";".$PSIblasthash{$hashkey}{C}.";".$PSIblasthash{$hashkey}{Q}.";".$PSIblasthash{$hashkey}{E}.";".$PSIblasthash{$hashkey}{G}.";".$PSIblasthash{$hashkey}{H}.";".$PSIblasthash{$hashkey}{I}.";".$PSIblasthash{$hashkey}{L}.";".$PSIblasthash{$hashkey}{K}.";".$PSIblasthash{$hashkey}{M}.";".$PSIblasthash{$hashkey}{F}.";".$PSIblasthash{$hashkey}{P}.";".$PSIblasthash{$hashkey}{S}.";".$PSIblasthash{$hashkey}{T}.";".$PSIblasthash{$hashkey}{W}.";".$PSIblasthash{$hashkey}{Y}.";".$PSIblasthash{$hashkey}{V}.";";
	push @SEQUENCE, $PSIblasthash{$hashkey}{residue};
    }
    
    for(my $j=0;$j<$winsize;$j++){
	push @PSIBLASTarray,"0;-1;-1;-1;-2;-1;-1;-1;-1;-1;-1;-1;-1;-1;-2;0;0;-2;-1;-1;";
    }
    
    my @temparray;
    my $tempinput;
    for(my $i=0;$i<scalar(@PSIBLASTarray)-18;$i++){
	$tempinput="";
	
	for(my $j=0;$j<19;$j++){
	    @temparray=split(/\;/,$PSIBLASTarray[$j+$i]);
	    foreach my $element (@temparray){
		$tempinput.=logistic($element)." ";
	    }
	}
#    print $tempinput,"\n";
	push @NET1, NN_fwd($tempinput,$net1file);
	push @NET2, NN_fwd($tempinput,$net2file);
	push @NET3, NN_fwd($tempinput,$net3file);
	push @NET4, NN_fwd($tempinput,$net4file);
	push @NET5, NN_fwd($tempinput,$net5file);
    }
}







if($mode eq "fasta"){
#######for fasta...
#####################################################
    my $net1file="$zpreddir/data/net1_fasta";
    my $net2file="$zpreddir/data/net2_fasta";
    my $net3file="$zpreddir/data/net3_fasta";
    my $net4file="$zpreddir/data/net4_fasta";
    my $net5file="$zpreddir/data/net5_fasta";

    my $fastasequence="";
    open(FILE,$fastafile) or die "can't open fastafile: $!\n$fastafile\n";
    while(<FILE>){
	chomp;
	if(/\>/){}
	else{
	    $fastasequence.=$_;
	    for(my $i=0;$i<length($_);$i++){
		push @SEQUENCE, substr($_,$i,1);
	    }
	}
    }
    close FILE or die;
    my @INDATA=read_in_sequence($fastasequence);
    
    foreach my $sparse_fasta (@INDATA){
	push @NET1, NN_fwd($sparse_fasta,$net1file);
	push @NET2, NN_fwd($sparse_fasta,$net2file);
	push @NET3, NN_fwd($sparse_fasta,$net3file);
	push @NET4, NN_fwd($sparse_fasta,$net4file);
	push @NET5, NN_fwd($sparse_fasta,$net5file);
    }
}
######################################################

###parse differently whether profile_tmhmm or fasta/profile is used
###draw the topology...

open(OUTFILE,">$outputfile") or die "can't open outputfile: $!\n$outputfile\n";
for(my $i=0;$i<scalar(@NET1);$i++){
    my $temp;
    $temp=($NET1[$i]+$NET2[$i]+$NET3[$i]+$NET4[$i]+$NET5[$i])/5;
    push @MEAN_NET, $temp;
    print OUTFILE $temp,"\t";
    print OUTFILE $NET1[$i],"\t",$NET2[$i],"\t",$NET3[$i],"\t",$NET4[$i],"\t",$NET5[$i],"\t",$SEQUENCE[$i],"\n";
}
close OUTFILE or die;

### Plot stuff: we don't need that for webservices
# my $gnufile="";
# if($mode eq "profile_topology"){
#     my $gnutext=make_gnuplot_topology($outputfile,$topology);
#     $gnufile=$outputfile.".gnu";
#     open(GNUFILE,">$gnufile") or die "can't open gnufile: $!\n$gnufile\n";
#     print GNUFILE $gnutext;
#     close GNUFILE or die;
# }
# else{
#     my $gnutext=make_gnuplot($outputfile);
#     $gnufile=$outputfile.".gnu";
#     open(GNUFILE,">$gnufile") or die "can't open gnufile: $!\n$gnufile\n";
#     print GNUFILE $gnutext;
#     close GNUFILE or die;
# }

# `gnuplot $gnufile`;
# my $from=$outputfile.".eps";
# my $temp=$outputfile."2.eps";
# my $to=$outputfile.".png";

# `epsffit 0 0 600 600 $from $temp`;
# `convert $temp $to`;
# `cp $to /var/www/html/cbr/zpred/images/.`;

# # copy datafile to webserver
# `cp $outputfile /var/www/html/cbr/zpred/datafiles/.`;

###--------------------
###     SUBROUTINES
###--------------------

sub read_in_sequence{
    my (@FASTADATA,@SPARSEDATA);
    my @winsize=9;

    my ($fastasequence)=@_;
    for(my $i=1;$i<=length($fastasequence);$i++){
	push @FASTADATA, get_win($winsize,$winsize,$i,$fastasequence);
    }
    foreach my $element (@FASTADATA){
	push @SPARSEDATA, convert_to_sparse_encoding($element);
    }
    return(@SPARSEDATA);
}


sub get_win {
    my ($output,$tempseq);
    my ($leftwin, $rightwin, $position, $seq) = @_;

    for(my $i=0;$i<=$leftwin;$i++){
        $tempseq.="X";
    }
    $seq=$tempseq.$seq;
    for(my $i=0;$i<=$rightwin;$i++){
        $seq.='X';
    }
    $output=substr($seq, $position, $leftwin+$rightwin+1);
    return $output;
}

sub convert_to_sparse_encoding{
  my $output;
  my ($seq)=@_;
  my %aaspe =('A', "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'R', "0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'N', "0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'D', "0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'C', "0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'Q', "0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'E', "0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'G', "0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'H', "0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ",
              'I', "0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ",
              'L', "0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ",
              'K', "0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ",
              'M', "0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ",
              'F', "0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ",
              'P', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ",
              'S', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ",
              'T', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ",
              'W', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ",
              'Y', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ",
              'V', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ",
              '-', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
              'X', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
             );
  for(my $i=0;$i<length($seq);$i++){
      $output.=$aaspe{substr($seq,$i,1)};
  }
#  $output.="\n";
  return $output;
}

sub numerically{$a<=>$b;}

sub logistic{
    my ($input)=@_;
    my $tempvalue;
    $tempvalue=(1/(1+exp(-$input)));
    return $tempvalue;
}

sub read_topology_file{
    my $TOPOLOGY="";
    my ($file)=@_;
    open(FILE,$file) or die "can't open file: $!\n$file\n";
    while(<FILE>){
	chomp;
	if(/^o/ or /^M/ or /^i/ or /^O/ or /^W/){
	    $TOPOLOGY.=$_;
	}
    }
    close FILE or die;
    return($TOPOLOGY);
}
sub convert_topology{
    my @TOPOarray;
    my ($topology)=@_;
    
    for(my $j=0;$j<9;$j++){ #The winsize...
	push @TOPOarray, "0 0 ";
    }
    for(my $i=0;$i<length($topology);$i++){
	if(substr($topology,$i,1) eq "i"){
	    push @TOPOarray, "1 0 ";
	}
	elsif((substr($topology,$i,1) eq "o") or (substr($topology,$i,1) eq "O")){
	    push @TOPOarray, "1 0 ";
	}
	elsif((substr($topology,$i,1) eq "M") or (substr($topology,$i,1) eq "W")){
	    push @TOPOarray, "0 1 ";
	}
	else{
	    print "Strange characters in topology file: ",substr($topology,$i,1),"\n";
	    exit(0);
	}
    }
    for(my $j=0;$j<9;$j++){ #The winsize...
	push @TOPOarray, "0 0 ";
    }
    return(@TOPOarray);
}

sub read_psiblast_profile{
    my ($file)=@_;
    my %PSIhash;
    my @fields;

    open(FILE,$file) or die "can't open file: $!\n$file\n";
    while(<FILE>){
	chomp;
	if(/^Last/ or /^Standard/ or /^PSI/ or /^           A/ or /^                      K/ or /^\n/){}
	elsif(/[0-9]/){
	    @fields=split;
#	    print $fields[0];
	    $PSIhash{$fields[0]}={
		residue => $fields[1],
		A => $fields[2],
		R => $fields[3],
		N => $fields[4],
		D => $fields[5],
		C => $fields[6],
		Q => $fields[7],
		E => $fields[8],
		G => $fields[9],
		H => $fields[10],
		I => $fields[11],
		L => $fields[12],
		K => $fields[13],
		M => $fields[14],
		F => $fields[15],
		P => $fields[16],
		S => $fields[17],
		T => $fields[18],
		W => $fields[19],
		Y => $fields[20],
		V => $fields[21],
		A_weighted => $fields[22],
		R_weighted => $fields[23],
		N_weighted => $fields[24],
		D_weighted => $fields[25],
		C_weighted => $fields[26],
		Q_weighted => $fields[27],
		E_weighted => $fields[28],
		G_weighted => $fields[29],
		H_weighted => $fields[30],
		I_weighted => $fields[31],
		L_weighted => $fields[32],
		K_weighted => $fields[33],
		M_weighted => $fields[34],
		F_weighted => $fields[35],
		P_weighted => $fields[36],
		S_weighted => $fields[37],
		T_weighted => $fields[38],
		W_weighted => $fields[39],
		Y_weighted => $fields[40],
		V_weighted => $fields[41],
		information => $fields[42],
		relative_weight => $fields[43],
	    };
	}
    }
    close FILE or die;
    return(%PSIhash);
}

sub make_gnuplot{
    my ($gnutext,$eps);
    my ($datafile)=@_;
    
    $gnutext="set encoding iso_8859_1\n";
    $gnutext.="set yrange [0:35]\n";
    $gnutext.="set ter pos sol col enh eps\n";
    $gnutext.="set xlabel \'Residue number\'\n";
    $gnutext.="set ylabel \'Absolute value of Z-coordinate ({\305})\'\n";
    $gnutext.="set out \'".$datafile.".eps\'\n";
    $gnutext.="set ytics \(\"0\" 0, \"5\" 5, \"10\" 10, \"15\" 15, \"20\" 20, \"25\" 25\)\n";
    
    $gnutext.="plot \'".$datafile."\' u \(\$0\+1\)\:1 w li lt 2 lw 3 t \'\'\n";
    $gnutext.="exit\n";
    return($gnutext);
}
sub make_gnuplot_topology{
    my ($gnutext,$eps);
    my ($datafile,$topology)=@_;
    my @fields;
    my @BOUNDARIES=get_boundaries_from_topology($topology);

    $gnutext="set encoding iso_8859_1\n";
    for(my $i=0;$i<scalar(@BOUNDARIES);$i++){
	@fields=split(/\t/,$BOUNDARIES[$i]);
	my $start=$fields[0]+1;
	my $end=$fields[1]+1;
	if($fields[2] eq "M"){
	    $gnutext.="set arrow from ".$start.",33 to ".$end.",33 nohead lt 1 lw 40\n";
	}
	if($fields[2] eq "i"){
	    $gnutext.="set arrow from ".$start.",32.4 to ".$end.",32.4 nohead lt 3 lw 10\n";
	}
	if($fields[2] eq "o"){
	    $gnutext.="set arrow from ".$start.",33.6 to ".$end.",33.6 nohead lt 4 lw 10\n";
	}
    }
    $gnutext.="set yrange [0:35]\n";
    $gnutext.="set ter pos sol col enh eps\n";
    $gnutext.="set xlabel \'Residue number\'\n";
    $gnutext.="set ylabel \'Absolute value of Z-coordinate ({\305})\'\n";
    $gnutext.="set out \'".$datafile.".eps\'\n";
    $gnutext.="set ytics \(\"0\" 0, \"5\" 5, \"10\" 10, \"15\" 15, \"20\" 20, \"25\" 25, \"Topology\" 33\)\n";
    
    $gnutext.="plot \'".$datafile."\' u \(\$0\+1\)\:1 w li lt 2 lw 3 t \'\'\n";
    $gnutext.="exit\n";
    return($gnutext);
}

sub get_boundaries_from_topology{
    my @BOUNDARIES;
    my ($topology)=@_;
    my $temp;
    
    while($topology=~m/i+/g){
	$temp=$-[0]."\t".$+[0]."\ti";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/o+/g){
	$temp=$-[0]."\t".$+[0]."\to";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/O+/g){
	$temp=$-[0]."\t".$+[0]."\to";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/M+/g){
	$temp=$-[0]."\t".$+[0]."\tM";
	push @BOUNDARIES, $temp;
    }
    return(@BOUNDARIES);
}
sub run_MODHMM{
    my ($fastafile)=@_;
    
    my $e_value_cutoff=1e-05;
    my $Fflag="F";
    my $MSAdir="$tmpdir/";
    my $MODdir="$tmpdir/";
    my $RESULTdir="$tmpdir/";
    my $SEQdir="$tmpdir/";

#     my $seqlist=$fastafile.".list";
#     my $tempfile=$fastafile.".fa";
#     my $modfile=$fastafile.".mod";
#     my $ucmodfile=$fastafile.".ucmod";

    my $shortname=`basename $fastafile`;
    chomp $shortname;

    my $seqlist="$tmpdir/$shortname.list";
    my $tempfile="$tmpdir/$shortname.fa";
    my $modfile="$tmpdir/$shortname.mod";
    my $ucmodfile="$tmpdir/$shortname.ucmod";
    
    my $fastasequence="";
    open(FASTA,$fastafile) or die "can't open fastafile: $!\n$fastafile\n";
    while(<FASTA>){
	chomp;
	if(/\>/){}
	else{
	    $fastasequence.=$_;
	}
    }
    close FASTA or die;

    open(TEMP,">$tempfile") or die "can't open tempfile for output: $!\n$tempfile\n";
    print TEMP ">",$shortname,"\n";
    print TEMP $fastasequence,"\n";
    close TEMP or die;

    open(SEQLIST,">$seqlist") or die "can't open sequence list for output: $!\n";
    print SEQLIST $shortname,"\n";
    close SEQLIST or die;
    
    my $blast_msa_command="$zpreddir/bin/modhmm/bin/blast_msa.pl";
    my $msa62mod_command="$zpreddir/bin/modhmm/bin/msa62mod.pl";   #msa62 has different output than msa2mod
    my $mod_upper_caser_command="$zpreddir/bin/modhmm/bin/mod_upper_caser.pl";
    my $hmm_command="$modhmmdir/bin/modhmms -m $tmhmmdir/HMM_FILES/PRODIV_TMHMM_0.92b.txt -f msa  -r $tmhmmdir/util/replacement_letter_multi.rpl -M GM -L -c 1 --max_d";
#    print "$res";
    
    
    my ($seqnames,$tempvariable,$tempvariable2,$name);
#    print "Blasting....\n";
    `$blast_msa_command $seqlist $e_value_cutoff $Fflag $MSAdir $SEQdir`;
#    print "Modding...\n";
    `$msa62mod_command $seqlist 0 $MSAdir $MODdir $SEQdir`;

#    print "Upper casing...\n";
    open(MOD,$modfile) or die "can't open modfile: $!\n$modfile\n";
    open(UCMOD,">$ucmodfile") or die "can't open ucmodfile: $!\n$ucmodfile\n";
    while(<MOD>){
	if($_ =~ '/') {
	}
	else {
	    $_ =~ tr/a-z/A-Z/;
	    print UCMOD "$_";
	}
    }
    close MOD or die;
    close UCMOD or die;

    `mv $ucmodfile $modfile`;
    
    open(SEQLIST,">$seqlist") or die "can't open sequence list for output: $!\n";
    print SEQLIST $modfile,"\n";
    close SEQLIST or die;

    my $TMHMM_spec_file="$zpreddir/bin/modhmm/bin/TMHMM2.1_msa_5_sj_best.txt";
    my $replacement_letter="$zpreddir/bin/modhmm/bin/util/replacement_letters.rpl";
#    print "Hmming...\n";
    my @ERROR=`$hmm_command -s $seqlist -o $tmpdir/ > $tmpdir/modhmmout.xml`;
#    print @ERROR;
    `$modhmmdir/bin/modhmmxml2res < $tmpdir/modhmmout.xml > $tmpdir/PRODIV_TMHMM_0.92b.hmg.res`;
    my $outfile=$modfile."hmm";
    `mv $tmpdir/PRODIV_TMHMM_0.92b.hmg.res $outfile`;
}
sub parse_MODHMM_results{
    my $topology="";
    my ($modhmmfile)=@_;

    my $flag=0;
    open(MODHMM,$modhmmfile) or die "can't open modhmmfile: $!\n$modhmmfile\n";
    while(<MODHMM>){
	chomp;
	if((/^o/ or /^i/ or /^O/) and $flag==0){
	    $flag=1;
	}
	if($flag==1){
	    $topology.=$_;
	}
    }
    close MODHMM or die;
    return($topology);
}
sub parse_PRODIV_results{
    my $topology="";
    my ($modhmmfile)=@_;

    my $flag=0;
    open(MODHMM,$modhmmfile) or die "can't open modhmmfile: $!\n$modhmmfile\n";
    while(<MODHMM>){
	chomp;
	
### I think this one is a little bit dangerous...
#	if(/^Labeling/ and $flag=0){
	if(/^Posterior/){
	    $flag=0;
	}
	if((/^o/ or /^i/ or /^O/) and $flag==0){
	    $flag=1;
	}
	if($flag==1){
	    $topology.=$_;
	}
    }
    close MODHMM or die;
    return($topology);
}

sub run_PSIBLAST{
    my ($fastafile)=@_;
    my $DBNAME="tmhmm_on_uniref50.fasta";
    my $DB="$zpreddir/DB/".$DBNAME;

#    my $tmpdir="$tmpdir/";
#`rsync --size-only  $DB $tmpdir`;
#`rsync --size-only  $DB.phr $tmpdir`;
#`rsync --size-only  $DB.pin $tmpdir`;
#`rsync --size-only  $DB.psq $tmpdir`;
    
    my $code2="-m 0 -h 1.e-5 -e 1 -F F -G 11 -E 1 -v 15000 -b 15000 -M BLOSUM62 ";
    
    my $BLAST2="$blastdir/bin/blastpgp ".$code2;
    
    #check whether there are 2 cpu's
    my $numcpu=`/bin/dmesg | /bin/grep CPU1`;
    if ($numcpu ne ''){
#	print "Using Dual CPU's\n";
	$BLAST2.=" -a 2";
    }
    
#    $ENV{'LD_LIBRARY_PATH'}="/afs/pdc.kth.se/i386_linux62/usr/local/vol/intel-compilers/compiler50/ia32/lib/";
    
    my @fields=split(/\//,$fastafile);
    
    my $dir="$tmpdir/";
    my $matrixname;
    $matrixname=$dir.$fields[scalar(@fields)-1];
    $matrixname.=".psi";
    `$BLAST2 -d $DB -i $fastafile -j 2 -Q $matrixname`;
}



sub NN_fwd
{
    
    my ($input,$netfile)=@_;
    my @input=split(/\s+/,$input);
    my ($nhidden,$counter,@temp,$nin);
    open(NET,"$netfile") or die "Cannot open $netfile\n$!\n";
    my @w1=();
    my @b1=();
    my @w2=();
    my $b2="";
    my $collect="";
    while(<NET>)
    {
	chomp;
	if(/^\D/ && not(/^-/))
	{
	    if(/^nhidden/)
	    {
		@temp=split(/\s+/);
		$nhidden=$temp[1];
	    }
	    elsif(/^nin/)
	    {
		@temp=split(/\s+/);
		$nin=$temp[1];
	    }
	    $collect=$_;
	    #print $collect."\n";
	    $counter=0;
	}
	else
	{
	    if($collect eq "w1")
	    {
		#print $_."\n";
		@temp=split(/\s+/);
		#print scalar @temp;
		#print "\n";
		for(my $i=0;$i<$nhidden;$i++)
		{
		    $w1[$counter][$i]=$temp[$i];
		    #print " $w1[$counter][$i]";
		}
		#print "\n";
		$counter++;
	    }
	    if($collect eq "w2")
	    {
		#print $_."\n";;
		push(@w2,$_);
	    }
	    if($collect eq "b1")
	    {
		@temp=split(/\s+/);
		@b1=@temp;
		#print @b1;
	    }
	    if($collect eq "b2")
	    {
		$b2=$_;
		#print $b2;
	    }

	}
    }
    my @layer=();
    my $sum=0;
    my $sum2=0;
    for(my $n=0;$n<$nhidden;$n++)
    {
	$sum=0;
	for(my $i=0;$i<$nin;$i++)
	{
	    $sum+=$w1[$i][$n]*$input[$i];
	    #print "$w1[$i][$n] $input[$i]\n";
	}
	$sum+=$b1[$n];
	#print "$sum"; 
	$sum=tanh($sum);
	#print "$sum\n";
	push(@layer,$sum);
    }
    #$sum2=0;
    for(my $n=0;$n<$nhidden;$n++)
    {
	$sum2+=$w2[$n]*$layer[$n];
	
    }
    $sum2+=$b2;
    return $sum2;
}
