#!/usr/bin/perl -w
use Cwd 'abs_path';
use File::Basename;

my $LWP_loaded=1;
my $HTML_loaded=0;

eval {
    require LWP::UserAgent; 
    LWP::UserAgent->import();
    1;
} or do {
   my $error = $@;
   $LWP_loaded=0

};

eval {
    require HTML::Parser;
    HTML::Parser->import();
    1;
} or do {
   my $error = $@;
   $HTML_loaded=0;
};


my $DB="/proj/wallner/users/x_bjowa/DB/uniref90.fasta";



if(scalar(@ARGV) % 2 != 0) {
    print STDERR "Usage: run_all_external.pl -pdb [pdbfile] -fasta [fastafile] -membrane 1 (if membrane protein) -overwrite 1 (if overwrite)\n";
    exit;
}
%param=@ARGV;
if(not(defined($param{-pdb})) &&
   not(defined($param{-fasta}))) {
    print STDERR "Usage: run_all_external.pl -pdb [pdbfile] -fasta [fastafile] -membrane 1 (if membrane protein) -overwrite 1 (if overwrite)\n";
    exit;
}
my $membrane=0;
my $overwrite=0;
if(defined($param{-overwrite})) {
    $overwrite=1;
}
if(defined($param{-membrane})) {
    $membrane=1;
}

if(defined($param{-pdb}))
{
    $infile=$param{-pdb};
} else {
    $infile=$param{-fasta};
}

if(defined($param{-db}))
{
    $DB=$param{-db};
}

my $cpu=8;
if(defined($param{-cpu}))
{
    $cpu=$param{-cpu};
}
#my $use_master=0;
#if(defined($param{-subsetof}))
#{
#    $master_fasta=$param{-subsetof};
#    $use_master=1;
#}

my $INSTALL_DIR=dirname(abs_path($0));
$INSTALL_DIR.="/.."; 
#my $DB=$INSTALL_DIR."/DB/uniref90.fasta";
my $full_path=abs_path($infile);
my $path=dirname($full_path);
my $pdb=basename($full_path); #the name is $pdb but it could be a fasta file...
chdir($path);
print "run all for $pdb in $path...\n";

#`mv -f $pdb $pdb.orig`;
#`grep ^ATOM $pdb.orig | $INSTALL_DIR/bin/kill_chain.pl  > $pdb`;
 
#fasta
$fasta="$pdb.fasta";
if(defined($param{-pdb})) {
    $seq=`$INSTALL_DIR/bin/aa321CA.pl $pdb`;
} else {
    $seq=`grep -v '>' $pdb`;
    $seq=~s/\n//g;
    $seq=~s/\s+//g;
}

chomp($seq);
if($overwrite || !-e $fasta) {
    print "Creating fasta file..\n";
    
    open(OUT,">$fasta");
    print OUT ">$pdb\n$seq\n";
    close(OUT);
}
#exit;
#print $seq."\n";
$profile_file="$pdb.psi";
#
if($overwrite) {
	unlink("$pdb.psi");
	unlink("$pdb.mtx");
	unlink("$pdb.chk");
	unlink("$pdb.ss2");
}
if(!-e $profile_file) {
    print "Creating profiles and predicting ss\n";
    
   # print "/Users/bjornw/Research/bin/create_profile.sh $fasta\n";
   # `/Users/bjornw/Research/bin/create_profile.sh $fasta`;
  #  if($use_master) {
#	print "using $master_fasta";
#	`$INSTALL_DIR/bin/profile_subset.pl 
#	    
 #   } else {

    print "$INSTALL_DIR/bin/create_profile.sh $fasta $DB $cpu\n";
#    exit;
    system("$INSTALL_DIR/bin/create_profile.sh $fasta $DB $cpu");
  #  }
}

if($membrane)
{
    if($LWP_loaded && $HTML_loaded) {
	$topcons="$pdb.topcons";
	$topcons_fa="$pdb.topcons.fa";
	if($overwrite || !-e $topcons) {
	    print "topcons...\n";
	    my $ua = new LWP::UserAgent;
	    my $response = $ua -> post('http://topcons.net',{'sequence' => $seq,'do' => 'Submit',});
	    $content=$response->content();
	    if($content=~/result\/(.+)\/topcons.txt/) {
		$id=$1;
		print $id."\n";
		my $response = $ua -> post("http://topcons.net/result/$id/topcons.txt");
		open(OUT,">$topcons");
		print OUT $response->content();
		close(OUT);
		#my $parser = new MyParser;
#	my $parsed=$parser->parse($content);
#	print $parsed."\n";
#	$string="sequence=$seq\&do=Submit";
#	$output=`echo "$string" | lynx -post_data http://topcons.net`;
#	    exit;
#	open(OUT,">$topcons");
#	print OUT $output;
#	close(OUT);
	    }
#`run_topcons.pl $fasta > $topcons`;
	}
	if($overwrite || !-e "$topcons.fa"){ 
	    print "Topcons to fasta...\n";
#	`$INSTALL_DIR/bin/parse_topcons.pl $topcons > $topcons.fa`;
	    my ($seq2,$topo2)=parse_topcons($topcons);
	    open(TOPCONS,">$topcons.fa");
	    print TOPCONS ">$topcons\n$topo2\n";
	    close(TOPCONS);
	}
	$spanfile="$pdb.span";
	if($overwrite || !-e $spanfile) {
	    my ($seq2,$topo2)=parse_topcons($topcons);
	    my $span=Mspan($topo2);
	    open(SPAN,">$spanfile");
	    print SPAN "Prediction from TOPCONS on $ARGV[0]\n";
	    print SPAN $span;
	    close(SPAN);
	    
	    
	}

	my $zpred = "$INSTALL_DIR/apps/zpred/zpred.pl"; #/afs/pdc.kth.se/home/k/kriil/vol_03/Programs/zpred/bin/zpred.pl";
	$zpred_file="$pdb.zpred";
	$temp_dir="/tmp/";
	if($overwrite || !-e $zpred_file) {
	    print "Zpred...\n";
	    print "$zpred -mode profile_topology -profile $profile_file  -topology $topcons_fa -prediction modhmm -out $zpred_file -tmpdir $temp_dir\n";
	    `$zpred -mode profile_topology -profile $profile_file  -topology $topcons_fa -prediction modhmm -out $zpred_file -tmpdir $temp_dir`;
#exit;
	}
	
	$mprap_file="$pdb.mpSA";
	if(-e $mprap_file) {
	    if(-s $mprap_file==0) {
		`rm $mprap_file`;
	    }
	    
	}
	if($overwrite || !-e $mprap_file) {
#	@temp=split(/\//,$fasta);
#	$outdir=join("/",@temp[0..$#temp-1]);
#print $outdir."\n";
	    # print "scripts/run_MPSA.py $fasta $outdir/\n";
#exit;
#	$outdir="." if(length($outdir)==0);
	    print "MPRAP\n";
#    print "/home/bjornw/afs/.vol/bjornw27/MPSA/run_MPSA.py $fasta $outdir\n";
#	print "$INSTALL_DIR/bin/MPSA/run_MPSA.py $fasta $outdir\n";
	    print "$INSTALL_DIR/apps/MPSA/run_MPSA.py $fasta $path $DB\n";
	    `$INSTALL_DIR/apps/MPSA/run_MPSA.py $fasta $path $DB | egrep ' E|B ' > $mprap_file`;
	    
	    
	}
    } else {
	print STDERR "Unable to load LWP::UserAgent or HTML::Parser, cannot run topcons server!\n";
    }
} else {
    $accfile="$pdb.acc";
    if($overwrite || !-e $accfile) {
	#pen(TEMP, ">$pdb.seq_file.tmp") || die "can't create temporary file.\n";
	#rint TEMP "1 20 3\n"; #create a title line required by ACCpro
	#rint TEMP "NONAME\n"; 
	#rint TEMP "$seq"; 
	#lose(TEMP);
	#`$INSTALL_DIR/apps/sspro4/script/process-blast.pl $pdb.fasta.blastpgp $pdb.msa_for_acc $pdb.fasta`;
	#print "$INSTALL_DIR/apps/sspro4/script/predict_seq_sa $INSTALL_DIR/apps/sspro4/model/accpro.model
	#
	#print "$INSTALL_DIR/apps/sspro4/script/homology_sa.pl $INSTALL_DIR/apps/sspro4/ $pdb.fasta $pdb.msa_for_acc $pdb.accpro\n";
	`$INSTALL_DIR/apps/sspro4/bin/predict_acc.sh $fasta $accfile.out`;
	`cat $accfile.out |tail -n 1 > $accfile`;

    }
}

sub parse_topcons
{
    my $file=shift;
    my $topcons="";
    my $seq="";
    my $octopus="";
    my $get_prediction=0;
    my %pred=();
    $pred{'TOPCONS'}="";
    open(FILE,$file);
    while(<FILE>)
    {
	chomp;
	
	if($get_sequence) {
	    if(length($_)==0) {
		$get_sequence=0;
	    } else {
		$seq.=$_;
	    }
	}
	if($get_prediction) {
	    if(length($_)==0) {
		$get_prediction=0;
	    } else {
		$pred{$method}.=$_;
	    }
	}
	if(/(.+)\spredicted topology:/) {
	    $method=$1;
	    $get_prediction=1;
	}
	$get_sequence=1 if(/Sequence:/);


    }
    close(FILE);
    $seq=~s/\s+//g;
 #   print $seq."\n";
    foreach $method(keys(%pred)) {
	$pred{$method}=~s/\s+//g;
#	print "$method\n$pred{$method}\n";
    }
    $pred{'TOPCONS'}=$pred{'OCTOPUS'} if(length($pred{'TOPCONS'})==0);
    if(length($topcons)==0) {
	$topcons=$seq;
	$topcons=~s/./i/g;
    }
   # "FILE $file\n";
    #print "SEQ     $seq\n";
    #print "TOPCONS $topcons\n";

    return($seq,$pred{'TOPCONS'});
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
