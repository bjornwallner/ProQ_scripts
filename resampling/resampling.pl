#!/usr/bin/perl -w
use Cwd 'abs_path';
use File::Basename;
my $SCRIPT_DIR=dirname(abs_path($0));

$rosetta_path="/home/x_bjowa/nobackup/git/rosetta_3.5/rosetta_source/bin/";
$rosetta_db="/home/x_bjowa/nobackup/git/rosetta_3.5/rosetta_database/";

$basename=$ARGV[0];
$repack_cmd="$rosetta_path/relax.linuxgccrelease  -database $rosetta_db -in:file:fullatom -out:file:silent_struct_type binary -relax:membrane -membrane:Membed_init -score:weights membrane_highres.wts -in:file:spanfile $basename.span -nstruct 10 -relax_script $SCRIPT_DIR/resampling/repack.script -in:file:s *.pdb -out:file:silent ProQM.repacked.silent";
$score_cmd="$rosetta_path/score.linuxgccrelease -database $rosetta_db -in:file:fullatom -ProQ:membrane -score:weights ProQM -in:file:silent ProQM.repacked.silent -out:file:scorefile ProQM.repacked.silent.score";
print $repack_cmd."\n";
print $score_cmd."\n";
system("$repack_cmd");
system("$score_cmd");
