#!/usr/bin/perl
use Cwd 'abs_path';
use File::Basename;
my $SCRIPT_DIR=dirname(abs_path($0));

#sspro
$sspro="$SCRIPT_DIR/apps/sspro4";
print $sspro."\n";
#exit;
chdir("$sspro");
symlink("../../apps/blast-2.2.18_x86_64/bin/","$sspro/blast2.2.8");
system("./configure.pl");


