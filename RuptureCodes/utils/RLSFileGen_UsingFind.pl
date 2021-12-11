#!/usr/bin/perl -s
#Author: Hunter Francoeur
# A recursive script to scan directories and generate RLS in the format Month_Day_Year_RLS.txt. Directory to start at is the first command arg. 

use Cwd;
use File::Find;

if (@ARGV != 5) {
  print "Usage: RLSFileGen_UsingFind.pl <Rupture Variation Directory> <GridFTP URL> <Pool> <ERF_ID> <Rup_Var_ID>\n";
  print "  For example:\n    RLSFileGen_UsingFind.pl /data/CyberShake/RuptureVariations gsiftp://tg-gridftp.sdsc.teragrid.org/ sdsc 35 4\n";
  exit;
}

$rupture_variation_directory = $ARGV[0];
$gridURL = $ARGV[1];
print $gridURL;
$pool = $ARGV[2];
$erf_id = $ARGV[3];
$rup_var_id = $ARGV[4];
#$gridURL = "gsiftp://tg-gridftp.sdsc.teragrid.org/gpfs/scec";
#$pool = "sdsc";

my ($startdir) = &cwd; # keep track of the beginning

# Create the RLS file and print it to the correct rupture directory
find(\&printRLSFile, $rupture_variation_directory);

sub printRLSFile{
    #if(/.rvm/ || /.txt.variation-/){
    if(/.rvm/ || /_event/){
	my ($sec,$min,$hour,$mday,$mon,$year,
          $wday,$yday,$isdst) = localtime time;
	unless(open(MAPFILE,">>" . $startdir . "/" . $pool . "-" .$mon . "_" . $mday . "_" . $year . ".rls")){#open the file handler
		die "Could't create output file: " . @nameparts[0] . ".rls";
	}
	@dirParts = split ("\\.",$File::Find::dir);
	$currentdir = $dirParts[$#dirParts];
	print MAPFILE "e". $erf_id . "_rv" . $rup_var_id . "_" . $_ . " " . $gridURL . $currentdir . "/" . $_ . " pool=" . $pool . "\n"; #write data
	close(MAPFILE); #close the file handler
	$flag=1;
    }
}
