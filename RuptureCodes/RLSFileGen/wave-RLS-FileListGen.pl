#!/usr/bin/perl -s
#Author: Hunter Francoeur
# A recursive script to scan directories and generate RLS in the format Month_Day_Year_RLS.txt. Directory to start at is the first command arg.
use Cwd;
use File::Find;
$gridURL = "gsiftp://wave.usc.edu";
my ($startdir) = &cwd; # keep track of the beginning
#Create the RLS file and print it to the correct rupture directory
find(\&printRLSFile, $ARGV[0]);
sub printRLSFile{
    if(/.rvm/ || /.txt.variation./){
	my ($sec,$min,$hour,$mday,$mon,$year,
          $wday,$yday,$isdst) = localtime time;
	unless(open(MAPFILE,">>" . $startdir . "/" . $mon . "_" . $mday . "_" . $year . "_RLS" . ".rls")){#open the file handler
		die "Could't create output file: " . @nameparts[0] . ".rls";
	}
	@dirParts = split ("\\.",$File::Find::dir);
	$currentdir = $dirParts[$#dirParts];
	print MAPFILE $_ . ", " . $gridURL . $currentdir . "/" . $_ . ", " . "wave\n";#write data
	close(MAPFILE);#close the file handler
	$flag=1;
    }
}

