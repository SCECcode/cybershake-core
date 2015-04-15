#!/usr/bin/perl -w

$count=0;
$bulk_query="";

while (<>) {
  chomp;
#  @field=split(", ");
#  $lfn=$field[0];
#  $pfn=$field[1];
#  $pool=$field[2];
  $lfn=$_;

  $count++;
  
#  print ("$lfn\n");
  $bulk_query .= "$lfn "; 

  # Bulk add every 100 entries, don't forget to add the remaining files at the end
  if ($count % 100 == 0) {
    system("globus-rls-cli bulk query lrc lfn $bulk_query rls://wave.usc.edu");

    # Clear out entries that have been added
    $bulk_query="";

  }
#  system("sleep 1");
}

system("globus-rls-cli bulk query lrc lfn $bulk_query rls://wave.usc.edu");
