#!/usr/bin/perl -w

$count=0;
$bulk_query="";

while (<>) {
  chomp;
#  @field=split(", ");
#  $lfn=$field[0];
#  $pfn=$field[1];
#  $pool=$field[2];
  $pfn=$_;

  $count++;
  
  $bulk_query .= "$pfn "; 
#  $bulk_attribute_add .= "$pfn pool pfn string $pool ";

  # Bulk add every 100 entries, don't forget to add the remaining files at the end
  if ($count % 100 == 0) {
    system("globus-rls-cli attribute bulk query pool pfn $bulk_query rls://wave.usc.edu");

    # Clear out entries that have been added
    $bulk_query="";

  }
#  system("sleep 1");
}

# Register the remaining files
system("globus-rls-cli attribute bulk query pool pfn $bulk_query rls://wave.usc.edu");
