#!/usr/bin/perl -w

$count=0;
$bulk_add="";
$bulk_attribute_add="";
while (<>) {
  chomp;
  @field=split(", ");
  $lfn=$field[0];
  $pfn=$field[1];
  $pool=$field[2];

  $count++;
  
  print ("$lfn - $pool\n");
  $bulk_add .= "$lfn $pfn "; 
  $bulk_attribute_add .= "$pfn pool pfn string $pool ";

  # Bulk add every 500 entries, don't forget to add the remaining files at the end
  if ($count % 500 == 0) {
    print ("Bulk adding at $count\n");
#    system("globus-rls-cli bulk create $bulk_add rls://wave.usc.edu\n");
    system("globus-rls-cli bulk add $bulk_add rls://shock.usc.edu");
    system("globus-rls-cli attribute bulk add $bulk_attribute_add rls://shock.usc.edu");

    # Clear out entries that have been added
    $bulk_add="";
    $bulk_attribute_add="";
  }
#  system("sleep 1");
}

# Register the remaining files
#system("globus-rls-cli bulk create $bulk_add rls://wave.usc.edu\n");
system("globus-rls-cli bulk add $bulk_add rls://shock.usc.edu");
system("globus-rls-cli attribute bulk add $bulk_attribute_add rls://shock.usc.edu");
print ("Finished adding $count files.\n");

