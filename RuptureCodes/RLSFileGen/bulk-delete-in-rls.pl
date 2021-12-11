#!/usr/bin/perl -w

$count=0;
$bulk_delete="";
while (<>) {
  chomp;
  @field=split(", ");
  $lfn=$field[0];
  $pfn=$field[1];

  $count++;
  
  print ("$lfn\n");
  $bulk_delete .= "$lfn $pfn "; 

  # Bulk delete every 500 entries, don't forget to delete the remaining files at the end
  if ($count % 500 == 0) {
    print ("Bulk delete at $count\n");
    system("globus-rls-cli bulk delete $bulk_delete rls://shock.usc.edu");

    # Clear out entries that have been deleted
    $bulk_delete="";
  }
}

# Delete the remaining entries
system("globus-rls-cli bulk delete $bulk_delete rls://shock.usc.edu");
print ("Finished deleting $count files.\n");
