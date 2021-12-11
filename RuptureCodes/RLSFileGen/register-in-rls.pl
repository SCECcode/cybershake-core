#!/usr/bin/perl -w

while (<>) {
  chomp;
  @field=split(", ");
  $lfn=$field[0];
  $pfn=$field[1];
  $pool=$field[2];

  print ("$lfn - $pool\n");
#  system("globus-rls-cli create $lfn $pfn rls://wave.usc.edu\n");
  system("globus-rls-cli add $lfn $pfn rls://wave.usc.edu");
  system("globus-rls-cli attribute add $pfn pool pfn string $pool rls://wave.usc.edu");

#  system("sleep 1");

}
