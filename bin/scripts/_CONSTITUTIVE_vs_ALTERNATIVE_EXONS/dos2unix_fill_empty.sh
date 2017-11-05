#!/bin/bash

file=$1;

r=$RANDOM;

less $file | sed 's//\n/g' > $r.$file
#less $r.$file | perl -e 'while(<>) { @t=split/\t/; die "Dieing $#t:@t\n" if ($#t!=15); for my $i(0..$#t) { $t[$i]="NA" if ($t[$i] eq ""); $t[$i]="NA\n" if ($t[$i] eq "\n"); print "$i^$t[$i]^\n"; } print join("\t",@t); }' 
less $r.$file | perl -e 'while(<>) { @t=split/\t/; if ($#t!=15) {print STDERR "Died $#t:@t"; next;}; for my $i(0..$#t) { $t[$i]="NA" if ($t[$i] eq ""); $t[$i]="NA\n" if ($t[$i] eq "\n"); } print join("\t",@t); }' 

rm $r.*;

