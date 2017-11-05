#!/bin/bash

#for i in ./data/*ces; do echo -e $i "\n" `less $i | perl -e '$i=36; $_=<>; while(<>){ @t=split; next if (length $t[1] < $i);$k=substr($t[1], $i, 4); for(1..$t[2]) { print $k, "\n";} }' | sort | uniq -c | sort -k 1,1gr` ; done  | less -S > max36nc.last4nc

for l in 32 35 36; do
	echo $l
	for i in ./data/*ces; do echo -e $i "\n" `less $i | perl -e '$i=shift; $_=<>; while(<>){ @t=split; next if (length $t[1] > $i);$k=substr($t[1], length($t[1])-4, 4); for(1..$t[2]) { print $k, "\n";} }' $l | sort | uniq -c | sort -k 1,1gr` ; done  | less -S > max35nc.last4nc
done
	exit;

for i in ./data/*ces; do echo -e $i "\n" `less $i | perl -e '$i=35; $_=<>; while(<>){ @t=split; next if (length $t[1] > $i);$k=substr($t[1], length($t[1])-4, 4); for(1..$t[2]) { print $k, "\n";} }' | sort | uniq -c | sort -k 1,1gr` ; done  | less -S > max35nc.last4nc

for i in ./data/*ces; do echo -e $i "\n" `less $i | perl -e '$i=32; $_=<>; while(<>){ @t=split; next if (length $t[1] > $i);$k=substr($t[1], length($t[1])-4, 4); for(1..$t[2]) { print $k, "\n";} }' | sort | uniq -c | sort -k 1,1gr` ; done  | less -S > max32nc.last4nc

#TODO 
#how much are there categories similar to each other; do we have a bias for long sequences where one kmer is always overrepresented?
