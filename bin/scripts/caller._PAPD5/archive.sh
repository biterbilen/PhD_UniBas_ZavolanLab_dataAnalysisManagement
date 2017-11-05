#!/bin/bash

extensions4backup=(.pl .R .sh .m)
indir=_PAPD5
odir=_$indir

mkdir -p $odir
for ext in ${extensions4backup[@]}; do
	find $indir -name "*$ext"
done | grep -v -P "0\.\d+\.sh" > $odir/files4backup 

less $odir/files4backup | while read i; do
	echo copying $i
	dir=_`perl -e ' $f = shift; $f =~ /(.*)\/([^\/]+)$/; print $1;' $i`;
	if [ ! -e $dir ]; then mkdir -p $dir; fi
	cp -prf $i _$i;
done
