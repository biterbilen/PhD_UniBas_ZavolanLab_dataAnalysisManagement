#!/bin/bash

outdir=Analysis

pushd $outdir

soutdir=Stupid_stats
mkdir $soutdir; pushd $soutdir;

#expressed

#clipped by different regions
cut=15
less ../EWSR1.annot | awk '$1 == "mRNA" && $NF>c && $10==2 && $(NF-1)>c { print $3;}' c=$cut | sort | uniq > samewindow12.expressed

less ../EWSR1.annot | awk '$1 == "mRNA" && $NF==0 && $10==1 && $(NF-1)>c { print $3;}' c=$cut | sort | uniq > 1.expressed
less ../EWSR1.annot | awk '$1 == "mRNA" && $(NF-1)==0 && $10==1 && $NF>c { print $3;}' c=$cut | sort | uniq > 2.expressed

#common gene names clipped by different regions
~bilebi00/_PAPD5/scripts/innerJoin.pl 1.expressed  2.expressed  1 1 1 > 12.expressed 





