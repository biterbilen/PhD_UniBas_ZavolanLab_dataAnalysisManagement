#!/bin/bash
SGE_TASK_ID=$1 #TODO doesn't work from clusters

outdir=Analysis
soutdir=RawData_mRNA
soutdir=ScaleNormalized_mRNA
soutdir=UniqueRawData_mRNA
soutdir=UniqueRawData_ALL
soutdir=UniqueRawData_ALL_woContamination

#TODO set
sid_prot_f=/import/bc2/home/zavolan/bilebi00/_DIS3L2/data/protein_sampleid_list
mkdir -p $outdir/$soutdir; pushd $outdir/$soutdir;

lib=`head -n $SGE_TASK_ID $sid_prot_f | tail -n 1 | cut -f 2`;

rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/DB_all_$lib.tar.gz .

tar -zxvf DB_all_$lib.tar.gz
rm -rf DB_all_$lib.tar.gz

echo $SGE_TASK_ID  $lib is DONE
