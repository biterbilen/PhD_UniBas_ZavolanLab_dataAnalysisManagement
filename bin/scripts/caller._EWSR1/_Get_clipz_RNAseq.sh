#!/bin/bash

wd=_EWSR1
outdir=data/clipz_RNAseq

mkdir -p $outdir; pushd $outdir

ids=(`less ~bilebi00/$wd/data/protein_sampleid_list | grep mRNAseq | cut -f 2`)

ps=(`less ~bilebi00/$wd/data/protein_sampleid_list | grep mRNAseq | cut -f 3`)
#ps=(mRNAseq_siCTRL-A_HeLa mRNAseq_siEWSR1_HeLa)

for SGE_TASK_ID in `seq 1 ${#ids[@]}`; do
	p=${ps[$((SGE_TASK_ID-1))]}
	if [ -e $p.exp ]; then
		continue;
	fi
	id=${ids[$((SGE_TASK_ID-1))]}
	rsync mirz@bc2-web08.bc2.unibas.ch:~mirz/clipzServer/data/samples/$id/transcript_expression $p.exp;
	echo $id
done	
