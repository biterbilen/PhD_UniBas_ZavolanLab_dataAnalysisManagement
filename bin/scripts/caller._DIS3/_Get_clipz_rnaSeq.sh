#!/bin/bash

outdir=data/clipz_RNAseq

mkdir -p $outdir; pushd $outdir

ids=(211 212 213   223 225 224   199 200   230 232 231)

ps=(mRNAseq-HEK-0 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-siDis3L2   totRNAseq-HEK-0 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2   totRNAseq-Hela-0 totRNAseq-Hela-siDis3L2   totRNAseq-HEK-0_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2)

for SGE_TASK_ID in `seq 1 ${#ids[@]}`; do
	id=${ids[$((SGE_TASK_ID-1))]}
	p=${ps[$((SGE_TASK_ID-1))]}
	scp mirz@bc2-web07.bc2.unibas.ch:~mirz/clipzServer/data/samples/$id/transcript_expression $p.exp;
done	
