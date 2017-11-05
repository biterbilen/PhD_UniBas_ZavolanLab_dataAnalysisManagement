#!/bin/bash

outdir=data/clipz_RNAseq
pwd=`pwd`
mkdir -p $outdir; pushd $outdir

ids=(376 377 634 635)

ps=(mRNAseq_Sam68_KO_cerebellum mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_KO_brainStem mRNAseq_Sam68_WT_brainStem)

for SGE_TASK_ID in `seq 1 ${#ids[@]}`; do
	id=${ids[$((SGE_TASK_ID-1))]}
	p=${ps[$((SGE_TASK_ID-1))]}
	#outfile is the second argumenta and supplied with the full path
	echo ~bilebi00/../rodak/progs/clipz/analysis/siteExtraction/generateTranscriptExpression.pl $id $pwd/$outdir/$p.exp 
	nohup ~bilebi00/../rodak/progs/clipz/analysis/siteExtraction/generateTranscriptExpression.pl $id $pwd/$outdir/$p.exp &> LOG.$id & #don't run it in the clusters; it is using 
	#scp mirz@bc2-web07.bc2.unibas.ch:~mirz/clipzServer/data/samples/$id/transcript_expression $p.exp;
done	
