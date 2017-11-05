#!/bin/bash
#$ -S /bin/bash
# $ -q fs_short@@qc_nehalem
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.Get_clip_rnaseq_intersect_byGeneID_$TASK_ID
#$ -t 1-1

#TODO writeme 
handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> LOG.Get_clip_rnaseq_intersect_byGeneID_$SGE_TASK_ID
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PERL5LIB=$HOME/lib:$PERL5LIB


project=_DIS3;
if [ $project == "_DIS3" ]; then
	moutdir=Stepanka
	sa1=(totRNAseq.HEK.0 totRNAseq.HEK.0 mRNAseq.HEK.0 mRNAseq.HEK.0 totRNAseq.Hela.0)
	sa2=(totRNAseq.HEK.oeDis3L2 totRNAseq.HEK.siDis3L2 mRNAseq.HEK.oeDis3L2 mRNAseq.HEK.siDis3L2 totRNAseq.Hela.siDis3L2)
elif [ $project == "_EWSR1" ]; then
	moutdir=Shivendra
	sa1=(mRNAseq_siCTRL-A_HeLa)
	sa2=(mRNAseq_siEWSR1_HeLa)
elif [ $project == "_SAM68" ]; then
	moutdir=Scheiffele
	sa1=(mRNAseq_Sam68_KO_brainStem mRNAseq_Sam68_KO_cerebellum mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_KO_cerebellum)
	sa2=(mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_KO_brainStem)
elif [ $project == "_KNOCKDOWNS" ]; then
	moutdir=_Zavolan
	sa1=(mRNASeq_siRNA_GFP mRNASeq_No_4SU_No_XL_rep_A mRNASeq_No_4SU_No_XL_rep_B)
	sa2=(mRNASeq_siRNA_HuR mRNASeq_siRNA_HuR mRNASeq_siRNA_HuR)
else
	echo Unknown project $project;
	exit;
fi

#TODO set
gi=gene; #gene; #isoform; #promoters; #splicing; #tss_group

outdir=$moutdir/CLIP_comparison_cuffdiff${trusted}_${src}_$gi

mkdir -p $outdir; pushd $outdir

s1=${sa1[$(($SGE_TASK_ID - 1))]}
s2=${sa2[$(($SGE_TASK_ID - 1))]}

p1=${s1}VS${s2}_$gi


~/_PAPD5/scripts/outerJoin.pl CLIP_comparison_clipz_RNAseq_qn_superClusters_gene/totRNAseq.HEK.0VStotRNAseq.HEK.oeDis3L2_gene.diff CLIP_comparison_cuffdiff_trusted_superClusters_gene/totRNAseq-HEK-oeDis3L2_2VStotRNAseq-HEK-oeDis3L2_gene.diff 3 3 '1-3,8-11,23-26' > clipz_RNAseq_qn-cuffdiff_trusted.totRNAseq-HEK-oeDis3L2_2VStotRNAseq-HEK-oeDis3L2_gene.diff_merged

~bilebi00/bin/R --no-save --args $i < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/.R 



