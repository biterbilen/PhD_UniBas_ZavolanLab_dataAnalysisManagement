#!/bin/bash
#$ -S /bin/bash
# $ -q fs_short@@qc_nehalem
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.Get_clip_rnaseq_intersect_byGeneID_w_cuffdiff_RNAseq_$TASK_ID
#$ -t 1-11
# $ -t 1-3
# $ -t 1-4
# $ -t 1-8
# $ -t 1-1

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> LOG.Get_clip_rnaseq_intersect_byGeneID_w_cuffdiff_RNAseq_$SGE_TASK_ID
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PERL5LIB=$HOME/lib:$PERL5LIB

project=_SAM68; #_DIS3; #_EWSR1; #_DIS3; #_EWSR1
if [ $project == "_DIS3" ]; then
	moutdir=Stepanka
	sa1=(mRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-0 mRNAseq-HEK-0       totRNAseq-HEK-0_2 totRNAseq-HEK-0_2 totRNAseq-HEK-0 totRNAseq-HEK-0 mRNAseq-HEK-0 mRNAseq-HEK-0 totRNAseq-Hela-0)
	sa2=(totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-0 totRNAseq-HEK-0_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-siDis3L2 totRNAseq-Hela-siDis3L2)
elif [ $project == "_EWSR1" ]; then
	moutdir=Shivendra
	sa1=(mRNAseq_siCTRL-A_HeLa)
	sa2=(mRNAseq_siEWSR1_HeLa)
elif [ $project == "_SAM68" ]; then
	moutdir=__Scheiffele; #TODO change
	sa1=(mRNAseq_Sam68_KO_brainStem mRNAseq_Sam68_KO_cerebellum mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_KO_brainStem)
	sa2=(mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_KO_cerebellum)
elif [ $project == "_KNOCKDOWNS" ]; then
	moutdir=_Zavolan
	sa1=(mRNASeq_siRNA_GFP mRNASeq_No_4SU_No_XL_rep_A mRNASeq_No_4SU_No_XL_rep_B)
	sa2=(mRNASeq_siRNA_HuR mRNASeq_siRNA_HuR mRNASeq_siRNA_HuR)
else
	echo Unknown project $project;
	exit;
fi

#TODO set it from the density plot
minfpkm=0; #0.01 #0.1 #(e^(-4))
ps=10

norm=""; #.qn; #.Hela; #.qn;
gi=gene;
src="ScaleNormMultNoiseDist";
trusted=Cuffdiff_trusted
#clipfolder=~bilebi00/$project/data/SP/*\.data;
#clipfolder=~bilebi00/$project/Analysis/DIS3*\.data
#clipfolder=~bilebi00/$project/4WEB_DIS3/PC*\.data

normt=`echo $norm | cut -b 2-`;
if [ "`echo $clipfolder | grep \"data/SP\"`" != "" ]; then
	src=superClusters
elif [ "`echo $clipfolder | grep Analysis`" != "" ]; then
	src=REP
fi

#--------------------
cuffdiff_gexpd_file=../${trusted}_RNAseq/${gi}_exp.diff
cuffdiff_gfpkm_file=../${trusted}_RNAseq/${gi}s.fpkm_tracking

if [ "$gi" == "promoters" ] || [ "$gi" == "splicing" ]; then
	cuffdiff_gexpd_file=../${trusted}_RNAseq/${gi}.diff
	cuffdiff_gfpkm_file="";
	minfpkm=-1
fi

outdir=$moutdir/Significant_${trusted}_${gi}
#outdir=$moutdir/CLIP_comparison_${trusted}_${normt}_${src}_${gi}
mkdir -p $outdir; pushd $outdir

s1=${sa1[$(($SGE_TASK_ID - 1))]}
s2=${sa2[$(($SGE_TASK_ID - 1))]}

p1=${s1}VS${s2}_$gi

less $cuffdiff_gexpd_file | perl -e '$s1=shift; $s2=shift; $minfpkm=shift; while(<>) { if ($_=~/\s$s1\s$s2\sOK\s/) { @t=split; print $_ if ($t[7]>$minfpkm && $t[8]>$minfpkm); }}' $s1 $s2 $minfpkm > $p1.diff.tmp

if [ "$gi" == "promoters" ] || [ "$gi" == "splicing" ]; then
	awk 'BEGIN{OFS="\t";} { print $0, $1;}' $p1.diff.tmp > $p1.diff
else
	#select only primary transcripts #TODO think; could it be a problem for EWSR1?
	#~bilebi00/_PAPD5/scripts/innerJoin.pl $p1.diff.tmp $cuffdiff_gfpkm_file 1 1 '' | awk '$20 ~ /TSS/ { print $0;}' | cut -f 1-14,20 > $p1.diff
	~bilebi00/_PAPD5/scripts/innerJoin.pl $p1.diff.tmp $cuffdiff_gfpkm_file 1 1 '' | cut -f 1-14,20 > $p1.diff
	~bilebi00/bin/R --no-save --args $p1.diff < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/cuffdiff_sig_scatter.R > /dev/null
	i=$p1.diff
	less $i | awk '$8>$9 && $14~/yes/{ print; }' > $i.pos; 
	less $i | awk '$8<$9 && $14~/yes/{ print; }' > $i.neg;  
	echo $i `less $i | awk '$14 ~ /yes/ { print; }' | wc -l` `less $i.pos | wc -l` `less $i.neg | wc -l`; 
fi

exit; #TODO Do clip intersects in a different folder; TODO rename the script
if [ $project == "_SAM68" ]; then
	exit; #TODO remove this if you compare with Splice array results
fi	

header="`head $cuffdiff_gexpd_file -n 1` joinedTSS chr source name beg end score score2 str frm foreground foreground2 background background2 annot_type annot_ids gene_ids";

region=EXON #TODO set
if [ $src == superClusters ]; then
	region='.'
fi
rm -rf $p1.aggregate
wh=1;
for f2 in $clipfolder; do
	p2=`basename $f2 .data`;
	out=${p1}_$p2
	echo $out
	echo $header | sed 's/ /\t/g' > $out 
	rand=$RANDOM;
	if [ $p2 == PC_DIS3L2_1vPC_DIS3L2_2andPC_DIS3L2_1vPC_DIS3L2_2 ]; then
		ln -sf $f2 $p2.$rand
	elif [ $src == REP ]; then
		minT2C=1
		less $f2 | awk '$6>cut { print;}' cut=$minT2C > $p2.$rand
	elif [ $src == ScaleNormMultNoiseDist ]; then
		minzvalue=1
		less $f2 | awk '$6>cut && $7>cut{print;}' cut=$minzvalue > $p2.$rand
	else
		ln -sf $f2 $p2.$rand
	fi
	~bilebi00/_PAPD5/scripts/leftJoin.pl $p1.diff $p2.$rand 3 16 '' 0 >> $out #genename join
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/verify_clip_zscoreBins_with_foldChange.pl $out $minfpkm $region $ps $p1 $p2 $wh >> $p1.aggregate 
	wh=0;
	rm -rf $out;
done
~bilebi00/bin/R --no-save --args $p1.aggregate < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/zscoreBins_v_foldChange.R > /dev/null

#rm -rf $p1.diff;
rm $p1.diff.tmp

echo DONE




