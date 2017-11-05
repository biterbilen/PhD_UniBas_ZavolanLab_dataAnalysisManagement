#!/bin/bash
#$ -S /bin/bash
# $ -q fs_short@@qc_nehalem
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.Get_clip_rnaseq_intersect_byGeneID_w_clipz_rnaSeq_$TASK_ID
#$ -t 1-11
# $ -t 1-1

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> LOG.Get_clip_rnaseq_intersect_byGeneID_w_clipz_rnaSeq_$SGE_TASK_ID
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PERL5LIB=$HOME/lib:$PERL5LIB

project=_DIS3; #_DIS3; #_EWSR1
if [ $project == "_DIS3" ]; then
	moutdir=Stepanka
	sa1=(mRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-0 mRNAseq-HEK-0       totRNAseq-HEK-0_2 totRNAseq-HEK-0_2 totRNAseq-HEK-0 totRNAseq-HEK-0 mRNAseq-HEK-0 mRNAseq-HEK-0 totRNAseq-Hela-0)
	sa2=(totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-0 totRNAseq-HEK-0_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-siDis3L2 totRNAseq-Hela-siDis3L2)
elif [ $project == "_EWSR1" ]; then
	moutdir=Shivendra
	sa1=(mRNAseq_siCTRL-A_HeLa)
	sa2=(mRNAseq_siEWSR1_HeLa)
else
	echo Unknown project $project;
	exit;
fi

#TODO set it from the density plot
minfpkm=0; #0.01 #0.1 #(e^(-4))
ps=1

norm=".qn"; #.qn; #.Hela; #.qn; #"";
gi=gene;
src="ScaleNormMultNoiseDist";
trusted=clipz_RNAseq
clipfolder=~bilebi00/$project/data/SP/*\.data;
#clipfolder=~bilebi00/$project/Analysis/DIS3*\.data
#clipfolder=~bilebi00/$project/Analysis/EXC/DIS3*\.data
#clipfolder=~bilebi00/$project/4WEB_DIS3/PC*\.data
rnaSeqFile=~bilebi00/$project/data/clipz_RNAseq/RNAseq.raw$norm.w_geneSymbol
gif=~bilebi00/_EWSR1/data/hg19_TR.info.all

normt=`echo $norm | cut -b 2-`;
if [ "$norm" == "" ]; then normt=raw; fi
if [ "`echo $clipfolder | grep \"data/SP\"`" != "" ]; then
	src=superClusters
elif [ "`echo $clipfolder | grep EXC`" != "" ]; then
	src=REP_exclusive	
elif [ "`echo $clipfolder | grep Analysis`" != "" ]; then
	src=REP	
fi

echo $clipfolder
outdir=$moutdir/CLIP_comparison_${trusted}_${normt}_${src}_${gi}
mkdir -p $outdir; pushd $outdir

s1=${sa1[$(($SGE_TASK_ID - 1))]}
s2=${sa2[$(($SGE_TASK_ID - 1))]}

p1=${s1}VS${s2}_$gi

i=`less $rnaSeqFile | head -n 1 | sed  's/\t/\n/g' | grep -w id -n | awk -F ":" '{ print $1; }'`;
gi=`less $rnaSeqFile | head -n 1 | sed  's/\t/\n/g' | grep -w gid -n | awk -F ":" '{ print $1; }'`;
i1=`less $rnaSeqFile | head -n 1 | sed  's/\t/\n/g' | grep -w $s1 -n | awk -F ":" '{ print $1; }'`;
i2=`less $rnaSeqFile | head -n 1 | sed  's/\t/\n/g' | grep -w $s2 -n | awk -F ":" '{ print $1; }'`;

if [ "$i1" == "" ]; then exit; fi

#for the file in cuffdiff format putting . in unknown regions
#logFC and testStat fields are in the same order as in cuffdiff
#teststat is the same as logFC 
less $rnaSeqFile | cut -f $i,$i1,$i2,$gi | perl -e '$s1=shift; $s2=shift; $minfpkm=shift; while(<>) { chomp; @t=split; if ( $t[2] ne "NA" && $t[1]>$minfpkm && $t[2]>$minfpkm) { print join("\t", $t[0], ".", $t[3], ".", $s1, $s2, "OK", $t[1], $t[2], log($t[2])-log($t[1]), log($t[1])-log($t[2]), ".", ".", ".", "."), "\n"; } }' $s1 $s2 $minfpkm > $p1.diff 

#~bilebi00/_PAPD5/scripts/innerJoin.pl $rnaSeqFolder/$s1.exp $rnaSeqFolder/$s2.exp 1 1 '1,4,8' | perl -e '$minfpkm=shift; while(<>) { @t=split; print $_ if ($t[1]>$minfpkm && $t[2]>$minfpkm); }' $minfpkm > $p1.diff.tmp
#~bilebi00/_DIS3/scripts/assign_gene_id.pl $p1.diff.tmp $gif 0 NA > $p1.diff

header="test_id gene_id gene locus sample_1 sample_2 status value_1 value_2 ln(fold_change) test_stat p_value q_value significant joinedTSS chr source name beg end score score2 str frm foreground foreground2 background background2 annot_type annot_ids gene_ids";
#header="`head $cuffdiff_gexpd_file -n 1` joinedTSS chr source name beg end score score2 str frm foreground foreground2 background background2 annot_type annot_ids gene_ids";

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
	elif [ "`echo $src | grep REP`" != "" ]; then
		minT2C=1
		if [ "`echo $p2 | grep DIS3L2`" != "" ]; then
			minT2C=4;
		fi
		less $f2 | awk '$6 > cut && $14 == annot && $16 != "NA" { print;}' cut=$minT2C annot=$region | sort -k 6,6gr > $p2.$rand
	elif [ $src == ScaleNormMultNoiseDist ]; then
		minzvalue=1
		less $f2 | awk '$6>cut && $7>cut{print;}' cut=$minzvalue > $p2.$rand
	else
		ln -sf $f2 $p2.$rand
	fi
	~bilebi00/_PAPD5/scripts/leftJoin.pl $p1.diff $p2.$rand 3 16 '' 0 1 >> $out #genename join
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/verify_clip_zscoreBins_with_foldChange.pl $out $minfpkm $region $ps $p1 $p2 $wh >> $p1.aggregate
	wh=0;
	rm -rf $out $p2.$rand;
done
~bilebi00/bin/R --no-save --args $p1.aggregate < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/zscoreBins_v_foldChange.R > /dev/null

#rm -rf $p1.diff;
rm $p1.aggregate.testStat.pdf

echo DONE




