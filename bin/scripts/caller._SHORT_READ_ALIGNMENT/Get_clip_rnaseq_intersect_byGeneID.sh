#!/bin/bash
#$ -S /bin/bash
# $ -q fs_short@@qc_nehalem
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.Get_clip_rnaseq_intersect_byGeneID_$TASK_ID
# $ -t 1-1
#$ -t 1-10
# $ -t 1-3
# $ -t 1-4
# $ -t 1-8

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

#TODO set
gi=gene; #splicing; #promoters; #gene; #isoform;
src="";
trusted=_trusted #"_combined"
clipfolder=~bilebi00/_DIS3/data/SP/*\.data;
#clipfolder=~bilebi00/_DIS3/4WEB_DIS3L2/PC*\.data
if [ "`perl -e '$cf=shift; $cf=~/data\/(SP)/; print "$1"; ' $clipfolder `" == "SP" ]; then 
	src=superClusters
fi
outdir=Stepanka/CLIP_comparison${trusted}_${src}_$gi
sa1=(totRNAseq-HEK-0_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-0_2 totRNAseq-HEK-0_2 totRNAseq-HEK-0 totRNAseq-HEK-0 mRNAseq-HEK-0 mRNAseq-HEK-0 totRNAseq-Hela-0)
sa2=(totRNAseq-HEK-0 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-siDis3L2 totRNAseq-Hela-siDis3L2)

#gi=splicing; #gene; #tss_group; #splicing; #promoters; #gene; #isoform;
#src="";
#trusted=_trusted; #"_combined"; # "_trusted"
#clipfolder=NOTEXIST; #TODO Could set Splicing ARRAY
#if [ "`perl -e '$cf=shift; $cf=~/data\/(SP)/; print "$1"; ' $clipfolder `" == "SP" ]; then 
#	src=superClusters
#fi
#outdir=Scheiffele/CLIP_comparison${trusted}_${src}_$gi
#sa1=(mRNAseq_Sam68_KO_brainStem mRNAseq_Sam68_KO_cerebellum mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_KO_cerebellum)
#sa2=(mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_KO_brainStem)

#gi=tss_group; #splicing; #promoters; #gene; #isoform;
#src="";
#trusted=_trusted; #"_combined"; # "_trusted"
#clipfolder=~bilebi00/_EWSR1/data/SP/*data
#if [ "`perl -e '$cf=shift; $cf=~/data\/(SP)/; print "$1"; ' $clipfolder `" == "SP" ]; then 
#	src=superClusters
#fi
#outdir=Shivendra/CLIP_comparison${trusted}_${src}_$gi
#sa1=(mRNAseq_siCTRL-A_HeLa)
#sa2=(mRNAseq_siEWSR1_HeLa)

#gi=splicing; #promoters; #gene; #isoform;
#src="";
#trusted=_trusted; #"_combined"; # "_trusted"
#clipfolder=~bilebi00/_KNOCKDOWNS/data/SP/*data
#if [ "`perl -e '$cf=shift; $cf=~/data\/(SP)/; print "$1"; ' $clipfolder `" == "SP" ]; then 
#	src=superClusters
#fi
#outdir=_Zavolan/CLIP_comparison${trusted}_${src}_$gi
#sa1=(mRNASeq_siRNA_GFP mRNASeq_No_4SU_No_XL_rep_A mRNASeq_No_4SU_No_XL_rep_B)
#sa2=(mRNASeq_siRNA_HuR mRNASeq_siRNA_HuR mRNASeq_siRNA_HuR)

#gi=splicing; #gene; #tss_group; #splicing; #promoters; #gene; #isoform;
#--------------------
cuffdiff_gexpd_file=../Cuffdiff${trusted}_RNAseq/${gi}_exp.diff
cuffdiff_gfpkm_file=../Cuffdiff${trusted}_RNAseq/${gi}s.fpkm_tracking
minfpkm=1
ps=10

if [ "$gi" == "promoters" ] || [ "$gi" == "splicing" ]; then
	cuffdiff_gexpd_file=../Cuffdiff${trusted}_RNAseq/${gi}.diff
	cuffdiff_gfpkm_file="";
	minfpkm=-1
fi

mkdir -p $outdir; pushd $outdir

s1=${sa1[$(($SGE_TASK_ID - 1))]}
s2=${sa2[$(($SGE_TASK_ID - 1))]}

p1=${s1}VS${s2}_$gi

less $cuffdiff_gexpd_file | perl -e '$s1=shift; $s2=shift; $minfpkm=shift; while(<>) { if ($_=~/\s$s1\s$s2\sOK\s/) { @t=split; print $_ if ($t[7]>$minfpkm && $t[8]>$minfpkm); }}' $s1 $s2 $minfpkm > $p1.diff.tmp

if [ "$gi" == "promoters" ] || [ "$gi" == "splicing" ]; then
	awk 'BEGIN{OFS="\t";} { print $0, $1;}' $p1.diff.tmp > $p1.diff
else
	#select only primary transcripts #TODO think; could it be a problem for EWSR1?
	~bilebi00/_PAPD5/scripts/innerJoin.pl $p1.diff.tmp $cuffdiff_gfpkm_file 1 1 '' | awk '$20 ~ /TSS/ { print $0;}' | cut -f 1-14,20 > $p1.diff
	~bilebi00/bin/R --no-save --args $p1.diff < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/cuffdiff_sig_scatter.R > /dev/null
	i=$p1.diff
	less $i | awk '$8>$9 && $14~/yes/{ print; }' > $i.pos; 
	less $i | awk '$8<$9 && $14~/yes/{ print; }' > $i.neg;  
	echo $i `less $i | awk '$14 ~ /yes/ { print; }' | wc -l` `less $i.pos | wc -l` `less $i.neg | wc -l`; 
fi

if [ "`echo $outdir | grep Scheiffele`" != "" ]; then
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
	~bilebi00/_PAPD5/scripts/leftJoin.pl $p1.diff $f2 3 16 '' 0 >> $out #genename join
#	~bilebi00/_PAPD5/scripts/innerJoin.pl $p1.diff $f2 3 16 '' >> $out #genename join
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/verify_clip_zscoreBins_with_foldChange.pl $out $minfpkm $region $ps $p1 $p2 $wh >> $p1.aggregate 
	wh=0;
	rm -rf $out;
done
~bilebi00/bin/R --no-save --args $p1.aggregate < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/zscoreBins_v_foldChange.R > /dev/null

#rm -rf $p1.diff;
rm $p1.diff.tmp

echo DONE




