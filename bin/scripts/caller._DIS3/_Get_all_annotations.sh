#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
# $ -o log_file
# $ -o LOG._Get_all_annotations_Normalization_PAPD5_T2C
# $ -o LOG._Get_all_annotations_Normalization_PAPD5
# $ -t 1-4
#$ -o LOG._Get_all_annotations_Normalization_PAPD5_copyAndT2C
#$ -t 1-2

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> log_file
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PERL5LIB=$HOME/lib:$PERL5LIB

#suff=.notfiltered
#filtfiles=(_04_14$suff _01_01$suff _02_12$suff _03_13$suff)
#1.
#odir=Normalization_PAPD5
#2.
#odir=Normalization_PAPD5_T2C
#3.
suff=.copy_t2c.notfiltered
filtfiles=(_04_04$suff _14_14$suff)
odir=Normalization_PAPD5_copyAndT2C

mkdir -p $odir; pushd $odir;

pextLen=3000
filtfile=${filtfiles[$((SGE_TASK_ID - 1))]}
minoverlap=20

tag=${filtfile/$suff/};
inp=$tag.inp
annot=${tag}.csv

trackfilespat=(~bilebi00/DATA/hg18_ucsc_tracks/refGene.gtf* ~bilebi00/DATA/hg18_ucsc_tracks/rnaGene.gtf* ~bilebi00/DATA/hg18_ucsc_tracks/snomir.gtf* ~bilebi00/DATA/hg18_ucsc_tracks/tRNAs.gtf* ~bilebi00/DATA/hg18_ucsc_tracks/refGenePromoterStart.gtf* ~bilebi00/DATA/hg18_ucsc_tracks/rmsk.gtf*)
#mrnafilespat=(~bilebi00/CDS_CLIP_ElMMo/data/coordinates/EXON/chr* ~bilebi00/CDS_CLIP_ElMMo/data/coordinates/INTRON/chr*)

rm -rf $annot;

f=''	
i=0
ps=11
~bilebi00/_DIS3/scripts/mergedzScoreFormat2clsFormat.pl $filtfile $tag 0 0 > $inp
for ann in ${trackfilespat[@]}; do
#for ann in refGene rnaGene snomir tRNAs rmsk refGenePromoterStart; do
	annbase=`basename $ann`;
	rs=(0);
	if [ `echo $ann | grep rmsk` ]; then rs=(0 1); fi
	for r in ${rs[@]}; do
		u=0;
		l=0;
		if [ `echo $ann | grep PromoterStart` ]; then u=$pextLen; fi
		if [ `echo $ann | grep TSR` ]; then u=$pextLen; fi
		out=${tag}_${annbase}_r${r};
		cp $inp $out.inp;
		atype=`~bilebi00/_DIS3/scripts/track2clsFormat.pl $ann 1 $u $l 0 $r 2>&1 >> $out.inp`;
		echo Doing $out
		~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $out.inp > $out.cls
		~bilebi00/_PAPD5/scripts/annot_from_cls.pl $tag $out.cls $minoverlap | awk 'BEGIN{OFS="\t";}{ if ($1 !~ /NULL/) print $2,atype,$1; }' atype=$atype > $out.annot;
		rm $out.inp $out.cls
		i=$((i+1))
		f=$f,$((i*3+ps))
		if [ -e $annot ] ; then
			~bilebi00/_PAPD5/scripts/outerJoin.pl $annot $out.annot 1 1 '' "NULL" > $annot.tmp
		else
			cp $out.annot $annot.tmp
		fi 
		rm $out.annot
		mv $annot.tmp $annot
	done
done

echo Merging...

echo -e "chromosome\tstrand\tbegin\tend\tz-score1\tz-score2\tUCSC RefSeq Genes Track\tUCSC Non-coding RNA Genes Track\tUCSC snoRNA and miRNA Genes Track\tUCSC tRNA Genes Track\tUCSC RefSeq Genes Promoters Track (${pextLen}KB)\tUCSC Repeat Masker Track\tUCSC Repeat Masker Track (Reverse)" > $annot.tmp
echo ~bilebi00/_PAPD5/scripts/outerJoin.pl $inp $annot 1 1 "2-7$f" "NULL"
~bilebi00/_PAPD5/scripts/outerJoin.pl $inp $annot 1 1 "2-7$f" "NULL" >> $annot.tmp
mv $annot.tmp $annot
rm $inp


