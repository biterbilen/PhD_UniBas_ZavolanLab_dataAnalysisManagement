#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.Get_web_prep$TASK_ID 
#$ -t 1-2

export PERL5LIB=$HOME/lib:$PERL5LIB

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

outdir=4WEB

db=mm9
mrnadir=~bilebi00/_SAM68/data/coordinates/
gif=~bilebi00/_SAM68/data/mm_TR.info.all
clsprog=~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel
pat="^chr"
src=BinomDist
dbin=4db #hard coded in der.sql
note="Differential expression of overlapping 100mer genomic windows";

zf=zf
tns=(tns)
tds=(tds)
libs=(libs)

mkdir -p $outdir; cd $outdir;

#TODO without filtering based on zScore the gene annotations might be wrong if the windows are all connected 
#TODO think a way how to eliminate this
#----
if [ $SGE_TASK_ID == 1 ]; then
	#cerebellum
	zf=~bilebi00/_SAM68/zScore/377_376.zscore
	tn="SAM68WTvKOcerebellum";
	td="SAM68 WT v KO in Mouse Cerebellum"
	libs="376 377 634 635"
else
#----
#brain stem
	zf=~bilebi00/_SAM68/zScore/635_634.zscore
	tn="SAM68WTvKObrainStem"
	td="SAM68 WT v KO in Mouse Brain Stem"
	libs="376 377 634 635"
fi

gtf=$tn.gtf
inp=$zf.inp;
annot=$inp.annot
~bilebi00/_SAM68/scripts/zScoreFormat2clsFormat.pl $zf > $inp 
~bilebi00/_SAM68/scripts/annotate_region.pl "$pat" $inp $mrnadir $gif $clsprog > $annot
#~bilebi00/_SAM68/scripts/formatted_zscore2gtf.pl "$annot" "$tn" "$td" "$src" "$note" "$libs" $db > $gtf #FIXME too BIG to render
~bilebi00/_SAM68/scripts/formatted_zscore2der_table.pl $dbin.data $dbin.metadata "$annot" "$tn" "$td" "$src" "$note" "$libs" $db

rm $inp $annot



