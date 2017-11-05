#!/bin/bash
#$ -S /bin/bash
# $ -q fs_long@@qc_nehalem
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_annotate_cls_format 
#$ -t 1-4

export PERL5LIB=$HOME/lib:$PERL5LIB

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> LOG._Get_annotate_cls_format 
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

prots=(EWSR1_1 EWSR1_2 EWSR1_1 EWSR1_2)
fts=(tag tag t2c t2c)
#TODO without filtering based on zScore the gene annotations might be wrong if the windows are all connected 
#TODO think a way how to eliminate this
i=$((SGE_TASK_ID - 1));
prot=${prots[$i]}
fn=${fts[$i]}
tn=${fn}_$prot

inp_src=~bilebi00/_EWSR1/data/clipz_CLIP/$prot/inputdata.$fn.sites.sorted 
pextLen=pextLen; #TODO set this
gif=~bilebi00/_EWSR1/data/hg19_TR.info.all
clsprog=~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel
trackfilespat=(~bilebi00/DATA/hg19_ucsc_tracks/rmsk.gtf* )
mrnafilespat=(~bilebi00/_EWSR1/data/hg19/coordinates/chr*)
outdir=CLIP_window_annot

mkdir -p $outdir; cd $outdir;

null=NA
inp=$tn.inp
annot=$tn.annot


#don't calculate if the datafile exists
if [ -e $dbin.annot ]; then
	exit;
fi

#Andreas' pipeline contains lots of one tag containing clusters
topN=`less $inp_src | wc -l | perl -e ' $_=<>;  print int($_ * 0.1), "\n"; '`;

less $inp_src | head -n $topN | awk '{ print tag"_"$0; }' tag=$tn > $inp;
allTracks=(${mrnafilespat[@]} ${trackfilespat[@]})
ac=0;
for ann in ${allTracks[@]}; do
	rs=(0);
	if [ `echo $ann | grep rmsk` ]; then rs=(0 1); fi
	for r in ${rs[@]}; do
		u=0;
		l=0;
		if [ `echo $ann | grep PromoterStart` ]; then u=$pextLen; fi
		if [ `echo $ann | grep TSR` ]; then u=$pextlen; fi
		cp $inp $ac.$tn.inp;
		atype=`~bilebi00/_DIS3/scripts/track2clsFormat.pl $ann 1 $u $l 0 $r 2>&1 >> $ac.$tn.inp`;
		echo Doing $ac $atype $ann
		~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $ac.$tn.inp > $ac.$tn.cls
		~bilebi00/_PAPD5/scripts/annot_from_cls.pl $tn $ac.$tn.cls 20 $null | awk 'BEGIN{OFS="\t";}{ if ($1 != null) print $2,atype,$1; }' atype=$atype null=$null > $ac.$tn.annot;
	rm $ac.$tn.inp $ac.$tn.cls
		ac=$((ac + 1))
	done
done
~bilebi00/_DIS3/scripts/select_annotation.pl $ac $tn.annot > final.$tn.annot
~bilebi00/_DIS3/scripts/assign_gene_id.pl final.$tn.annot $gif 2 NA > final2.$tn.annot
~bilebi00/_PAPD5/scripts/outerJoin.pl $inp final2.$tn.annot 1 1 '1-6,8-10' $null > $annot
rm *\.$tn.annot
rm $inp $annot

