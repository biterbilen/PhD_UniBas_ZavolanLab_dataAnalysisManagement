#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
# $ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-runtimes

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

zfs=(zfs)
tns=(tns)
tds=(tds)
alllibs=(libs)
zscorecut=zscoreCut
pextLen=pextLen;
indir=indir
db=db
gif=gif
clsprog=clsprog
src=src
trackfilespat=trackfilesspat
mrnafilespat=mrnafilespat
note=note
outdir=output_directory

mkdir -p $outdir; cd $outdir;

#TODO without filtering based on zScore the gene annotations might be wrong if the windows are all connected 
#TODO think a way how to eliminate this
i=$(($SGE_TASK_ID - 1));
zf=${zfs[$i]}
tn=${tns[$i]}
td="${tds[$i]}"
libs="${alllibs[$i]}"

null=NA
inp=$tn.inp
annot=$tn.annot
dbin=$tn

#don't calculate if the datafile exists
if [ -e $dbin.data ]; then
	exit;
fi

~bilebi00/_DIS3/scripts/mergedzScoreFormat2clsFormat.pl $indir$zf $tn $zscorecut -$zscorecut > $inp;
allTracks=(${trackfilespat[@]} ${mrnafilespat[@]})
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
~bilebi00/_PAPD5/scripts/outerJoin.pl $inp final2.$tn.annot 1 1 '1-11,13-15' $null > $annot
rm *.$tn.annot
~bilebi00/_DIS3/scripts/formatted_mergedzScore2der_table.pl $dbin.data $dbin.metadata $annot $tn "$td" $src "$note" "$libs" $db

rm $inp $annot

echo DONE
