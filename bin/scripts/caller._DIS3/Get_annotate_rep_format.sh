#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
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

allTracks=(annotfiles)
pextLen=pextLen;
clsprog=clsprog
tn=repf
outdir=output_directory

pushd $outdir;

null=NA

ac=$((SGE_TASK_ID - 1));
inp=$ac.$tn

~bilebi00/_EWSR1/scripts/rep2clsFormat.pl $tn $tn > $inp;
trackf=${allTracks[$ac]};

#rs=(0);
#if [ `echo $trackf | grep rmsk` ]; then rs=(0 1); fi
#for r in ${rs[@]}; do
r=0;
u=0;
l=0;
if [ `echo $trackf | grep PromoterStart` ]; then u=$pextLen; fi
if [ `echo $trackf | grep TSR` ]; then u=$pextlen; fi
atype=`~bilebi00/_DIS3/scripts/track2clsFormat.pl $trackf 1 $u $l 0 $r 2>&1 >> $inp`;
echo Doing $inp $atype $trackf
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $inp.cls
~bilebi00/_PAPD5/scripts/annot_from_cls.pl $tn $inp.cls 20 $null | awk 'BEGIN{OFS="\t";}{ if ($1 != null) print $2,atype,$1; }' atype=$atype null=$null > $inp.annot;
rm $inp $inp.cls

