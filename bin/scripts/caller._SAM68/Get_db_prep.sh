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

indir=indir
zfs=(zfs)
tns=(tns)
tds=(tds)
alllibs=(libs)
db=db
gif=gif
clsprog=clsprog
src=src
mrnadir=mrnadir
regionlinepat=reglinepat
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

gtf=$tn.gtf
inp=$tn.inp;
annot=$tn.annot
dbin=$tn

#don't calculate if the datafile exists
if [ -e $dbin.data ]; then
	exit;
fi

~bilebi00/_SAM68/scripts/zScoreFormat2clsFormat.pl $indir$zf > $inp 
~bilebi00/_SAM68/scripts/annotate_region.pl "$regionlinepat" $inp $mrnadir $gif $clsprog > $annot
#~bilebi00/_SAM68/scripts/formatted_zscore2gtf.pl "$annot" "$tn" "$td" "$src" "$note" "$libs" $db > $gtf #FIXME too BIG to render
~bilebi00/_SAM68/scripts/formatted_zscore2der_table.pl $dbin.data $dbin.metadata $annot $tn "$td" $src "$note" "$libs" $db

rm $inp $annot



