#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file

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

outdir=output_directory
shuffleN=shuffleN
maxw=maxw
tag=tag

pushd $outdir/MEME

odir=$tag


for i in `seq 0 $shuffleN`; do
	mdir=$odir
	if [ $i -gt 0 ]; then mdir=$odir.$i; fi
	chmod 755 $mdir
	echo `less $mdir/meme.html | grep "^letter-probability matrix" -A $maxw`
#TODO Integrate but don't run motif searh again (???) :
#less $mdir/meme.txt | perl -e '$d=0; $header=shift; $h=0; while(<>) { $h=1 if ($_ =~ /$header/); next unless ($h); $d=1 if ($_ =~ /^seq/); last if ($d==1 and $_=~/^\-/);  if($_=~/^seq/) { (@c) = ($_=~/(\[)/g); $h{scalar @c}++; } } for (sort {$h{$b} <=> $h{$a}} keys %h) { print $_, "\t", $h{$_}, "\n"; } ' "Motif 1 block diagrams"
#topN=200
#1 66
#2 21
#4 1
#3 1
	if [ $i -gt 0 ]; then rm -rf $mdir; fi
done > $odir/allstats 
awk '{ print $8, $10; }' $odir/allstats > $odir/stats
~bilebi00/bin/R --no-save --args $odir/stats < ~bilebi00/_EWSR1/scripts/plot_MEME_stats.R > /dev/null 
/usr/bin/convert $odir/logo1.eps $odir/tmp.logo.pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$odir/tmp.logo.pdf $odir/logo1.eps 
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$odir/summary.pdf $odir/tmp.logo.pdf $odir/stats.pdf
rm $odir/tmp*

#put in the web folder
#chmod 755 $odir
#pushd ~bilebi00/www/MEME/EWSR1/;
#ln -sf $pwd/$outdir/$soutdir/$odir .
#-------

echo DONE;
exit;
