#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_seed_counts_and_cuffdiffDE.sh 
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

PATH=$HOME/bin:$PATH
export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=Guo; pushd $outdir;

mirs=(hsa-miR-1 hsa-miR-155);

mir=${mirs[$((SGE_TASK_ID - 1))]}

if [ `echo $mir| grep miR-155` ]; then 
	miralias=mir155; 
elif [ `echo $mir| grep miR-1` ]; then
	miralias=mir1
else
	echo Unknown mir: $mir;
	exit;
fi

pairs=("mrna_${miralias}_12hr mrna_mock_12hr" "mrna_${miralias}_32hr mrna_mock_32hr")

for pair in "${pairs[@]}"; do
pair="mrna_${miralias}_12hr mrna_mock_12hr"
	for cat in isoforms genes; do 
	cat=isoforms
		out=$miralias.`echo $pair | sed -e 's/ /VS/g'`.$cat
		echo $out;
		infile=Cuffdiff/`echo $cat | sed 's/s$//g'`_exp.diff;
		infile2=Cuffdiff/`echo $cat | sed 's/s$//g'`_exp.diff;
		tmp=tmp$RANDOM
		less $infile | perl -e '$g1=shift; $g2=shift; $_=<>; print $_; while(<>) { @t=split; next if ($t[6] ne "OK" or $t[10] ne "OK"); if ($t[4] eq $g1 and $t[5] eq $g2) { print $_;} }' $pair > $tmp; 
		cat $tmp >> $out;
		echo "R --no-save --args $out < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/FPKM_scatter.R" 
		R --no-save --args $out < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/FPKM_scatter.R 
		rm $tmp
	done
	#TODO redo with barchart with errorbar using ggplot library
	R --no-save --args $mir.seeds.isoforms < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/real_bartel_data_plots.R
done




