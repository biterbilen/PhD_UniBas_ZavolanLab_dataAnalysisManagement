#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-runs


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
mkdir -p $outdir; cd $outdir;

#infiles="DB_69-DB_70 DB_59-DB_60";
#exts=(2   3  5 8  10 15 20);
#begs=(18 17 15 12 10  5  0);
infiles="in_files";
exts=(exts);
begs=(begs);

#TODO split here
i=$(($SGE_TASK_ID -1))
for infile in $infiles; do
	folders=`echo $infile | sed -e 's/-/ /g'`;
	fieldindex=12;
	for f in $folders; do
		~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/count_substring.pl $infile.rep $fieldindex ${begs[$i]} ${exts[$i]} 2 > $f.freq${exts[$i]}
		fieldindex=$(($fieldindex + 1));
	done
done
nmer=$((${exts[$i]} * 2 + 1));
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/get_combined_motif_freqs.pl . .freq${exts[$i]} > combined.freq${nmer}mer
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_plot_correlation_dendogram.sh  /import/bc2/soft/app/matlab/current/Linux/ combined.freq${nmer}mer $nmer
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_call_plot_correlation_dendogram.sh  /import/bc2/soft/app/matlab/current/Linux/ combined.freq${nmer}mer $nmer
rm *.freq${exts[$i]}


