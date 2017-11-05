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
pushd $outdir;

soutdir=Reproducibility
mkdir -p $soutdir; pushd $soutdir;

infiles="in_files";
exts=(exts);
begs=(begs);

seqindex=2
libcountindex=6
maxlibcount=1; #if set to one check only in singleton positions not reproduced in other replicated libraries
#TODO split here
i=$(($SGE_TASK_ID -1))
for infile in $infiles; do
	folders=`echo $infile | sed -e 's/-/ /g'`;
	fieldindex=12;
	for f in $folders; do
		~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/count_substring.pl $infile.rep $fieldindex ${begs[$i]} ${exts[$i]} $seqindex $libcountindex $maxlibcount > $f.freq${exts[$i]}
		fieldindex=$(($fieldindex + 1));
	done
done

nmer=$((${exts[$i]} * 2 + 1));
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/get_combined_motif_freqs.pl . .freq${exts[$i]} > combined.freq${nmer}mer
rm *.freq${exts[$i]}

echo DONE

