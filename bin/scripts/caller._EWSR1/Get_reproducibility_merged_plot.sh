#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
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
exts=(exts);

pushd $outdir;

soutdir=Reproducibility

##~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_call_plot_correlation_dendrogram.sh  /import/bc2/soft/app/matlab/current/Linux/ combined.freq${nmer}mer $nmer

#plot the positional reproducibility
fti=12; lci=6; mti=5; s=""; h=1; 
for i in $soutdir/*rep; do 
	s="$s $i"; 
	~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/prep_positional_reproducibility_of_reps_forplotting.pl $i $h $fti $lci $mti $sid_prot_f; 
	h=0; 
done > $soutdir.posstat

~bilebi00/bin/R --no-save --args $soutdir.posstat < ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/positional_rep.R > /dev/null

for e in ${exts[@]}; do 
	nmer=$(($e * 2 + 1)); 
	~bilebi00/bin/R --no-save --args $soutdir/combined.freq${nmer}mer ${nmer}mer < ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/plot_correlation_matrix_dendrogram.R > /dev/null; 
done

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${soutdir}_combined_freq_kmer.pdf $soutdir/combined.freq*pdf > /dev/null

#`gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=OUT/_freq_dendrogram_$outdir.pdf $outdir/*freq_dendrogram.pdf`;
#`gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=OUT/_nc_content_$outdir.pdf $outdir/_nc_content*.pdf`;
#`gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=OUT/_positional_reproducibility_$outdir.pdf $outdir/_positional_reproducibility*.pdf`;

echo DONE

