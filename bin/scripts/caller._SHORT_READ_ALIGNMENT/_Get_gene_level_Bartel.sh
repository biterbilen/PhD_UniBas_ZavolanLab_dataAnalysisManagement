#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_gene_level_Bartel.sh 

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

outdir=Guo

pushd $outdir

#1.
g1=mrna_mir1_32hr
g2=mrna_mock_32hr
mir=hsa-miR-1

#2.
g1=mrna_mir155_12hr
g2=mrna_mock_12hr
mir=hsa-miR-155

less Cuffdiff_mrna/gene_exp.diff | grep $g1 | grep $g2 > $g1$g2 
~bilebi00/_PAPD5/scripts/innerJoin.pl ../data/Guo_et_al_2010_nature_GSE21992/*${g2}_quantification.txt.gz $g1$g2 2 3 '3,7-20' > $g2-$g1$g2
~bilebi00/_PAPD5/scripts/innerJoin.pl ../data/Guo_et_al_2010_nature_GSE21992/*${g1}_quantification.txt.gz $g2-$g1$g2 2 4 '1-3,7-21' > $g1$g2-$g1$g2 
header="representative_transcript gene_id bartelRPKM1 bartelRPKM2";
header="$header "`head -n 1 Cuffdiff_mrna/gene_exp.diff`;
header="$header "`head -n 1 $mir.seed_count`

echo $header | sed -e 's/ /\t/g' > $g1$g2-$g1$g2.seeds;
~bilebi00/_PAPD5/scripts/innerJoin.pl $g1$g2-$g1$g2 $mir.seed_count 1 1 '' >> $g1$g2-$g1$g2.seeds

/import/bc2/home/zavolan/bilebi00/bin/R --no-save --args $g1$g2-$g1$g2.seeds ${g1}VS${g2} < /import/bc2/home/zavolan/bilebi00/_SHORT_READ_ALIGNMENT/scripts/gene_level_Bartel.R > /dev/null

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${g1}VS${g2}.pdf ${g1}VS${g2}_*.pdf

rm $g1$g2 $g2-$g1$g2 $g1$g2-$g1$g2 
