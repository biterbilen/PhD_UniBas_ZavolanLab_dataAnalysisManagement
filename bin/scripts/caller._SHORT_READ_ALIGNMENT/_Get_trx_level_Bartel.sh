#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_trx_level_Bartel.sh 

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

outdir=Guo/Trx_level
mkdir -p $outdir; pushd $outdir

tag=isoform
#1.
g1=mrna_mir1_12hr
g2=mrna_mock_12hr
mir=hsa-miR-1

#2.
g1=mrna_mir155_32hr
g2=mrna_mock_32hr
mir=hsa-miR-155

less ../Cuffdiff_mrna/${tag}_exp.diff | grep $g1 | grep $g2 > _$g1$g2 
~bilebi00/_PAPD5/scripts/innerJoin.pl _$g1$g2 ../Cuffdiff_mrna/${tag}s.fpkm_tracking 1 1  '1-4,16,17,7-14' | awk '$5 ~ /=/ { print; }' > $g1$g2

~bilebi00/_PAPD5/scripts/innerJoin.pl ~bilebi00/_SHORT_READ_ALIGNMENT/data/Guo_et_al_2010_nature_GSE21992/*${g2}_quantification.txt.gz $g1$g2 2 3 '3,7-20' > $g2-$g1$g2
~bilebi00/_PAPD5/scripts/innerJoin.pl ~bilebi00/_SHORT_READ_ALIGNMENT/data/Guo_et_al_2010_nature_GSE21992/*${g1}_quantification.txt.gz $g2-$g1$g2 2 4 '1-3,7-21' > $g1$g2-$g1$g2 
header="representative_transcript gene_id bartelRPKM1 bartelRPKM2";
header="$header "`head -n 1 ../Cuffdiff_mrna/gene_exp.diff`;
header="$header "`head -n 1 ../$mir.seed_count`

echo $header | sed -e 's/ /\t/g' > $g1$g2-$g1$g2.seeds;
tagi=10; #1 in gene case
~bilebi00/_PAPD5/scripts/innerJoin.pl $g1$g2-$g1$g2 ../$mir.seed_count $tagi 1 '' >> $g1$g2-$g1$g2.seeds

/import/bc2/home/zavolan/bilebi00/bin/R --no-save --args $g1$g2-$g1$g2.seeds $tag.${g1}VS${g2} < /import/bc2/home/zavolan/bilebi00/_SHORT_READ_ALIGNMENT/scripts/gene_level_Bartel.R > /dev/null

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$tag.${g1}VS${g2}.pdf $tag.${g1}VS${g2}_*cufflinks*.pdf

rm _$g1$g2 $g1$g2 $g2-$g1$g2 $g1$g2-$g1$g2 
