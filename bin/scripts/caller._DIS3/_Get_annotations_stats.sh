#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
# $ -o log_file
#$ -o LOG._Get_annotations_stats_Normalization_PAPD5_T2C
# $ -o LOG._Get_annotations_stats_Normalization_PAPD5
#$ -t 1-4
# $ -o LOG._Get_annotations_stats_Normalization_PAPD5_copyAndT2C
# $ -t 1-2

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

export PERL5LIB=$HOME/lib:$PERL5LIB

#1.
suff=.notfiltered
cuts=(1 1.5 2)
cuts2=(1 1.5 2)
odir=Normalization_PAPD5_T2C
filtfiles=(_04_14$suff _01_01$suff _02_12$suff _03_13$suff)

#2.
#suff=.notfiltered
#cuts=(3 4 5 6)
#cuts2=(3 4 5 6)
#odir=Normalization_PAPD5
#filtfiles=(_04_14$suff _01_01$suff _02_12$suff _03_13$suff)

#3.
#suff=.copy_t2c.notfiltered
#cuts=(4    5)
#cuts2=(1.5 1.5)
#odir=Normalization_PAPD5_copyAndT2C
#filtfiles=(_04_04$suff _14_14$suff)

mkdir -p $odir; pushd $odir;

filtfile=${filtfiles[$((SGE_TASK_ID - 1))]}

tag=${filtfile/$suff/};
annot=${tag}.csv

for i in `seq 0 $((${#cuts[@]}-1))`; do
	c=${cuts[$i]}
	c2=${cuts2[$i]}
	if [ -e $annot ]; then 
		echo -e Total Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($12 != "NULL" || $13 != "NULL")) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_rRNA Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($12 != "NULL" || $13 != "NULL") && ($12 ~ /rRNA/ || $13 ~ /rRNA/ || $7 ~ /rRNA/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_SINE Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($12 ~ /Alu/ || $12 ~ /MIR/ || $13 ~ /Alu/ || $13 ~ /MIR/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_LINE Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($12 ~ /^L[1234]/ || $12 ~ /HAL/ || $13 ~ /^L[1234]/ || $13 ~ /HAL/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" mir Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($9 ~ /^hsa/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" snoRNA Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && (($8 ~ /^U[0-9]/ || $8 ~ /^HB/ || $9 ~ /^U[0-9]/ || $9 ~ /^HB/) && $12 ~ /NULL/ && $13 ~ /NULL/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" mRNA_PROMPT Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($11 ~ /NM_/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_tRNA Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($10 ~ /tRNA/ && ($12 ~ /tRNA/ || $13 ~ /tRNA/))) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" tRNA Enriched `less $annot | awk 'NR>1{ if ($5>c && $6>c2 && ($10 ~ /tRNA/ && ($12 ~ /NULL/ || $13 ~ /NULL/))) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n"Total Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($12 != "NULL" || $13 != "NULL")) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_rRNA Depleted  `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($12 != "NULL" || $13 != "NULL") && ($12 ~ /rRNA/ || $13 ~ /rRNA/ || $7 ~ /rRNA/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_SINE Depleted  `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($12 ~ /Alu/ || $12 ~ /MIR/ || $13 ~ /Alu/ || $13 ~ /MIR/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_LINE Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($12 ~ /^L[1234]/ || $12 ~ /HAL/ || $13 ~ /^L[1234]/ || $13 ~ /HAL/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" mir Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($9 ~ /^hsa/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" snoRNA Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && (($8 ~ /^U[0-9]/ || $8 ~ /^HB/ || $9 ~ /^U[0-9]/ || $9 ~ /^HB/) && $12 ~ /NULL/ && $13 ~ /NULL/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" mRNA_PROMPT Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($11 ~ /NM_/)) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" Repeat_and_tRNA Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($10 ~ /tRNA/ && ($12 ~ /tRNA/ || $12 ~ /tRNA/))) print $0; }' c=$c c2=$c2 | wc -l` \
			"\n" tRNA Depleted `less $annot | awk 'NR>1{ if ($5<-c && $6<-c2 && ($10 ~ /tRNA/ && ($12 ~ /NULL/ || $13 ~ /NULL/))) print $0; }' c=$c c2=$c2 | wc -l` \
			> $annot.$c.$c2.stat
		echo Stats are in $annot.$c.stat
		~bilebi00/bin/R --no-save --args $annot.$c.$c2.stat < ~bilebi00/_DIS3/scripts/fisher.ex.t.R  > /dev/null
	else
		echo $annot does not exist!
	fi 
done


