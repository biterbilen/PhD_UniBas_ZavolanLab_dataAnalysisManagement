#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_seed_counts.sh 
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

less ~bilebi00/CDS_CLIP_ElMMo/data/hsa.fa | grep -w -A 1 ">$mir" > $mir.fa
for dom in CDS 3UTR; do
~bilebi00/CDS_CLIP_ElMMo/scripts/count_mir_fam_seed_matches_in_fasta.pl $mir.fa ~bilebi00/CDS_CLIP_ElMMo/data/hg_${dom}_GMAP.all.fa > $mir.$dom.seeds 
done

tmp=tmp$RANDOM
~bilebi00/_PAPD5/scripts/innerJoin.pl $mir.3UTR.seeds $mir.CDS.seeds 1 1 "1,2,4" | awk -F "|" 'BEGIN{OFS="\t";}{if(NR>1){print $1,$2;}}' > $tmp
echo "id coord 3UTR.seedCount CDS.seedCount seedType" | sed -e 's/ /\t/g' > $mir.seed_count;
less $tmp | ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/setSeedType.pl >> $mir.seed_count

#clean
for dom in CDS 3UTR; do
	rm $mir.$dom.seeds
done

rm $mir.fa $tmp




