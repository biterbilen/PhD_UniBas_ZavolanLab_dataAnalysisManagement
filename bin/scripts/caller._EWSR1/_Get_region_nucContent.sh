#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-2
#$ -o LOG._Get_region_nucContent.sh$TASK_ID

export PATH=$HOME/bin:$PATH

outdir=Analysis/region_nucContent
mkdir -p $outdir; pushd $outdir;

if [ $SGE_TASK_ID == 1 ]; then
	db=hg19
	genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa

	intronf=~bilebi00/_EWSR1/Analysis/xlinkEnrichment/regions/hg19_GMAP_GENEintrons_refseq.gtf.gz
	clipf=~bilebi00/_EWSR1/Analysis/xlinkEnrichment/annotation/EWSR1/pernuc_merged_XL

	bedtools intersect -wa -s -a $intronf -b $clipf | sort | uniq > introns.clipped.gtf
	bedtools nuc -s -seq "-fi" $genomefa -bed introns.clipped.gtf > introns.clipped.nucbed

	n=`less introns.clipped.nucbed | wc -l`;
	bedtools subtract -s -a $intronf -b introns.clipped.gtf | head -n $n > introns.notclipped.gtf
	bedtools nuc -s -seq "-fi" $genomefa -bed introns.notclipped.gtf > introns.notclipped.nucbed
	touch prepdone
else
	while [ ! -e prepdone ]; do
		echo waiting for prep
		sleep 1m;
	done
fi

kmer=$SGE_TASK_ID
~bilebi00/_EWSR1/scripts/count_kmers.pl introns.clipped.nucbed $kmer gtf count > ${kmer}mers 2> ${kmer}mers.onames
~bilebi00/_EWSR1/scripts/count_kmers.pl introns.notclipped.nucbed $kmer gtf count > ${kmer}mers.background 2> /dev/null
R --no-save --args ${kmer}mers ${kmer}mers.onames ${kmer}mers.background < ~bilebi00/_EWSR1/scripts/kmer_enrichment_w_background_mu_ref.R > /dev/null
R --no-save --args ${kmer}mers ${kmer}mers.onames < ~bilebi00/_EWSR1/scripts/kmer_enrichment_w_background_mu_ref.R > /dev/null
for i in ${kmer}mers.*pbeta; do
	less $i | sort -k 5,5g -k 2,2gr > ${i}_
	mv ${i}_ $i;
done

if [ $SGE_TASK_ID ]; then
	rm -rf prepdone
fi

echo $kmer DONE
