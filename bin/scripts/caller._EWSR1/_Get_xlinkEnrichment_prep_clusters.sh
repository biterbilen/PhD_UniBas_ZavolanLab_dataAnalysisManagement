#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long
# $ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-8
#$ -o LOG._Get_xlinkEnrichment_prep_clusters.sh$TASK_ID

export PATH=$HOME/bin:$PATH
#----------------
indir=~bilebi00/_ARE/Analysis/UniqueRawData_ALL_woContamination
sid_prot_f=~bilebi00/_ARE/data/protein_sampleid_list
#----------------
mincoverage=1 #take all cluster where the depth is more then 1
tags=(`less $sid_prot_f | grep -w PARCLIP | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w PARCLIP | awk '$1 == tag { print $2; }' tag=$tag | sort `)
#XXX hard coded for 2 replicates
id1=${ids[0]}
id2=${ids[1]}

echo Doing $tag $id1 $id2
#----------------
outdir=Analysis/xlinkEnrichment/regions
mkdir -p $outdir; pushd $outdir;

#----------------
#clusters tags with bedtools merge --book-ended
for s in "-" "+"; do
	cat $indir/DB_all_$id1/copies_$s $indir/DB_all_$id2/copies_$s | \
		awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
		bedtools merge -s -i stdin -scores max -nms | \
		awk 'BEGIN{OFS="\t"; c=0; }{ if ($5>depth) { c=c+1; $4=c strand; print $1,src,fea,$2+1,$3,$5,strand,".","gene_id \""$4"\"; transcript_id \""$4"\";"; } }' src=$tag fea=cluster strand=$s depth=$mincoverage
done | gzip -c > ${tag}_clusters.gtf.gz

#----------------
#cluster length stats
echo -e "tag\ttype\tmu\tsd" > $tag.stats
zless ${tag}_clusters.gtf.gz  | awk '{ print tag"\t"type"\t"$5-$4+1; }' tag=$tag type=cluster_length | \
	  bedtools groupby -g 2 -c 3,3 -o mean,stdev -full | cut -f 1-2,4- >> $tag.stats

#----------------
#pernuc regions
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
bedtools nuc -fi $genomefa -bed ${tag}_clusters.gtf.gz -s -seq | cut -f 1-9,19 | awk 'NR>1{ print; }' | \
	gzip -c > ${tag}_clusters.gtf.nucbed.gz
#----------------
nucs="A C G T"
otag=_${tag}_clusters.gtf
~bilebi00/_EWSR1/scripts/split_gtf_nucbed.pl ${tag}_clusters.gtf.nucbed.gz "$nucs" $otag 0  
for nuc in $nucs; do
	gzip -f $nuc$otag;
done

echo DONE
#----------------

