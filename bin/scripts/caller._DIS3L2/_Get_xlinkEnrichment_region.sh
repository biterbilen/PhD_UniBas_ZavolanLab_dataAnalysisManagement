#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-24
#$ -o LOG._Get_xlinkEnrichment_region.sh.clusters._tagorder.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_region.sh.introns._tagorder.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_region.sh.introns._tagorder.$TASK_ID

#HOWTO for submission TODO set
#c=exons; s=_Get_xlinkEnrichment_region.sh; for i in `seq 1 8`; do sed -e "s/_tagorder/$i/" $s > $s$i.$c; qsub $s$i.$c; echo $i; done
#c=clusters; s=_Get_xlinkEnrichment_region.sh; for i in `seq 1 8`; do sed -e "s/_tagorder/$i/" $s > $s$i.$c; qsub $s$i.$c; echo $i; done
#c=introns; s=_Get_xlinkEnrichment_region.sh; for i in `seq 1 8`; do sed -e "s/_tagorder/$i/" $s > $s$i.$c; qsub $s$i.$c; echo $i; done

#HOWTO check runs in a batch?
#lines=7; grep DONE LOG._Get_xlinkEnrichment_region.sh* | awk  -F ":" '{ print $1; }' | while read i; do wc -l $i; done | grep "^$lines " | awk '{ print $2; }' | while read i; do rm $i; done

#XXX hard coded for 2 parclip as foreground and 4 rnaseq as background
#for i in LOG._Get_xlinkEnrichment.sh.289.pernuc.*; do echo $i `less $i | grep -w waiting -v | wc -l`; done | grep -v -w 7
#-----------------
export PATH=$HOME/bin:$PATH
#-----------------

tagorder=_tagorder; #1..8; 4.EWSR1 8.TIA1

sid_prot_f=~bilebi00/_DIS3L2/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((tagorder-1))]}

#TODO set
c=exons_refseq; regionfile=../../regions/hg19_GMAP_GENE${c}.gtf.gz
c=clusters; regionfile=../../regions/${tag}_$c.gtf.gz
#c=introns_refseq; regionfile=../../regions/hg19_GMAP_GENE${c}.gtf.gz

#-----------------
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk '$1==tag { print $2; }' tag=$tag | sort` \
 	`less $sid_prot_f | grep -w mRNAseq | grep -w background | cut -f 2 | sort`)
nucs=(A C G T); 

#-----------------
outdir=Analysis/xlinkEnrichment/basedOn$c/$tag
mkdir -p $outdir; pushd $outdir;
#-----------------
i=`echo $((SGE_TASK_ID-1)) / ${#ids[@]} % ${#nucs[@]} | bc`;
nuc=${nucs[$i]};

i=`echo $((SGE_TASK_ID-1)) / 1          % ${#ids[@]}  | bc`;
id=${ids[$i]};
#-----------------
r=$id.$c.$nuc
#-----------------

for mut in xlink mut; do
	if [ $mut == xlink ] && [ $nuc != T ]; then
		echo `date` $mut for $nuc skipped
		echo `date` Copies and $mut statistics skipped
		continue;
	fi 
	if [ ! -e $tag.$r.${mut}_count.gtf.gz ]; then
		f2=../../basedOnNucleotide/$tag/$tag.$id.clusters.$nuc.${mut}_count.gtf.gz
		bedtools intersect -wao -s -a $regionfile -b $f2 | awk -F "\t" '$19 > 0 {print;}' | \
			sort -k 1,9 | sed 's/copies.*"\([0-9]*\)"; mut_count.*"\([0-9]*\)";/\t\1\t\2/' | \
			bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 19,20 -o sum,sum | \
			awk -F "\t" 'BEGIN {OFS="\t";} { $9=$9" copies_count \""$10"\"; mut_count \""$11"\";"; if ($10 > 0) { print; } }' | \
			cut -f 1-9 | \
			gzip -c > $tag.$r.${mut}_count.gtf.gz
		echo `date` $tag.$r.${mut}_count.gtf.gz done
	else
		echo `date` $tag.$r.${mut}_count.gtf.gz exists
	fi

	if [ ! -e $tag.$r.${mut}_count.freqstat ]; then
		zless $tag.$r.${mut}_count.gtf.gz | \
			perl -e '$tag=shift; $type=shift; while(<>){ $_ =~ /copies_count "(\d+).*mut_count "(\d+)";$/; print join("\t", $tag, $type, $2/$1),"\n"; }' $tag.$r.$mut freq | \
			bedtools groupby -g 1,2 -c 3,3,3 -o mean,stdev,count > $tag.$r.${mut}_count.freqstat
		echo `date` Copies and $mut statistics done
	else
		echo `date` $tag.$r.${mut}_count.freqstat exists
	fi
done

echo `date` Cleaning done
echo DONE


