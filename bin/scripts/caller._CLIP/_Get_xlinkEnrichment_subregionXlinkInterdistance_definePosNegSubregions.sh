#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 4-5
#$ -o LOG._Get_interdistance_of_xlinks_of_subregion.$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1";
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=$HOME/bin:$PATH
#-----------------
#XXX set
subregion=exon
subregion=intron
topN=10000
c=xlinkEnrichment
exptag=CTRL
exptag=HEK293
cell=HEK293; #""
subregionf=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/regions/hg19_GMAP_GENE${subregion}s_refseq.gtf.gz 

#-----------------

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") { print $2; }' tag=$tag cell="$cell" | sort`)

outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_${c}_topN${topN}_cell${cell}/$subregion
mkdir -p $outdir; pushd $outdir;

less ../venn.$subregion.txt | perl -e 'while(<>){ chomp;@t=split; $g=shift @t; $a=0; for (@t) { $a+=$_;} print "$g\n" if ($a == @t); }' > common_subregion_genes

suff=bed.w_common_subregion_genes
for id in ${idsf[@]}; do
	f="../$id.bed.w_type_genename.gz"
	~bilebi00/_xlinkEnrichment/scripts/grep_f_w_like_w_column.pl common_subregion_genes $f 7 | grep $subregion > $id.$suff
done

fisherf=subregion_xlinkCount.fisher_stats
echo -e "lib\ttype\tvalue" > $fisherf
statf=subregion_xlinkCount.stats
echo -e "lib\ttype\tvalue\tN" > $statf
for id in ${idsf[@]}; do
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	#	cut=6
	#	n=`bedtools closest -io -s -d -t "first" -a $id.$suff -b $id.$suff | awk '$17<cut{print}' cut=$cut | wc -l`
	#	N=`less $id.$suff | wc -l`;
	#	echo | awk 'BEGIN{OFS="\t"}{ print lib,type,n/N*100;}' lib="$idn (N=$N)" type=intradistance$cut n=$n N=$N >> $statf 
	#interdistance between the samples
	for id2 in ${idsf[@]}; do
		if [ $id2 -le $id ]; then continue; fi
		idn2="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id2`";
		lib="$idn and $idn2"
		s1=`less $id.$suff | wc -l`;
		s2=`less $id2.$suff | wc -l`;
		u=`cat $id.$suff $id2.$suff | cut -f 1-3,6 | sort | uniq | wc -l`;
		Ns=$s2
		Nb=$s1
		for dist in 0 10 100 1000; do 
			n=`bedtools closest -s -d -t "first" -a $id2.$suff -b $id.$suff | awk '$17<=cut{print}' cut=$dist | wc -l`
			n2=`bedtools closest -s -d -t "first" -a $id.$suff -b $id2.$suff | awk '$17<=cut{print}' cut=$dist | wc -l`
			if [ $n2 -lt $n ]; then
				n=$n2;
			fi
			echo | awk 'BEGIN{OFS="\t";}{ print lib,type,n/N*100,N;}' lib="$lib" type="Distance=$dist" n=$n N=$Ns >> $statf;
#			u=10000 
			echo -e $lib" d=<$dist\t"`R  --no-save --args  $n $Ns $Nb $u < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"` >> $fisherf 
		done
	done
done

topN=100
stack=0
R --no-save --args $statf $topN $stack < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null

#TODO write clipped unCLIPPED CODE here
for id in ${idsf[@]}; do
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	#Find crosslinked regions using top crosslinked positions
	~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl $subregionf common_subregion_genes | \
			bedtools intersect -s -u -a stdin -b $id.$suff > $id.CLIPPED
	#~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl $subregionf common_subregion_genes | \
	#		bedtools intersect -s -v -a stdin -b $id.$suff > $id.NOTCLIPPED

	#Find not clipped regions on clipped genes using all crosslink positions
	~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl $subregionf common_subregion_genes | \
			bedtools intersect -s -v -a stdin -b ../$tag.$id.$id.T.xlink_mu_refmax.bed.gz | awk '$5-$4>200{print}' > $id.NOTCLIPPED
	echo negative and positive sets prep
done

echo prepared CLIPPED and NOTCLIPPED regions
for suff in CLIPPED NOTCLIPPED; do
	files="";
	for id in ${idsf[@]}; do
		files="$files $id.$suff"
	done
	cat $files | sort | uniq -D | uniq > joint.$suff
	echo joint.$suff exists
done

echo prepared joint CLIPPED and NOTCLIPPED regions

echo DONE
