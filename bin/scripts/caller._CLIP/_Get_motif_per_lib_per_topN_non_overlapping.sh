#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-3
#$ -o LOG._Get_motif_per_lib_per_topN_non_overlapping.sh$TASK_ID

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

#TODO set tag and LOG file above to not to mix the results
tag=EWSR1
#tag=EWSR1_FlpIn
tag=DIS3L2
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
N=2000
n=20
#TODO integrate this option in the script
omit_overlapping_dont_sum_score=0
b=neighbour
cell="" #HEK293

region_len=$((n*2+1))
idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && (cell == "" || $7==cell) {print $2;}' tag=$tag cell="$cell" | sort | uniq`)

if [ $SGE_TASK_ID -gt ${#idsf[@]} ]; then
	exit;
fi

date
echo ${idsf[@]}

#TODO delete singleton
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/motif_topN_non_overlapping_singleton
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/motif_topN_non_overlapping_representative_skewed_score_double
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/motif_topN_non_overlapping_representative_strong_skewed_score
outdir=Analysis/Project_$tag/xlinkEnrichment/$tag/motif_topN_non_overlapping_representative_strong_skewed_score
mkdir -p $outdir; pushd $outdir;

#-----------------DATA PREP
id=${idsf[$((SGE_TASK_ID-1))]}
file=../xlinkEnrichment/$tag.$id.$id.T.xlink_mu_refmax.bed.gz
ln -sf $file .

allCCRs=all$id.CCRs.bed.gz
no_allCCRs=no_$allCCRs #non-overlapping
extrafields=$((`zless $file | awk 'NR==1{ print NF;}'`-6))

if [ ! -e $allCCRs ]; then
	bedtools slop -b $n -i $file -g $genomelen | gzip -c > $allCCRs
else
	echo $allCCRs exists
fi

if [ ! -e $no_allCCRs ]; then
	if [ $omit_overlapping_dont_sum_score == 1 ]; then
		bedtools intersect -wao -s -a $allCCRs -b $allCCRs | \
			~bilebi00/_xlinkEnrichment/scripts/get_nonoverlapping_from_top.pl $omit_overlapping_dont_sum_score | \
			sort -k 5,5g | gzip -c > $no_allCCRs
	else 
		bedtools intersect -wao -f 0.5 -s -a $allCCRs -b $allCCRs | \
			~bilebi00/_xlinkEnrichment/scripts/get_nonoverlapping_from_top.pl | \
			sort -k 5,5g | gzip -c > $no_allCCRs
	fi
else
	echo $no_allCCRs exists
fi

if [ ! -e foreground$id.nucbed.gz ]; then
	zless $no_allCCRs | awk '$3-$2 == len { print; }' len=$region_len | \
		head -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin | gzip -c > foreground$id.nucbed.gz
else
	echo foreground exists
fi

if [ ! -e background$id.nucbed.gz ]; then
	if [ $b == "bottom" ]; then
		zless $no_allCCRs | awk '$3-$2 == len { print; }' len=$region_len | \
			tail -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin gzip -c > background$id.nucbed.gz
	elif [ $b == "neighbour" ]; then
		zless $allCCRs | bedtools flank -b $region_len -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$region_len | \
			bedtools intersect -v -s -a stdin -b $allCCRs -f 0.05 | head -n $N | \
			bedtools nuc -s -seq "-fi" $genomefa -bed stdin | gzip -c > background$id.nucbed.gz
	fi
else
	echo background exists
fi

echo Doing motifs
mkdir -p $id; pushd $id
zless ../foreground$id.nucbed.gz | head -n $N | awk 'BEGIN{c=0;} { if($1 ~ /chr/){ c=c+1; print ">seq"c"\n"$(16+ef); }}' ef=$extrafields > foreground.fa
zless ../background$id.nucbed.gz | awk 'BEGIN{c=0;} { if($1 ~ /chr/){ c=c+1; print ">seq"c"\n"$(16+ef); }}' ef=$extrafields > background.fa

N=$((`zless foreground.fa | grep "^>" | wc -l`));

nmotifs=2
minw=5
maxw=10
mod=zoops

minsites=`echo "$N * 30 / 100" | bc`;
maxsites=$N;

#tried p=50:10:80 and the motif GAAGAA or GG[AT]GG[AT] for subregion=exon
#if [ $tag == EWSR1 ]; then mod=anr; minsites=`echo "$N * 80 / 100" | bc`; fi
# if [ $tag == EWSR1 ]; then minsites=`echo "$N * 80 / 100" | bc`; fi
#if [ $mod == "anr" ]; then maxsites=$((N * 2)); fi

m=$(($minw-4))
less background.fa | fasta-get-markov -m $m -norc > markov$m

echo -e "\n" meme foreground.fa -p 8 -maxsize 1000000 -sf $tag -dna -bfile markov$m -mod $mod -nmotifs $nmotifs -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc ./ -nostatus
meme foreground.fa -p 8 -maxsize 1000000 -sf $tag -dna -bfile markov$m -mod $mod -nmotifs $nmotifs -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc ./ -nostatus > /dev/null
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$tag.logo.pdf logo1.eps
popd

touch runfinished.$SGE_TASK_ID

if [ $SGE_TASK_ID == 1 ]; then
	while [ `ls runfinished* | wc -l` -lt ${#idsf[@]} ]; do
		echo waiting for runs;
		sleep 1m
	done
	#check crosslink distance
	for tag in ${idsf[@]}; do 
		idn="`less $sid_prot_f | awk -F "\t" 'id~$2{print $6; }' id=$tag`"; 
		less $tag/meme.txt | perl -e '$tag=shift; $start=0; while(<>){ if ($_ =~ /^Sequence name\s*Start/) { $_=<>; $start=1; next; } if ($_=~/^\-/ && $start==1) {exit; } if ($start) { @t=split/\W+/;print $t[1]."\t".$tag."\n";} }' "$idn" | \
			sort | bedtools groupby -g 1 -c 1 -full -o count -i stdin | sort -k 1,1g | \
	 		awk -F "\t" 'BEGIN{OFS="\t";}{ if (c==0){c=$1+1;} for(i=c;i<$1;i++) { print c,$2,0; c=c+1; } print $0; c=$1+1; }'; 
	done | awk -F "\t" 'BEGIN{OFS="\t";} { $1=$1-region_len; print $0; }' region_len=$((n+1)) > motif_distance
	rows=2
	if [ ${#idsf[@]} -lt 3 ]; then
		rows=1
	fi
	R --no-save --args motif_distance $rows < ~bilebi00/_ARE/scripts/plot_profile_4Groups_4Motif.R > /dev/null

	if [ $tag == EWSR1 ]; then
		less motif_distance | grep -P "80kDa|120kDa" > motif_distance.288.289
		R --no-save --args motif_distance.288.289 2 1 < ~bilebi00/_ARE/scripts/plot_profile_4Groups_4Motif.R > /dev/null
	fi
fi

date
echo DONE
