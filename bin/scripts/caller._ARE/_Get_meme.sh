#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-8
#$ -o LOG._Get_meme.sh$TASK_ID

#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.mutrate.pdf */all.*.mutrate.pdf
#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.freqstats.pdf */all.*.freqstats.pdf

export PATH=$HOME/bin:/import/bc2/soft/bin/hmmer:$PATH

#SGE_TASK_ID=4
#-----------------
#TODO set
#c=clusters; f="5 / 100";
c=Nucleotide; n=20; N=2000; f="1 /100"; #top 1 percent
#b=bottom
b=neighbour
#--
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#-----------------
sid_prot_f=~bilebi00/_ARE/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w PARCLIP | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w PARCLIP | grep -w $tag | cut -f 2 | sort`)
idsXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep XL | cut -f 2 | sort`)
idswoXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep -v XL | cut -f 2 | sort`)

#-----------------
outdir=Analysis/xlinkEnrichment/MEME_${c}_${b}/$tag
mkdir -p $outdir; pushd $outdir;

#-----------------
nuc=T
id1=${ids[0]}
id2=${ids[1]}

libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
mu_ref_XL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
mu_ref_woXL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
merged_XL=../../basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.*.$nuc.xlink_mutrate$mu_ref_XL.bed

#T=`less $merged_XL | wc -l` 
#N=`echo "$T * $f" | bc`;

if [ $c == "clusters" ]; then
	less $merged_XL | head -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > foreground.nucbed
	if [ $b == "bottom" ]; then
		less $merged_XL | tail -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background.nucbed
	elif [ $b == "neighbour" ]; then
		cluster_len=`less $merged_XL | awk 'BEGIN{OFS="\t";} {print tag,$3-$2}' tag=$tag | bedtools groupby -g 1 -c 2 -o mean | awk '{ print int($2)+1;}'`;
		less $merged_XL | bedtools flank -b $cluster_len -i stdin -g $genomelen | bedtools intersect -v -s -a stdin -b $merged_XL -f 0.1 | head -n $N | \
			bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background.nucbed
	fi
elif [ $c == "Nucleotide" ]; then
	cluster_len=$((n*2+1))
	less $merged_XL | bedtools slop -b $n -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$cluster_len > all.CCRs
	less all.CCRs | head -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > foreground.nucbed
	if [ $b == "bottom" ]; then
		less all.CCRs |	tail -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background.nucbed
	elif [ $b == "neighbour" ]; then
		less all.CCRs | bedtools flank -b $cluster_len -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$cluster_len | \
			bedtools intersect -v -s -a stdin -b all.CCRs -f 0.1 | head -n $N | \
			bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background.nucbed
	fi
	rm all.CCRs
fi

less foreground.nucbed | head -n $((N+1)) | awk 'BEGIN{c=0;} { if (NR>1){ c=c+1; print ">seq"c"\n"$16; }}' > foreground.fa
less background.nucbed | head -n $((N+1)) | awk 'BEGIN{c=0;} { if (NR>1){ c=c+1; print ">seq"c"\n"$16; }}' > background.fa

nmotifs=1
minw=4
maxw=8
mod=anr
mod=zoops

meme=~bilebi00/bin/meme

minsites=`echo "$N * 95 / 100" | bc`;
maxsites=$N;
if [ $mod == "anr" ]; then maxsites=$((N * 2)); fi

m=2
less background.fa | fasta-get-markov -m $m -norc > markov$m

echo -e "\n" $meme foreground.fa -maxsize 1000000 -sf $tag -dna -bfile markov2 -mod $mod -nmotifs $nmotifs -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc ./ -nostatus
$meme foreground.fa -maxsize 1000000 -sf $tag -dna -bfile markov$m -mod $mod -nmotifs $nmotifs -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc ./ -nostatus > /dev/null 
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$tag.logo.pdf logo1.eps 

echo DONE
exit;

#XXX Get motif distance to the crosslink site
ll | grep "^d" | awk '{ print $9; }' | grep -v "^_" | grep -v -P "\.|EWSR1|CFIm68" | while read tag; do 
	less $tag/meme.txt | \
		perl -e '$tag=shift; $start=0; while(<>){ if ($_ =~ /^Sequence name\s*Start/) { $_=<>; $start=1; next; } if ($_=~/^\-/ && $start==1) {exit; } if ($start) { @t=split/\W+/;print $t[1]."\t".$tag."\n";} }' $tag | \
		sort | bedtools groupby -g 1 -c 1 -full -o count -i stdin | sort -k 1,1g; \
	done | awk 'BEGIN{OFS="\t";} { $1=$1-21;print $0; }' | sed -e 's/hnRNPQ/hnRNPQ1/g' > motif_distance
R --no-save --args motif_distance < ~bilebi00/_ARE/scripts/plot_profile_4Groups_4Motif.R > /dev/null

