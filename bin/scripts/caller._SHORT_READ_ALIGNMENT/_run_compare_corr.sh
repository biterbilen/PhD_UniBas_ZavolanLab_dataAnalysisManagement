#!/bin/bash

outdir=Zavolan
outdir=Guo
outdir=Rep_chr19
pushd $outdir

#for inp in cufflinks_GTF_GMAP_mfpt8_*/transcripts.gtf; do 
for inp in Cufflinks_*/transcripts.gtf; do 
	echo $inp
	#outf=`perl -e '$_=shift; $_=~/cufflinks_GTF_GMAP_mfpt8_(\S+)\/transcripts.gtf/; print $1;' $inp`.exp
	#outf=`perl -e '$_=shift; $_=~/Cufflinks_(\S+)\/transcripts.gtf/; print $1;' $inp`.high_score.exp
	outf=`perl -e '$_=shift; $_=~/Cufflinks_(\S+)\/transcripts.gtf/; print $1;' $inp`.fields.exp
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/filter_gtf.pl 0 transcript 'transcript_id,gene_id,chr,score,FPKM' $inp | grep -v CUFF > $outf
	#less $inp | perl -e 'while(<>){ next if ($_=~/\texon\t/); $_=~/transcript_id "([\w\.\d]+)".*FPKM "([\d.]+)/; print "$1\t$2\n" if (defined $1 and defined $2); }' | grep -v CUFF > $outf;
done

rand=$RANDOM
out=transcripts.EXP
files=(`ls *.exp`);
cp ${files[0]} $rand.0; 
header="transcript_id gene_id chr ${files[0]}_score ${files[0]}";
ii='1,2,3,4,5';
n=5
for i in `seq 1 $((${#files[@]} - 1))`; do
	header="$header ${files[$i]}_score ${files[$i]}";
	j=$(($((i + 1)) * $n))
	ii=$ii",$((j - 1)),$j";
	echo "~bilebi00/_PAPD5/scripts/outerJoin.pl $rand.$(($i - 1)) ${files[$i]} 1 1 '' 0 > $rand.$i;"
	~bilebi00/_PAPD5/scripts/outerJoin.pl $rand.$(($i - 1)) ${files[$i]} 1 1 '' 0 > $rand.$i;
done

echo $header | sed 's/ /\t/g' > $out
cat $rand.$i | cut -f $ii >> $out
rm $rand*;

