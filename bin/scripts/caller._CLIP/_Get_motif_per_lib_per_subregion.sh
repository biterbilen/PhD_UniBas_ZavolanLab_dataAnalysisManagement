#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 4-4
#$ -o LOG._Get_motif_per_lib_per_subregion.sh$TASK_ID

#script checks kmer freqs of top $N regions of libraries of $type annotated, if specified, extended $b nucleotides and uses the $b region for the lower bound of the observed motif proportion; uses mu_refmax 

#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.mutrate.pdf */all.*.mutrate.pdf
#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.freqstats.pdf */all.*.freqstats.pdf

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
#SGE_TASK_ID=4
#-----------------
#TODO set
subregion=intron; distance=1000; restrictN=112;
subregion=exon; distance=100; restrictN=516;
restrictN=smallest
restrictN=1000;
kmers=(1 2 3 4 5 6)
posterior=-10
c=xlinkEnrichment
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
n=20; 
b=bottom
b=neighbour
#-----------------
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk '$1==tag { print $2; }' tag=$tag | sort -r`)
region_len=$((n*2+1))

#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_posterior${posterior}_$c/motif_perlib_background${b}_n${n}_$subregion
mkdir -p $outdir; pushd $outdir;

extrafields=0
#-----------------DATA PREP
for id in ${idsf[@]}; do
	cg=../$subregion/common_subregion_genes
	ifile=../$subregion/$id.bed.w_common_subregion_genes
	bfile=../all.$id.bed.w_type_genename.gz
	allCCRs=common_subregion_genes.all$id.CCRs.bed.gz
	extrafields=$((`less $ifile | awk 'NR==1{ print NF;}'`-6))

	if [ ! -e foreground$id.nucbed.gz ]; then
		zless $ifile | bedtools slop -b $n -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$region_len | \
			bedtools nuc -s -seq "-fi" $genomefa -bed stdin | gzip -c > foreground$id.nucbed.gz
	else
		echo foreground exists
	fi
	N=$((`zless foreground$id.nucbed.gz | wc -l`-1));
	if [ ! -e $allCCRs ]; then
		~bilebi00/_xlinkEnrichment/scripts/grep_f_w_like_w_column.pl $cg $bfile 7 | grep -w $subregion | \
			bedtools slop -b $n -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$region_len | gzip -c > $allCCRs
		if [ $b == "bottom" ]; then
			zless $allCCRs | tail -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin | \
				gzip -c > background$id.nucbed.gz
		elif [ $b == "neighbour" ]; then
			zless $allCCRs | bedtools flank -b $region_len -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$region_len | \
				bedtools intersect -v -s -a stdin -b $allCCRs -f 0.05 | head -n $N | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin | gzip -c > background$id.nucbed.gz
		fi
	else
		echo background exists
	fi
done
echo data prep finished

popd
#new folder
boutdir=`basename $outdir`;
Noutdir=${outdir}_N$restrictN
mkdir -p $Noutdir; pushd $Noutdir;

#-----------------DATA PREP IF THE OPTION for `restrictN` IS SET
#echo <<'END'
if [ "$restrictN" == smallest ]; then
	#find the minimum sequence having file and its id
	minfile="";
	minid="";
	N=100000;
	#search for the minimum sequence having file
	for f in `ls ../$boutdir/fore*gz`; do
		fN=`zless $f | wc -l`;
		of=`basename $f .gz`;
		if [ $fN -lt $N ]; then N=$fN; minfile=$f; minid=`echo $of | sed 's/foreground//' | sed 's/.nucbed//'`; fi
	done

	#prepare foreground files by extracting `N` sequences from the previous foreground file that is closest to the sites of the `minfile` sites and has the highest score; 
	#softlinks the minfile foregorund file here
	for f in `ls ../$boutdir/fore*gz`; do
		of=`basename $f .gz`;
		if [ -e $of.gz ]; then echo $restrictN files exists; break; fi
		if [ "`echo $of | grep foreground$minid.nucbed`" != "" ]; then 
			ln -sf $f .; continue; 
		fi
		NF=`zless $f | awk 'NR==1{print NF;}'`
		NF2=$((NF*2+1))
		zless $f | head -n 1 > $of
		bedtools closest -s -d -t first -a $f -b $minfile | awk '$NF<distance{ print; }' distance=$distance | sort -k $NF2,${NF2}g -k 5,5g | \
			head -n $((N-1)) | cut -f 1-$NF >> $of
		gzip $of
	done
	#prepare background files by extracting `N` sequences from the previous background file that is closest to the sites of the `minfile` sites and has the highest score; 
	#softlinks the minfile background file here
	for f in `ls ../$boutdir/back*gz`; do
		of=`basename $f .gz`;
		if [ -e $of.gz ]; then echo $restrictN files exists; break; fi
		if [ "`echo $of | grep background$minid.nucbed`" != "" ]; then 
			ln -sf $f .; continue; 
		fi
		NF=`zless $f | awk 'NR==1{print NF;}'`
		NF2=$((NF*2+1))
		#		ffile=`echo $of | sed -e 's/background/foreground/'`.gz
		zless $f | head -n 1 > $of
		bedtools closest -s -d -t first -a $f -b $minfile | awk '$NF<distance{ print; }' distance=$(($distance*5)) | sort -k $NF2,${NF2}g -k 5,5g | \
			head -n $((N-1)) | cut -f 1-$NF >> $of
		gzip $of
	done
elif [[ $restrictN =~ "^[0-9]+$" ]]; then
	idsf=()
	for f in `ls ../$boutdir/fore*gz`; do
		if [ `zless $f | wc -l` -ge $(($restrictN+1)) ]; then
			of=`basename $f`;
			zless $f | head -n $(($restrictN+1)) | gzip -c > $of; 
			fb=`echo $f | sed -e 's/foreground/background/'`
			ofb=`echo $of | sed -e 's/foreground/background/'`
			zless $fb | head -n $(($restrictN+1)) | gzip -c > $ofb; 
			#add to idsf array
			id=`echo $of | sed -e 's/foreground//' | sed -e 's/.nucbed.gz//'`
			idsf=(${idsf[@]} $id)
		else
			echo "skipping $f:#positions <$((restrictN+1))"
			continue;
		fi
	done
else
	echo Error: unknown restrictN=$restrictN!
	exit;
fi
#END

#-----------------MOTIF SEARCH
for id in ${idsf[@]}; do
	#-----------------KMER
	mkdir -p kmer; pushd kmer
	for kmer in ${kmers[@]}; do
		~bilebi00/_EWSR1/scripts/count_kmers.pl ../foreground$id.nucbed.gz $kmer bed count $tag ${kmer}mer $extrafields > ${kmer}mers$id 2> ${kmer}mers$id.onames
		~bilebi00/_EWSR1/scripts/count_kmers.pl ../background$id.nucbed.gz $kmer bed count $tag ${kmer}mer $extrafields > ${kmer}mers$id.background 2> /dev/null
		R --no-save --args ${kmer}mers$id ${kmer}mers$id.onames ${kmer}mers$id.background < ~bilebi00/_EWSR1/scripts/kmer_enrichment_w_background_mu_ref.R > /dev/null
		R --no-save --args ${kmer}mers$id ${kmer}mers$id.onames < ~bilebi00/_EWSR1/scripts/kmer_enrichment_w_background_mu_ref.R > /dev/null
		for m in ${kmer}mers$id.*pbeta; do
			less $m | sort -k 5,5g -k 2,2gr > ${m}_
			mv ${m}_ $m;
		done
	done
	echo kmers ${kmers[@]} finished
	popd
	#-----------------MEME
	mkdir -p meme/$id; pushd meme/$id	
	#prepinput
	zless ../../foreground$id.nucbed.gz | awk 'BEGIN{c=0;} { if (NR>1){ c=c+1; print ">seq"c"\n"$(16+extrafields); }}' extrafields=$extrafields > foreground.fa
	zless ../../background$id.nucbed.gz | awk 'BEGIN{c=0;} { if (NR>1){ c=c+1; print ">seq"c"\n"$(16+extrafields); }}' extrafields=$extrafields > background.fa

	N=`zless ../../foreground$id.nucbed.gz | wc -l`;
	nmotifs=1
	minw=6
	maxw=6
	mod=anr
	mod=zoops

	minsites=`echo "$N * 95 / 100" | bc`;
	maxsites=$N;

	#tried p=50:10:80 and the motif GAAGAA or GG[AT]GG[AT] for subregion=exon
	if [ $tag == EWSR1 ]; then mod=anr; minsites=`echo "$N * 80 / 100" | bc`; fi
	#	if [ $tag == EWSR1 ]; then minsites=`echo "$N * 80 / 100" | bc`; fi
	if [ $mod == "anr" ]; then maxsites=$((N * 2)); fi

	m=$(($minw-4))
	less background.fa | fasta-get-markov -m $m -norc > markov$m

	echo -e "\n" meme foreground.fa -p 8 -maxsize 1000000 -sf $tag -dna -bfile markov2 -mod $mod -nmotifs $nmotifs -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc ./ -nostatus
	meme foreground.fa -p 8 -maxsize 1000000 -sf $tag -dna -bfile markov$m -mod $mod -nmotifs $nmotifs -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc ./ -nostatus > /dev/null
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$tag.logo.pdf logo1.eps
	popd
done

pushd kmer
topn=8
for mfix in backgroundFile_mu_ref uniform_mu_ref; do
	for kmer in ${kmers[@]}; do
		idns=()
		files=()
		for i in `seq 1 ${#idsf[@]}`; do
			id1=${idsf[$((i-1))]};
			id1n="`less $sid_prot_f | awk -F "\t" '$2==id{print $6;}' id=$id1`"
			idns=("${idns[@]}" "'$id1n'")
			files=(${files[@]} ${kmer}mers$id1.$mfix.pbeta)
			for j in `seq $((i+1)) ${#idsf[@]}`; do
				id2=${idsf[$((j-1))]};
				id2n="`less $sid_prot_f | awk -F "\t" '$2==id{print $6;}' id=$id2`"
				R --no-save --args ${kmer}mers$id1.$mfix.pbeta ${kmer}mers$id2.$mfix.pbeta ${kmer}mer.$id1.$id2.$mfix $topn "$id1n" "$id2n" < ~bilebi00/_DIS3L2/scripts/kmer_rep.R > /dev/null
			done
		done
		R --no-save --args logpbeta all.${kmer}mer.$mfix ${#files[@]} ${files[@]} "${idns[@]}" < ~bilebi00/_CLIP/scripts/kmer_dengrogram.R > /dev/null
	done
done
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="kmerheatmaps_`basename $outdir`.pdf" all*mer*.pdf;
rm -rf all*.heatmap.pdf
popd

pushd meme
ls -latr | grep "^d" | awk '$9 !~ /^\./ { print $9; }' | while read tag; do
idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$tag`";
less $tag/meme.txt | \
	perl -e '$tag=shift; $start=0; while(<>){ if ($_ =~ /^Sequence name\s*Start/) { $_=<>; $start=1; next; } if ($_=~/^\-/ && $start==1) {exit; } if ($start) { @t=split/\W+/;print $t[1]."\t".$tag."\n";} }' "$idn" | \
	sort | bedtools groupby -g 1 -c 1 -full -o count -i stdin | sort -k 1,1g | \
	awk -F "\t" 'BEGIN{OFS="\t";}{ if (c==0){c=$1+1;} for(i=c;i<$1;i++) { print c,$2,0; c=c+1; } print $0; c=$1+1; }'; \
done | awk -F "\t" 'BEGIN{OFS="\t";} { $1=$1-region_len; print $0; }' region_len=$((n+1)) > motif_distance
rows=2
if [ ${#idsf[@]} -lt 3 ]; then
	rows=1
fi
R --no-save --args motif_distance $rows < ~bilebi00/_ARE/scripts/plot_profile_4Groups_4Motif.R > /dev/null
popd

echo DONE
