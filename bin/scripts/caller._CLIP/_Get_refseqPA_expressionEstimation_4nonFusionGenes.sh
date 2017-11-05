#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
# $ -q mpi8_very_long
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-4
# $ -o LOG._Get_refseqPA_expressionEstimation_4nonFusionGenes.sh.EWSR1.$TASK_ID
#$ -o LOG.EWSfuncsh_Get_refseqPA_expressionEstimation_4nonFusionGenes.sh.EWSR1.$TASK_ID
# $ -o LOG.EWSfuncsi1_Get_refseqPA_expressionEstimation_4nonFusionGenes.sh.EWSR1.$TASK_ID
# $ -o LOG.EWSfuncsi2_Get_refseqPA_expressionEstimation_4nonFusionGenes.sh.EWSR1.$TASK_ID
#TODO change queue type in bayseq mode

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

#TODO get data if not done before
#rsync -t mirz@web08:~bilebi00/BITER/UniqueGenomicAlignments/$id.bedplus.gz $indir/.
#exit;

project=DIS3L2; ids4GM=(212 1019); DEmethod=bayseq; cellLine=null; otag="";
project=EWSR1; ids4GM=(1344 1219); DEmethod=baySeq; cellLine=HeLa; otag="";
project=ARE; ids4GM=(1218 1219); DEmethod=DESeq; backgrounds=(siCTRL siUntreated); cellLine=null; otag="";
project=EWSfuncsh; ids4GM=(1037 1039); DEmethod=DESeq; backgrounds=(shCTRL); cellLine=null; otag="";
#project=EWSfuncsi1; ids4GM=(271 272); DEmethod=DESeq; backgrounds=(siCTRL); cellLine=null
#project=EWSfuncsi2; ids4GM=(651 647); DEmethod=DESeq; backgrounds=(siCTRL); cellLine=null
project=EWSR1; ids4GM=(); DEmethod=baySeq; cellLine=A673; otag=$cellLine;

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
indir=~bilebi00/_CLIP/data/clipz_UniqueGenomicAlignments
RefGeneFa=~bilebi00/DATA/hg19_ucsc_tracks/refGene.fa.gz
RefGene_3UTR_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_3UTR.gtf.gz
RefGene_exon_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_exon.gtf.gz
RefGene_gene_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_gene.gtf.gz
minReadLen=25

ids=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Expr" && (cl=="null" || cl==$5) { print $2;}' project=$project cl=$cellLine`)
annots=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Expr" && (cl=="null" || cl==$5)  { print $10;}' project=$project cl=$cellLine`)
annot=${annots[0]}

if [ $SGE_TASK_ID -gt ${#ids[@]} ]; then
	exit;
fi

if [ `echo ${annots[@]} | sed 's/ /\n/g' | sort | uniq | wc -l` -gt 1 ]; then
	echo inconsistent annotation categories for project $project ids:${ids[@]};
	exit;
fi

annotStr=" -wP $annot"
if [ $annot == ALL ]; then
	annotStr=" -wvP bacterial|fungus|vector|Markers_Adaptors|viral"
fi

id=${ids[$((SGE_TASK_ID-1))]}

echo Doing $id $indir for $annotStr

outdir=Analysis/Project_$project/RNAseq$otag

mkdir -p $outdir; pushd $outdir

#select representative transcripts according to the project name
if [ $SGE_TASK_ID == 1 ]; then
	#precedence for representative transcript selection is done by 
	# 1. NM_ or NR_ / NM_ if project=ARE
	# 2. RefSeq status(Reviewed or Validated or Inferred or Provisional or Predicted) 
	# 3. length of gene locus (takes longest)
	# 4. length of exonic region (takes longest) / length of 3'UTR region if project=ARE (takes longest)
	# for mRNA and ncRNA
	# omitting fusion genes (i.e. have a dash in the gene name)
	if [ ! -e representative.id.gid ]; then
		if [ $project == ARE ]; then
			zless $RefGene_3UTR_regions
		else 
			zless $RefGene_exon_regions
		fi | \
			awk -F "\t" 'BEGIN{OFS="\t"}{ print $1, $9, $5-$4+1}' | \
			sort -k 1,2 | bedtools groupby -g 1,2 -c 3 -o sum > length		

		~bilebi00/_PAPD5/scripts/outerJoin.pl $RefGene_gene_regions length 1,9 1,2 '1-9,12' NULL | \
			#corrected on 20121110 to include those non-fusion genes with a dash (-) in gene symbol
#			awk 'BEGIN{OFS="\t";}{ if ($6 < 6 && ($3 == "mRNA" || $3 == "ncRNA") && $10 !~ /-/) { k=1; if ($0 ~ /NR_/) { k=2; } print $10,$12,k,$6,($5-$4+1),$13} } ' | \
		awk 'BEGIN{OFS="\t";}
			{ 
				if (project == "ARE") {
					k=1;
					if ($6 < 6 && $3 == "mRNA" && ($10 ~ /-[1-9]/ || $10 !~ /-/)) {
						print $10,$12,k,$6,($5-$4+1),$13;
					} 
				} else {
					if ($6 < 6 && ($3 == "mRNA" || $3 == "ncRNA") && ($10 ~ /-[1-9]/ || $10 !~ /-/)) {
						k=1; if ($0 ~ /NR_/) { k=2; } print $10,$12,k,$6,($5-$4+1),$13;
					} 
				}
			}' project=$project | \
			sed 's/";*//g' | sort -k 1,1 -k 3,3g -k 4,4g -k 5,5gr -k 6,6gr > refgene.gid.id.bioType.status.locusLength.regionLength 
		less refgene.gid.id.bioType.status.locusLength.regionLength | \
			perl -e '$g=""; while(<>){ @t=split; print $t[1],"\t",$t[0],"\n" if ($g ne $t[0]); $g=$t[0];}' \
			> representative.id.gid
		~bilebi00/_xlinkEnrichment/scripts/grep_fa_given_id.pl representative.id.gid $RefGeneFa 0 0 1 1 | \
			awk '{ print ">"$1"\n"$2 }' > representative.fa

		rm length
	else
		echo representative.id.gid exists
	fi

	#get exonic regions
	if [ ! -e representative.bed.gz ]; then
		~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_exon_regions representative.id.gid 0 0 transcript_id | \
			~bilebi00/_CLIP/scripts/gtf2bed.pl | \
			sort -k 1,1 -k 2,2g | \
			gzip -c > representative.bed.gz
	else
		echo representative.bed.gz exists
	fi

	touch prep
fi

while [ ! -e prep ]; do
	echo waiting for prep
	sleep 1m
done

#select reads in mRNAseq based on minReadLen and annotStr
if [ ! -e $id.bed.gz ]; then
	zless $indir/$id.bedplus.gz | \
		awk '$NF > minReadLen{ print }' minReadLen=$minReadLen | \
		grep $annotStr | cut -f 1-6 | \
		sort -k 1,1 -k 2,2g | \
		gzip -c > $id.bed.gz
else
	echo $id.bed.gz exists
fi

#map mRNAseq reads onto representative genes exons
if [ ! -e $id.exp ]; then
	zless representative.bed.gz | \
		bedtools map -null 0 -s -o sum -a stdin -b $id.bed.gz | \
		awk 'BEGIN{OFS="\t"}{print $0,$3-$2}' | \
		sort -k 4,4 | bedtools groupby -g 4 -c 8,7 -o sum,sum \
		> $id.exp
else
	echo $id.exp exists
fi

touch run.$SGE_TASK_ID

if [ $SGE_TASK_ID == 1 ]; then
	#wait all runs to finish
	while [ `ls run* | wc -l` -lt ${#ids[@]} ]; do
		echo waiting for runs;
		sleep 1m
	done
	rm -rf prep run*

	# estimate expressed genes with selected samples by fitting a Gaussian Mixture model to the data
	soutdir=GM
	outf=RNAseq
	if [ ! -e $soutdir/$outf.w_geneSymbol ] && [ ${#ids4GM[@]} -gt 0 ]; then
		mkdir -p $soutdir; pushd $soutdir
		for id4GM in ${ids4GM[@]}; do
		  ln -sf ../$id4GM.exp .
		done
		files4GM=(`ls *exp`);
		ext=.mclust
		fcount=0;
		incr=3
		rand=$RANDOM
		for f in ${files4GM[@]}; do
			#get expressed genes first
			echo R --no-save --args $f $ext "<" ~bilebi00/_CLIP/scripts/classificationMixtureModel_4column.R  ">" /dev/null
			R --no-save --args $f $ext < ~bilebi00/_CLIP/scripts/classificationMixtureModel_4column.R  > /dev/null
			ff=$f$ext
			if [ "$f" == ${files4GM[0]} ]; then
				cp $ff $rand
				fields='1,2,3'
				i=3
				continue;
			fi
			~bilebi00/_PAPD5/scripts/outerJoin.pl $rand $ff 1,2 1,2 '' 0 > $rand.tmp;
			i=$((i + $incr));
			fields="$fields,$i"
			mv $rand.tmp $rand
			fcount=$((fcount + 1));
		done

		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=exp.classification.pdf *\.exp.classification.pdf
		rm -rf *\.exp.classification.pdf

		echo "id length ${files4GM[@]}" | sed 's/\s/\t/g' | sed 's/\.exp//g' > $outf
		awk 'NR>1{print;}' $rand | cut -f $fields >> $outf
		~bilebi00/_PAPD5/scripts/leftJoin.pl RNAseq ../representative.id.gid 1 1 '' | \
			awk 'BEGIN{OFS="\t"} { if ($NF == 0) { $NF="gid" } $(NF-1) = $NF; $NF=null; print }' > $outf.w_geneSymbol

		less $outf.w_geneSymbol | awk 'BEGIN{OFS="\t"}{ if (NR>1) { print $1,$NF}}' > ../expressed.id.gid
		rm $rand
		popd
	else
		echo $soutdir/$outf.w_geneSymbol exists
	fi
	#estimate differentially expressed genes with baySeq
	if [ $DEmethod == "baySeq" ]; then
		echo running bayseq...
		#TODO
		#merge counts
		#DIS3L2
		if [ $project == DIS3L2 ]; then
			ls *exp | while read i; do 
				id=`basename $i .exp`; 
				less $sid_prot_f | \
					#	awk -F "\t" '$2==id && $9==project && $3 !~ /HeLa/ && $3 !~ /oeDIS/ { print $2"\t"$3}' id=$id project=$project; 
				awk -F "\t" '$2==id && $9==project && $3 !~ /HeLa/ && $3 !~ /siDIS/ { print $2"\t"$3}' id=$id project=$project; 
			done
		elif [ $project == EWSR1 ] && [ $cellLine == A673 ]; then
			ls *exp | while read i; do 
				id=`basename $i .exp`; 
				less $sid_prot_f | \
					awk -F "\t" '$2==id && $9==project && (cl=="null" || $5==cl) { print $2"\t"$3}' id=$id project=$project cl=$cellLine; 
			done
		else
			ls *exp | while read i; do 
				id=`basename $i .exp`; 
				less $sid_prot_f | \
					awk -F "\t" '$2==id && $9==project && (cl=="null" || $5==cl) && $3 ~ /^si/ { print $2"\t"$3}' id=$id project=$project cl=$cellLine; 
			done
		fi | sort -k 2,2 > selected

		less selected | cut -f 1 | while read id; do
			idn=`less selected | grep -w $id | cut -f 2`;
			echo -n " $idn.id $idn.len $idn.exp"
		done | sed 's/ /\t/g' | sed 's/\t//' > input_4_bayseq

		echo >> input_4_bayseq
		files=`less selected | cut -f 1 | sed 's/$/.exp/' `
		paste $files >> input_4_bayseq 

		gidf=representative.id.gid
		nullData=T
		bootStraps=100
		ncl=8
		date
		R --no-save --args input_4_bayseq $gidf $nullData $bootStraps $ncl < ~bilebi00/_CLIP/scripts/baySeq_dave.R  &> log.bayseq
		date

		R --no-save --args input_4_bayseq_DE.txt .pernucExp < ~bilebi00/_CLIP/scripts/brandNewPrettyLookingScatter_4baySeq_wColumnSelection.R > /dev/null
	elif [ $DEmethod == DESeq ]; then
		gidf=expressed.id.gid
		#design matrix
		echo -e "id\tcondition" > selected
		for id in ${ids[@]}; do
			less $sid_prot_f | \
				awk -F "\t" '$2==id && $9==project { idn=$3; sub(/_[0-9]*$/,"",$3); print idn"\t"$2"\t"$3; }' id=$id project=$project
		done >> selected

		#data file
		files=`less selected | grep -v condition | cut -f 2 | sed 's/$/.exp/' `
		paste $files > input_4_deseq.tmp

		less selected | grep -v condition | cut -f 2 | while read id; do
			idn=`less selected | grep -w $id | cut -f 1`;
			echo -n " $idn.id $idn.len $idn.exp"
		done | sed 's/ /\t/g' | sed 's/\t//' > input_4_deseq
		echo >> input_4_deseq
		#select expressed
		nsamples=$((${#ids[@]} * 3))
		~bilebi00/_PAPD5/scripts/innerJoin.pl input_4_deseq.tmp $gidf 1 1 "'1-$nsamples'" >> input_4_deseq

		rm input_4_deseq.tmp

		padj=0.2
		for b in ${backgrounds[@]}; do
			R --no-save --args input_4_deseq selected $gidf $b $padj < ~bilebi00/_CLIP/scripts/DESeq_dave.R &> log.deseq.$b
		done

		if [ $project == ARE ]; then
			less selected | grep -vP "siCTRL|siUntreated|condition" | cut -f 2 | while read id; do
				idn=`less selected | grep -w $id | cut -f 1`;
				files=`ls input_4_deseq*_$idn.genenames`; 
				nms=`echo $files  | sed "s/.genenames//g" | sed "s/input_4_deseq_//g"`
				R --no-save --args `echo $files|sed 's/ /,/g'` "" venn.$idn $nms < ~bilebi00/_DIS3/scripts/venn.R > /dev/null
			done
		fi
	fi
fi

echo DONE
