#!/bin/bash

#TODO set
#1.
project=DIS3L2
expr=GM
#exprJoin is expr=GM specific parameter to get the HeLa and HEK common genes
exprJoin=inner
exprJoin=outer
#referenceList=../GeneExpression_GMExp_scaleNorminner_Untreated/RNAseq.w_geneSymbol
referenceList=../GeneExpression_GMExp_scaleNormouter_oeDIS3L2_cellHEK293/RNAseq.w_geneSymbol
#expr=common #GM : Gaussian mixture

#2.
norm=quantile
norm=scale
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list

cell=""
cell=HeLa
cell=HEK293

#tag=siCTRL
#tag=Untreated
#tag=Rep
tag=oeDIS3L2
tag=si
tag=1si
#outdir=Analysis/clipz_RNAseq/Project_$project/GeneExpression_${expr}Exp_${norm}Norm${exprJoin}_$tag
outdir=Analysis/Project_$project/clipz_RNAseq/GeneExpression_${expr}Exp_${norm}Norm${exprJoin}_$tag

if [ "$cell" != "" ]; then
	outdir=${outdir}_cell$cell
fi
mkdir -p $outdir; pushd $outdir

ids=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ /Expr/ && $0 ~ tag && (cell=="" || $5==cell){ print $2}' project=$project cell="$cell" tag="$tag"`)

for id in ${ids[@]}; do
	idn="`less $sid_prot_f | awk -F "\t" '$2==id && $9==project{print $3; }' project=$project id=$id`";
	if [ ! -e ../$idn.exp ]; then
		echo Getting $id $idn
		rsync -t mirz@bc2-web08.bc2.unibas.ch:~mirz/clipzServer/data/samples/$id/transcript_expression ../$idn.exp;
	fi
	ln -sf ../$idn.exp 
done	

files=(`ls *exp`);

outf=RNAseq
gif=~bilebi00/_EWSR1/data/hg19_TR.info.all

if [ $expr == common ]; then
	cat ${files[@]} | cut -f 1 | grep -v "^#" | sort | uniq -c | awk '$1==n{print $2;}' n=${#files[@]} > common_ids 
	#join files                                                                                       
	header="id";
	is="2";
	k=1
	for f in ${files[@]}; do
		k=$((k+4)); 
		is="$is,$k";
		header="$header `basename $f`";  
		~bilebi00/_PAPD5/scripts/leftJoin.pl common_ids $f 1 1 '' > $outf
		mv $outf common_ids
	done      
	echo -e "$header" | sed 's/.exp//g' | sed 's/ /\t/g' > $outf;
	less common_ids | cut -f $is >> $outf
	rm common_ids
elif [ $expr == GM ]; then 
	ext=.mclust
	fcount=0;
	incr=4
	rand=$RANDOM
	for f in ${files[@]}; do 
		#get expressed genes first
		~bilebi00/bin/R --no-save --args $f $ext < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/classificationMixtureModel.R > /dev/null
		ff=$f$ext
		if [ "$f" == ${files[0]} ]; then
			cp $ff $rand
			fields='1,4'
			i=4
			continue;
		fi
		echo $f
		if [ $exprJoin == "inner" ]; then
			~bilebi00/_PAPD5/scripts/innerJoin.pl $rand $ff 1 1 '' > $rand.tmp;
		elif [ $exprJoin == "outer" ]; then
			~bilebi00/_PAPD5/scripts/outerJoin.pl $rand $ff 1 1 '' 0 > $rand.tmp;
		else 
			echo $exprJoin not defined
			exit
		fi
		i=$((i + $incr));
		fields="$fields,$i"
		mv $rand.tmp $rand
		fcount=$((fcount + 1));
	done 

	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=exp.classification.pdf *\.exp.classification.pdf
	rm -rf *\.exp.classification.pdf

	echo "id ${files[@]}" | sed 's/\s/\t/g' | sed 's/\.exp//g' > $outf
	awk 'NR>1{print;}' $rand | cut -f $fields >> $outf

	if ! echo $referenceList | grep --quiet "^$" ;then
		cc=`less $referenceList | awk -F "\t" '{print NF+1; exit; }'`
		~bilebi00/_PAPD5/scripts/leftJoin.pl $referenceList $outf 1 1 '' 0 | cut -f $cc- > $rand.tmp;
		mv $rand.tmp $outf
	fi	
	rm $rand*
fi

if [ $norm == quantile ]; then
	echo quantile normalization
	R --no-save --args $outf < ~bilebi00/_KNOCKDOWNS/scripts/qn.R > /dev/null
	mv $outf.qn > $outf
fi

~bilebi00/_DIS3/scripts/assign_gene_id.pl $outf $gif 0 NA 1 > $outf.w_geneSymbol
~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf 1-${#files[@]} exp > $outf.df

rm -rf $outf.splom.pdf $outf.heatmap.pdf $outf.density.pdf $outf.bw.pdf
R --no-save --args $outf $outf.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$outf.pdf $outf\.*\.pdf
rm -rf $outf.splom.pdf $outf.heatmap.pdf $outf.density.pdf $outf.bw.pdf

zless $outf.w_geneSymbol | awk 'NR==1 || ($1 ~ /NM_/ && $NF !~ /\-/) { print }' > $outf.refseq.mRNA.single

echo DONE

