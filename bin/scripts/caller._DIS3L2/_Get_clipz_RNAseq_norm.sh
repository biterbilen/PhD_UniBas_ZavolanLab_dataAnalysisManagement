#!/bin/bash

#XXX set
task=2
type=expressed_quantile; #expressed transcripts and quantile normalization
type=common_scale; #scale normalization

maind=`pwd`;
outdir=Analysis/clipz_RNAseq
mkdir -p $outdir; pushd $outdir

if [ $task == 1 ]; then
	tag=mRNAseq
	tag2=Dima

	#XXX set
	ids=(`less $maind/data/protein_sampleid_list | grep $tag | grep $tag2 | cut -f 2`)
	ps=(`less $maind/data/protein_sampleid_list | grep $tag | grep $tag2 | cut -f 3`)

	for SGE_TASK_ID in `seq 1 ${#ids[@]}`; do
		p=${ps[$((SGE_TASK_ID-1))]}
		if [ -e $p.exp ]; then
			continue;
		fi
		id=${ids[$((SGE_TASK_ID-1))]}
		rsync mirz@bc2-web08.bc2.unibas.ch:~mirz/clipzServer/data/samples/$id/transcript_expression $p.exp;
		echo $id
	done	
fi

if [ $task == 2 ]; then
	#XXX set
	tag=""
	outf=${tag}_RNAseq_$type.raw
	gif=~bilebi00/_EWSR1/data/hg19_TR.info.all

	if [ $type == common_scale ]; then
		cat *exp | cut -f 1 | grep -v "^#" | sort | uniq -c | awk '$1==n{print $2;}' n=`ls *exp | wc -l` > common_ids 
		#join files                                                                                       
		header="id";
		is="2";
		k=1
		for f in *exp; do
			k=$((k+4)); 
			is="$is,$k";
			header="$header $f";  
			~bilebi00/_PAPD5/scripts/leftJoin.pl common_ids $f 1 1 '' > $outf
			mv $outf common_ids
		done      
		echo -e "$header" | sed 's/.exp//g' | sed 's/ /\t/g' > $outf;
		less common_ids | cut -f $is >> $outf
		rm common_ids

		~bilebi00/_DIS3/scripts/assign_gene_id.pl $outf $gif 0 NA 1 > $outf.w_geneSymbol
		~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf 1-`ls *exp | wc -l` exp > $outf.df

		rm -rf $outf.splom.pdf $outf.heatmap.pdf $outf.density.pdf $outf.bw.pdf
		R --no-save --args $outf $outf.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$outf.pdf $outf\.*\.pdf
		rm -rf $outf.splom.pdf $outf.heatmap.pdf $outf.density.pdf $outf.bw.pdf
	elif [ $type == expressed_quantile ]; then
		files=(`ls *$tag*exp`);
		incr=4
		rand=$RANDOM
		echo ${files[@]}

		#export HOME=/import/bc2/home/zavolan/

		ext=.mclust
		fcount=0;
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
			~bilebi00/_PAPD5/scripts/outerJoin.pl $rand $ff 1 1 '' 0 > $rand.tmp;
			i=$((i + $incr));
			fields="$fields,$i"
			mv $rand.tmp $rand
			fcount=$((fcount + 1));
		done 	
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=exp.classification.pdf *\.exp.classification.pdf
		rm -rf *\.exp.classification.pdf

		echo "id ${files[@]}" | sed 's/\s/\t/g' | sed 's/\.exp//g' > $outf
		awk 'NR>1{print;}' $rand | cut -f $fields >> $outf
		rm $rand*

		R --no-save --args $outf < ~bilebi00/_KNOCKDOWNS/scripts/qn.R > /dev/null
		~bilebi00/_DIS3/scripts/assign_gene_id.pl $outf $gif 0 NA 1 > $outf.w_geneSymbol
		~bilebi00/_DIS3/scripts/assign_gene_id.pl $outf.qn $gif 0 NA 1 > $outf.qn.w_geneSymbol

		~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf 1-$((fcount+1)) exp > $outf.df
		~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf.qn 1-$((fcount+1)) exp > $outf.qn.df

		#plot and merge plots
		rm -rf $outf.splom.pdf $outf.heatmap.pdf $outf.density.pdf $outf.bw.pdf
		R --no-save --args $outf $outf.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$outf.pdf $outf\.*\.pdf
		rm -rf $outf.splom.pdf $outf.heatmap.pdf $outf.density.pdf $outf.bw.pdf

		rm -rf $outf.qn\.*\.pdf
		R --no-save --args $outf.qn $outf.qn.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$outf.qn.pdf $outf.qn\.*\.pdf
		rm -rf $outf.qn\.*\.pdf

		rm -rf *\.df $outf $outf.qn
	fi
	echo DONE
fi
