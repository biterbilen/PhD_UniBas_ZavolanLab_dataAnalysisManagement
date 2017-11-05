#!/bin/bash

#XXX set
task=2
tag=mRNAseq2

maind=`pwd`;
outdir=Analysis/clipz_RNAseq
mkdir -p $outdir; pushd $outdir

if [ $task == 1 ]; then
	#XXX set
	ids=(`less $maind/data/protein_sampleid_list | grep $tag | cut -f 2`)
	ps=(`less $maind/data/protein_sampleid_list | grep $tag | cut -f 3`)

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
	tag=HeLa
	tag=HEK

	gif=~bilebi00/_EWSR1/data/hg19_TR.info.all
	files=(`ls *$tag*exp`);
	incr=4
	outf=${tag}_RNAseq.raw
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

	echo "id ${files[@]}" | sed 's/\s/\t/g' | sed 's/\.exp//g' > $outf
	awk 'NR>1{print;}' $rand | cut -f $fields >> $outf
	rm $rand*

	R --no-save --args $outf < ~bilebi00/_KNOCKDOWNS/scripts/qn.R > /dev/null
	~bilebi00/_DIS3/scripts/assign_gene_id.pl $outf $gif 0 NA 1 > $outf.w_geneSymbol
	~bilebi00/_DIS3/scripts/assign_gene_id.pl $outf.qn $gif 0 NA 1 > $outf.qn.w_geneSymbol

	#head -n 1 $outf.qn.w_geneSymbol > dima_clip_targets_in_Hela 
	#for i in DIS3L2 DIS3 DIS3L HSPA1B ACTG1 PMAIP1 MAZ BIRC5 RNF44 SCD EEF2; do grep -w $i $outf.qn.w_geneSymbol; done >> dima_clip_targets_in_Hela

	~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf 1-$((fcount+1)) exp > $outf.df
	~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf.qn 1-$((fcount+1)) exp > $outf.qn.df

	R --no-save --args $outf $outf.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null
	R --no-save --args $outf.qn $outf.qn.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null

fi
