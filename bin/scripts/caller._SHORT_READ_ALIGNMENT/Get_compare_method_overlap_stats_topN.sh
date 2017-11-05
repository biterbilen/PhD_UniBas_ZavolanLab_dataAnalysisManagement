#!/bin/bash

export PERL5LIB=$HOME/lib:$PERL5LIB


#TODO set
N=100;
src=_superClusters
skipped=""; #HEK.siDis3L2
moutdir=Stepanka; #Scheiffele; #Shivendra
gi=_gene
#TODO set
trusted1=_cuffdiff_combined
trusted2=_cuffdiff_trusted
trusted1=_clipz_RNAseq_raw
trusted2=_cuffdiff_trusted
trusted1=_clipz_RNAseq_qn
trusted2=_cuffdiff_trusted
#trusted1=_clipz_RNAseq_raw
#trusted2=_cuffdiff_combined
#trusted1=_clipz_RNAseq_qn
#trusted2=_cuffdiff_combined

outdir=$moutdir/OVERLAP_STATS_merged_top${N}${trusted1}${trusted2}${gi}_skipped$skipped; 
mkdir -p $outdir; pushd $outdir;

indir1=../OVERLAP_STATS_top${N}${trusted1}${src}${gi}_skipped$skipped; 
indir2=../OVERLAP_STATS_top${N}${trusted2}${src}${gi}_skipped$skipped; 

# Stabilized genes
#appear in the same direction for overexpression and siRNA knockdown
#TODO set
stabfiles1=`ls $indir1/*oeDis3L2*\.neg $indir1/*siDis3L2*\.pos`;
stabfiles2=`ls $indir2/*oeDis3L2*\.neg $indir2/*siDis3L2*\.pos`;

# Destabilized genes
#appear in the same direction for overexpression and siRNA knockdown
#TODO set
destabfiles1=`ls $indir1/*oeDis3L2*\.pos $indir1/*siDis3L2*\.neg`;
destabfiles2=`ls $indir2/*oeDis3L2*\.pos $indir2/*siDis3L2*\.neg`;

#-----------------------
outfiles=(stabilized.matrix destabilized.matrix)

for outfile in ${outfiles[@]}; do
	echo $outfile;
	files1="$destabfiles1";
	files2="$destabfiles2";
	if [ $outfile == "stabilized.matrix" ]; then
		files1="$stabfiles1"
		files2="$stabfiles2"
	fi

	echo -n " " > $outfile
	echo $files1 >> $outfile
	for f1 in $files2; do
		echo -n $f1;
		for f2 in $files1; do
			echo -n " "`~bilebi00/_PAPD5/scripts/innerJoin.pl $f1 $f2 3 3 3 | sort | uniq | wc -l`;	
		done
		echo
	done >> $outfile
	sed -i 's/.gene.diff//g' $outfile
	sed -i 's/\.\.\/'"`basename $indir1`"'/'$trusted1'/g' $outfile
	sed -i 's/\.\.\/'"`basename $indir2`"'/'$trusted2'/g' $outfile
	~bilebi00/bin/R --no-save --args $outfile scaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null
	~bilebi00/bin/R --no-save --args $outfile unscaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null
done

exit;
#-----------------------

outfile=distrusted.matrix

echo -n " " > $outfile
for f1 in $destabfiles; do
	echo -n $f1;
	for f2 in $stabfiles; do
		echo -n " "`~bilebi00/_PAPD5/scripts/innerJoin.pl $f1 $f2 1 1 1 | wc -l`;	
	done
	echo
done >> $outfile
sed -i 's/.gene_exp.diff//g' $outfile
sed -i 's/\.\.\/'"`basename $indir1`"'/'$trusted1'/g' $outfile
sed -i 's/\.\.\/'"`basename $indir2`"'/'$trusted2'/g' $outfile

~bilebi00/bin/R --no-save --args $outfile scaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null
~bilebi00/bin/R --no-save --args $outfile unscaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null

exit;
