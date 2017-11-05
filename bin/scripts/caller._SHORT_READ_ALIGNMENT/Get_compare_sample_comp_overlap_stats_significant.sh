#!/bin/bash

function run {
tag=$1
c1=$2
c2=$3
t1=$4
t2=$5
s1=$6
s2=$7
tr1=$8
tr2=$9
two=${10}
two2=${11}
~bilebi00/_PAPD5/scripts/innerJoin.pl ${t1}RNAseq-$c1-0${two}VS${t1}RNAseq-$c1-${tr1}Dis3L2$two.gene_exp.diff.$s1 ${t2}RNAseq-$c2-0${two2}VS${t2}RNAseq-$c2-${tr2}Dis3L2$two2.gene_exp.diff.$s2 1 1 3 | sort | uniq > $tag-$c1-$c2-$t1-$t2-$s1-$s2-$tr1-$tr2-$two-$two2; 
echo "#" $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 rep2:$two rep2:$two2: `less ${t1}RNAseq-$c1-0${two}VS${t1}RNAseq-$c1-${tr1}Dis3L2$two.gene_exp.diff.$s1 | cut -f 3 | sort | uniq | wc -l` `less ${t2}RNAseq-$c2-0${two2}VS${t2}RNAseq-$c2-${tr2}Dis3L2$two2.gene_exp.diff.$s2 | cut -f 3 | sort | uniq | wc -l` `less $tag-$c1-$c2-$t1-$t2-$s1-$s2-$tr1-$tr2-$two-$two2 | wc -l`;

return 0	
}

#TODO set
skipped=""; #"HEK.siDis3L2"; #""
#trusted=_cuffdiff_trusted_superClusters_gene;
trusted=_cuffdiff_combined_superClusters_gene;
outdir=Stepanka/OVERLAP_STATS_significant${trusted}_skipped$skipped; 

#trusted="_trusted_superClusters"; #"_combined"; #"_trusted";
#outdir=Shivendra/OVERLAP_STATS$trusted; 

mkdir -p $outdir; pushd $outdir;

rm -rf sigCounts;
for d in ../CLIP_comparison$trusted/*0VS*diff ../CLIP_comparison$trusted/*0_2VS*diff; do 
	if [ "$skipped" != "" ] && [ "`echo $d | grep \"$skipped\"`" != "" ]; then continue; fi
	i=`basename $d`; 
	less $d | awk '$8>$9 && $14~/yes/{ print; }' > $i.pos; 
	less $d | awk '$8<$9 && $14~/yes/{ print; }' > $i.neg;  
	echo "# $i `less $d | awk '$14~/yes/{ print; }' | wc -l` `less $i.pos | wc -l` `less $i.neg | wc -l`" >> sigCounts; 
done

# Stabilized genes
#appear in the same direction for overexpression and siRNA knockdown
#TODO set
stabfiles=`ls *oeDis3L2*\.neg *siDis3L2*\.pos`;
cat $stabfiles | cut -f 3 | sort | uniq -c | sort -k 1,1gr | awk 'OFS="\t"{ print $1,$2; }'  > all_stabilized
less all_stabilized  | awk '$1>1{print $2;}' > all_stabilized_trusted
less all_stabilized | awk '{ print $2; }' | sed 's/,/\n/g' | sort | uniq > stabilized

#-----------------------
# Destabilized genes
#appear in the same direction for overexpression and siRNA knockdown
#TODO set
destabfiles=`ls *oeDis3L2*\.pos *siDis3L2*\.neg`;
cat $destabfiles | cut -f 3 | sort | uniq -c | sort -k 1,1gr | awk 'OFS="\t"{ print $1,$2; }' > all_destabilized
less all_destabilized  | awk '$1>1{print $2;}' > all_destabilized_trusted
less all_destabilized | awk '{ print $2; }' | sed 's/,/\n/g' | sort | uniq > destabilized

#-----------------------
~bilebi00/_PAPD5/scripts/innerJoin.pl all_stabilized all_destabilized 2 2 '' > both_stab_unstab
~bilebi00/_PAPD5/scripts/innerJoin.pl all_stabilized_trusted all_destabilized_trusted 1 1 1 > both_stab_unstab_trusted

#-----------------------
outfiles=(stabilized.matrix destabilized.matrix)

for outfile in ${outfiles[@]}; do
	echo $outfile;
	files="$destabfiles";
	if [ $outfile == "stabilized.matrix" ]; then
		files=`ls $stabfiles`;
	fi

	echo -n " " > $outfile
	echo $files >> $outfile
	for f1 in $files; do
		echo -n $f1;
		for f2 in $files; do
			echo -n " "`~bilebi00/_PAPD5/scripts/innerJoin.pl $f1 $f2 1 1 1 | wc -l`;
		done
		echo
	done >> $outfile
	sed -i 's/.gene.diff//g' $outfile

	~bilebi00/bin/R --no-save --args $outfile scaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null
	~bilebi00/bin/R --no-save --args $outfile unscaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null
done
#-----------------------
outfile=distrusted.matrix

echo -n " " > $outfile
echo $stabfiles >> $outfile
for f1 in $destabfiles; do
	echo -n $f1;
	for f2 in $stabfiles; do
		echo -n " "`~bilebi00/_PAPD5/scripts/innerJoin.pl $f1 $f2 1 1 1 | wc -l`;
	done
	echo
done >> $outfile
sed -i 's/.gene_exp.diff//g' $outfile

~bilebi00/bin/R --no-save --args $outfile scaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null
~bilebi00/bin/R --no-save --args $outfile unscaled $N < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/intersectionmatrix.R > /dev/null

exit;

#get maximum clip score
for sig_gene_file in all_destabilized all_stabilized; do
	echo $sig_gene_file;
	for clipf in ~bilebi00/_DIS3/4WEB_DIS3L2/PC*\.data; do
		outf=`basename $clipf .data`.$sig_gene_file;
		echo $outf.$sig_gene_file;
		less $sig_gene_file | while read i; do if [ $i == "-" ]; then continue; fi; less $clipf | grep -w $i | grep -w EXON | sort -k 7,7gr | head -n 1 ; done > $outf
	done
done

header="chr source name beg end score score2 str frm foreground foreground2 background background2 annot_type annot_ids annot_id"

echo $header | sed 's/ /\t/g' > PC.all_stabilized; 
cat PC_DIS3L2_*all_stabilized >> PC.all_stabilized

echo $header | sed 's/ /\t/g' > PC.all_destabilized; 
cat PC_DIS3L2_*all_destabilized >> PC.all_destabilized

#run a.R in the same folder

R --no-save --args PC.all_stabilized  < a.R
R --no-save --args PC.all_destabilized  < a.R
