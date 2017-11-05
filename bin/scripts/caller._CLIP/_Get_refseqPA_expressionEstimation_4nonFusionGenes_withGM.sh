#TODO set
#run after _Get_refseqPA_expressionEstimation_4nonFusionGenes.sh

#Deprecated; done in _Get_refseqPA_expressionEstimation_4nonFusionGenes.sh

exit

ids=(212 1019)
project=DIS3L2

ids=(1344 1219)
project=EWSR1

ids=(1218 1219)
project=ARE


outdir=Analysis/Project_$project/RNAseq/GM
mkdir -p $outdir; pushd $outdir

for id in ${ids[@]}; do
	ln -sf ../$id.exp .
done

files=(`ls *exp`);

outf=RNAseq

ext=.mclust
fcount=0;
incr=3
rand=$RANDOM
for f in ${files[@]}; do
	#get expressed genes first
	echo R --no-save --args $f $ext "<" ~/_CLIP/scripts/classificationMixtureModel_4column.R  ">" /dev/null
	R --no-save --args $f $ext < ~/_CLIP/scripts/classificationMixtureModel_4column.R  > /dev/null
	ff=$f$ext
	if [ "$f" == ${files[0]} ]; then
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

echo "id length ${files[@]}" | sed 's/\s/\t/g' | sed 's/\.exp//g' > $outf
awk 'NR>1{print;}' $rand | cut -f $fields >> $outf
~/_PAPD5/scripts/leftJoin.pl RNAseq ../representative.id.gid 1 1 '' | \
	awk 'BEGIN{OFS="\t"} { if ($NF == 0) { $NF="gid" } $(NF-1) = $NF; $NF=null; print }' > $outf.w_geneSymbol

less $outf.w_geneSymbol  | awk '{ print $NF}' | grep -v gid | sort | uniq -c | awk '{ print $2"\t"$1}' > expressed.gid

rm $rand
