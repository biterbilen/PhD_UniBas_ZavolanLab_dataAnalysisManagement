outdir=Analysis/Project_DIS3L2/nonCodingRNAExpEst/
mkdir -p $outdir; pushd $outdir

id=1598
lens=(0 22)
types=(tRNA rRNA snoRNA snRNA miRNA)
for type in ${types[@]}; do
	echo -n Doing $type
	zless ../annot_regions/hg19_UCSC_$type.gtf.gz | sort -k 1,1 -k 4,4g | gzip -c > hg19_UCSC_$type.gtf.gz
	for len in ${lens[@]}; do
		echo -e "\t" $id $len
		less ~/_CLIP/data/clipz_UniqueGenomicAlignments/$id.bedplus.gz | \
			awk '$9==type && $3-$2>len{ print }' type=$type len=$len | cut -f 1-6 | \
			sort -k 1,1 -k 2,2g | \
			bedtools map -s -o sum -null 0 -b stdin -a hg19_UCSC_$type.gtf.gz | \
			awk -F "\t" 'BEGIN{OFS="\t"}{ $6 = $10; print }' | cut -f 1-9 | \
			gzip -c > minlen$len.hg19_UCSC_$type.gtf.gz
		if [ $len == 100 ]; then
			zless minlen$len.hg19_UCSC_$type.gtf.gz | awk 'BEGIN{OFS="\t"}{ if ($6>0) { print $10, $5-$4+1, $6}}' | \
				sed 's/";*//g' > $type.exp; 
			R --no-save --args $type.exp .mclust < ~/_CLIP/scripts/classificationMixtureModel_4column.R > /dev/null
		fi
	done
	rm hg19_UCSC_$type.gtf.gz
done

popd
