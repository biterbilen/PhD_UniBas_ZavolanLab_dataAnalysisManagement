#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-16
#$ -o LOG._Get_xlinkEnrichment_summary.Nucleotide.sh.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_summary.clusters.sh.$TASK_ID

#XXX read
#each posterior estimate is done in a different core: 13 posteriors profiled -t parameter do it accordingly with protein count 26 for two proteins
export PATH=$HOME/bin:$PATH
#SGE_TASK_ID=2
#-----------------
#TODO set
#c=introns_refseq
#c=exons_refseq
#c=clusters
#c=Nucleotide
c=mutatedInCLIP;

annot=mRNA
extension=20 #window extension
db=hg18
#-----------------
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#-----------------
posteriors=(0.00001 0.0001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
posteriors=(0.00001 0.0001 0.001 0.01 0.02 0.03 0.04 0.05)
posteriors=(500 200 100 50 25 10 5 2 1)
jobcount=${#posteriors[@]}
sid_prot_f=~bilebi00/_xlinkEnrichment/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
#---------------------
i=`echo $((SGE_TASK_ID-1)) / ${#posteriors[@]} % ${#tags[@]} | bc`;
tag=${tags[$i]};

posteriori=`echo $((SGE_TASK_ID-1)) / 1 %  ${#posteriors[@]} | bc`;
posterior=${posteriors[$posteriori]};
#---------------------
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_${annot}_${c}/$tag/summary
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_${annot}_${c}/$tag/mutationCount/summary
mkdir -p $outdir; pushd $outdir;

#-----------------
id1=${ids[0]}
id2=${ids[1]}

#---------------------
ifs=(../$tag.$id1.$id1.T.xlink_mutratemax.bed ../$tag.$id2.$id2.T.xlink_mutratemax.bed)

#prepare buffer file with the most permissive posterior cut-off
if [ $posteriori == 0 ]; then
	for i in `seq 1 ${#ifs[@]}`; do
		ifile=${ifs[$(($i-1))]}
		ln -sf $ifile .
		less $ifile | awk '$5<posterior && $5 != "NA" {print}' posterior=${posteriors[$((${#posteriors[@]}-1))]} | sort -k 5,5g | gzip -c > buffer.$i.bed.gz
	done
	touch buffer.finished
else
	while [ ! -e buffer.finished ]; do
		echo waiting for buffer preperation
		sleep 1m;
	done
fi

echo Doing $SGE_TASK_ID $tag $posterior

#---------------------
#generate posterior profile
tmp=tmp.$posterior

percent=20
distance=0
for i in `seq 1 ${#ifs[@]}`; do
	zless buffer.$i.bed.gz | awk '$5<posterior{print}' posterior=$posterior | bedtools slop -b $distance -i stdin -g $genomelen > $tmp.slop.$i
	echo $tmp.slop.$i generated
done
#get intersecting ones within a distance
b=`bedtools closest -s -d -t "first" -a $tmp.slop.1 -b $tmp.slop.2 | awk '$13==0{print}' | wc -l`
s1=`less $tmp.slop.1 | wc -l`;
s2=`less $tmp.slop.2 | wc -l`;
u=`cat $tmp.slop.1 $tmp.slop.2 | sort | uniq | wc -l`;
echo $tag $distance $posterior $b $s1 $s2 $u > $tmp.posterior.profile

touch $tmp.posteriorProfile.finished

if [ $posteriori == 0 ]; then
	while [ "`ls tmp*posteriorProfile.finished | wc -l`" -lt $jobcount ]; do
		ls tmp*posterior.profile | wc -l
		echo waiting for posterior.profile files
		sleep 1m;
	done
	#merge
	#mark posterior cutoff candidates as the reproduced position percent drops below $percent
	#get high scoring positions and extend
	echo "proteinname distance posterior common sample1 sample2 union common_over_union_percent" | sed -e 's/ /\t/g' > $tag.posterior.profile
	for posterior in ${posteriors[@]}; do
		cat tmp.$posterior.posterior.profile
	done |sed -e 's/ /\t/g' | awk '{k=100*$4/$7; if (k<percent) { print"*"$0"\t"k; } else { print $0"\t"k; } }' percent=$percent >> $tag.posterior.profile
	#clean
	rm -rf tmp* *finished

	#get min of marked posteriors as the cutoff; then plot and generate ofiles for further use
	posterior=`less $tag.posterior.profile | grep "^*" | bedtools groupby -g 1 -c 3 -o min -i stdin | cut -f 2`
	R --no-save --args $tag.posterior.profile $posterior < ~bilebi00/_xlinkEnrichment/scripts/plot_posterior_profile.R > /dev/null
	ofs=(sample1_XL.postcut.sorted.bed sample2_XL.postcut.sorted.bed)
	for i in `seq 1 ${#ofs[@]}`; do
		ofile=${ofs[$(($i-1))]}
		#get high scoring positions
		zless buffer.$i.bed.gz | awk '$5<posterior{print}' posterior=$posterior > $ofile
		#get non-overlapping positions
		less $ofile | bedtools cluster -s -d $((extension-1)) -i stdin | sort -k 1,1 -k 7,7 -k 5,5g | \
			bedtools groupby -i stdin -g 1,7 -c 5 -o max -full | cut -f 1-6 | sort -k 5,5g > repr_$ofile
	done

	#---------------------
	#get overlap in top bins to make comparison with Kishore et al. 2011 
	echo "proteinname sample1name sample2name distance topN union intersection intersection_over_topN_percent" | sed -e 's/ /\t/g' > $tag.overlap.stats
	for type in "" repr_; do
		id1_XL=$type${ofs[0]}
		id2_XL=$type${ofs[1]}
		for distance in 1 5 20; do
			for topN in 1000 2000 3000 4000 5000 10000 15000; do
				#extend regions
				less $id1_XL | head -n $topN > tmp.1
				less $id2_XL | head -n $topN > tmp.2
				#get intersecting ones within a distance
				b=`bedtools closest -s -d -t "first" -a tmp.1 -b tmp.2 | awk '$13 > -distance && $13 < distance {print}' distance=$distance | wc -l`
				u=`cat tmp.1 tmp.2 | sort | uniq | wc -l`;
				echo $tag $id1_XL $id2_XL $distance $topN $u $b
			done
		done
	done | sed -e 's/ /\t/g' | awk '{k=100*$7/$5; print $0"\t"k; }' >> $tag.overlap.stats
	R --no-save --args $tag.overlap.stats < ~bilebi00/_xlinkEnrichment/scripts/positional_reproducibility.R > /dev/null
	rm -rf tmp*
fi
#---------------------

echo DONE
exit;




