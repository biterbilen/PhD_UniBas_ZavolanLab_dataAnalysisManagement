#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-100
#$ -o LOG._Get_PhastCons_4nonOverlappingInterval.sh$TASK_ID

export PATH=$HOME/bin:$PATH

outdir=Analysis/regions/
mkdir -p $outdir; pushd $outdir;

#TODO set
gtff=hg19_GMAP_GENEintrons_refseq.gtf.gz; #takes more than an hour
op=median
windowsize=1000
op=mean
windowsize=1

strands=(- +);
chrs=(`zless $gtff | cut -f 1 | sort | uniq`)
tag=`basename $gtff .gtf.gz`;

if [ $SGE_TASK_ID -gt ${#chrs[@]} ]; then
	touch $tag.tmp.runfinished.$SGE_TASK_ID
	echo exiting;
	exit;
fi

chr=${chrs[$((SGE_TASK_ID-1))]}
ofile=$tag.tmp_$chr
rm -rf $ofile.cons

for strand in "+" "-"; do
	echo Doing $chr $strand
	zless $gtff | awk '$1==chr && $7==strand {print }' chr=$chr strand=$strand > $ofile
	~bilebi00/_CLIP/scripts/parse_wigFile_4gtf_w_nonOverlappingIntervals.pl $ofile $op $windowsize >> $ofile.cons;
done

touch $tag.tmp.runfinished.$SGE_TASK_ID

if [ $SGE_TASK_ID == 1 ]; then
	#wait
	while [ "`ls $tag.tmp.runfinished* | wc -l`" -lt 100 ]; do
		echo waiting for other jobs
		sleep 1m
	done
	#merge results
	cat $tag*.cons | gzip -c > ${tag}_wPhastCons_op${op}_ws${windowsize}.gtf.gz
	#clean
	rm -rf $tag.tmp*;
fi

echo DONE
