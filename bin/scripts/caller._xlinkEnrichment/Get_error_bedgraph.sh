#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-40
#$ -o LOG.Get_error_bedgraph._tagindex.sh.$TASK_ID

#HOWTO for submission ad check runs
#s=Get_error_bedgraph.sh; for i in `seq 0 7`; do sed -e "s/_tagindex/$i/" $s > _$s$i; qsub _$s$i; echo $i; done
#tag='LOG.Get_error_bedgraph.4.sh.'; wc -l $tag* | sort -k 1,1gr; grep DONE $tag* | wc -l;

echo Does not work
echo Chris''s files are missing
exit;

SGE_TASK_ID=$1
#TODO set
db=hg18
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
sid_prot_f=/import/bc2/home/zavolan/bilebi00/_xlinkEnrichment/data/protein_sampleid_list
tags=(Ago2_MNase_Clip42 Ago2_MNase_Clip43 HuR_MNase_CLIP_42 HuR_MNase_Clip43 mRNASeq_4SU_365_XL_rep_A_Clip13 mRNASeq_4SU_365_XL_rep_B_Clip13 mRNASeq_4SU_No_XL_rep_A_Clip13 mRNASeq_4SU_No_XL_rep_B_Clip13)
tagi=_tagindex
annot=mRNA
annot=ALL
#------------------------
tag=${tags[$tagi]}
sid=`less $sid_prot_f | grep -w $tag | cut -f 2`;

export PATH=/import/bc2/home/zavolan/bilebi00/bin:$HOME/bin:$PATH

indir=~bilebi00/_xlinkEnrichment/data/Methods
outdir=Analysis/UniqueRawData_${annot}_woContamination/DB_all_$sid/
mkdir -p $outdir; pushd $outdir;


#extract unique ALL mappers
if [ $SGE_TASK_ID == 1 ]; then

	ln -sf $indir/$tag data

	#select unique mappers
	~bilebi00/_xlinkEnrichment/scripts/selectAnnot.pl data/mapped_sequences data/genome_mappings $annot 1 1 1 \
	 	2> mapped_sequences | sort -k 2,2 > genome_mappings.sorted
	less mapped_sequences | sort -k 1,1 > mapped_sequences.sorted

	if [ ! -e data/genome_mappings.DECODED ]; then
		echo "Generating DECODED"
		less genome_mappings.sorted  | awk 'BEGIN{OFS="\t";}{ print $3,$7-1,$8,$2,0,$6,$9,$1;}' > tmp
		bedtools nuc -s -seq "-fi" $genomefa -bed tmp > tmp1;
		less tmp1 | awk 'NR>1{ print; }' | sort -k 4,4  > genome_mappings.nucbed.sorted

		join -t $'\t' genome_mappings.nucbed.sorted mapped_sequences.sorted -1 4 -2 1 |  cut -f 1-8,18-20 > all
		~bilebi00/_xlinkEnrichment/scripts/decode.pl all | sort -k 2,2 > genome_mappings.sorted.DECODED

		~bilebi00/_xlinkEnrichment/scripts/selectAnnot.pl mapped_sequences.sorted genome_mappings.sorted.DECODED $annot 1 1 2> size > /dev/null

		rm -rf tmp tmp1 mapped_sequences all
	else
		~bilebi00/_xlinkEnrichment/scripts/selectAnnot.pl data/mapped_sequences data/genome_mappings.DECODED $annot 1 1 \
			2> size | sort -k 2,2 > genome_mappings.sorted.DECODED
	fi

	join -t $'\t' -1 1 -2 1 genome_mappings.sorted genome_mappings.sorted.DECODED > joined

	less joined | cut -f 17 | sort -r | uniq | grep -v -w "P" > types

	touch prepfinished
	echo `date` prepfinished

	~bilebi00/_xlinkEnrichment/scripts/density.pl joined copies
	for strand in "+" "-"; do
		less copies | perl ~bilebi00/www/MARA/scripts/extend_bed_fromScoreField.pl 1 | sort -k 1,1 > tmp
		bedtools genomecov -i tmp -g ~bilebi00/aux/human.$db.genome -bg -strand $strand > copies_$strand
		echo `date` bg for copies_$strand finished
	done
	rm -rf tmp;
	echo `date` copies finished
else
	while [ ! -e prepfinished ]; do
		echo Waiting for prepfinished.
		sleep 30s;
	done
fi

tc=`less types | wc -l`;

if [ $SGE_TASK_ID -gt $tc ]; then
	touch run.finished.$SGE_TASK_ID
	echo DONE
	exit;
fi

type=`less types | head -n $SGE_TASK_ID | tail -n 1`;

~bilebi00/_xlinkEnrichment/scripts/density.pl joined $type 2> $type.decodedtype

decodedtype=`cat $type.decodedtype`;
bedtools merge -d -1 -nms -s -scores sum -i $decodedtype > $decodedtype.mergeBed

for strand in "+" "-"; do
#	less $decodedtype.mergeBed | awk 'BEGIN{OFS="\t";}{ if ($6 == strand) { print $1,$2+1,$3,$5; }}' strand="$strand" > ${decodedtype}_$strand
	less $decodedtype.mergeBed | awk 'BEGIN{OFS="\t";}{ if ($6 == strand) { print $1,$2,$3,$5; }}' strand="$strand" > ${decodedtype}_$strand
done

rm -rf $type.decodedtype $decodedtype.mergeBed $decodedtype

touch run.finished.$SGE_TASK_ID

if [ $SGE_TASK_ID == 1 ]; then
	while [ `ls run.finished.* | wc -l` -lt 40 ]; do
		echo Waiting for 40 "run.finished.*".
		sleep 30s;
	done
	rm -rf prepfinished run.finished.* genome_mappings* copies mapped_sequences* joined
fi

echo `date` $decodedtype finished
echo DONE

