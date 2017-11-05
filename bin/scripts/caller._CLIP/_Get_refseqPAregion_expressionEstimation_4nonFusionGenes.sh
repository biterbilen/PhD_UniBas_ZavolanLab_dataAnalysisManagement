#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q mpi8_very_long
# $ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-10
#$ -o LOG._Get_refseqPAregion_expressionEstimation_4nonFusionGenes.sh.EWSR1.$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1";
}
trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=$HOME/bin:$PATH
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
project=ARE
project=DIS3L2
project=EWSR1

region=intron
region=exon
ids=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Expr" { print $2;}' project=$project`)

if [ $SGE_TASK_ID -gt ${#ids[@]} ]; then
	exit;
fi

id=${ids[$((SGE_TASK_ID-1))]}

echo Doing $id $indir for $annotStr

outdir=Analysis/Project_$project/RNAseq/$region
mkdir -p $outdir; pushd $outdir

#preps coordinates of for the $region of expressed genes estimated by Bayseq for the siCTRL and siEWSR1 samples 
if [ $SGE_TASK_ID == 1 ]; then
	if [ ! -e ${region}s.bed.gz ]; then
		is=`less ../input_4_bayseq_DE.txt | awk 'NR==1{ print "1,"NF;}'`
		~/_PAPD5/scripts/leftJoin.pl ../input_4_bayseq_DE.txt ../representative.bed.gz 1 4 $is | \
			awk 'NR>1{ print }' | sort | uniq > expressed.id.gid
		~/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl ~/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_$region.gtf.gz expressed.id.gid 0 0 transcript_id | \
			~/_CLIP/scripts/gtf2bed.pl transcript_id | \
			sort -k 1,1 -k 2,2g | \
			gzip -c > ${region}s.bed.gz
	fi
	touch prep
fi

while [ ! -e prep ]; do
	echo waiting for prep
	sleep 1m
done

if [ ! -e $id.exp ]; then
	zless ${region}s.bed.gz | \
		bedtools map -null 0 -s -o sum -a stdin -b ../$id.bed.gz | \
		awk 'BEGIN{OFS="\t"}{print $4"|"$1":"($2+1)"-"$3$6,$3-$2,$7}' \
		> $id.exp
else
	echo $id.exp exists
fi

touch run.$SGE_TASK_ID

if [ $SGE_TASK_ID == 1 ]; then
	while [ `ls run* | wc -l` -lt ${#ids[@]} ]; do
		echo waiting for runs;
		sleep 1m
	done
	rm -rf prep run*

	echo running bayseq...
	#TODO
	#merge counts
	#DIS3L2
	if [ $project == DIS3L2 ]; then
		ls *exp | while read i; do 
			id=`basename $i .exp`; less $sid_prot_f | \
	#	awk -F "\t" '$2==id && $9==project && $3 !~ /HeLa/ && $3 !~ /oeDIS/ { print $2"\t"$3}' id=$id project=$project; 
		awk -F "\t" '$2==id && $9==project && $3 !~ /HeLa/ && $3 !~ /siDIS/ { print $2"\t"$3}' id=$id project=$project; 
		done
	else
		ls *exp | while read i; do 
			id=`basename $i .exp`; less $sid_prot_f | \
			awk -F "\t" '$2==id && $9==project && $3 ~ /^si/ { print $2"\t"$3}' id=$id project=$project; 
		done
	fi | sort -k 2,2 > selected

	less selected | cut -f 2 | while read id; do
		echo -n " $id.id $id.len $id.exp"
	done | sed 's/ /\t/g' | sed 's/\t//' > input_4_bayseq
	echo >> input_4_bayseq
	files=`less selected | cut -f 1 | sed 's/$/.exp/' `
	paste $files >> input_4_bayseq 

	#estimate DE
	gidf=null
	nullData=T
	bootStraps=100;
	bootStraps=1;
	ncl=8
	date
	R --no-save --args input_4_bayseq $gidf $nullData $bootStraps $ncl < ~/_CLIP/scripts/baySeq_dave.R  &> log.bayseq
	date

	#add gid -pufff
	less input_4_bayseq_DE.txt | sed "s/ids/ids|tag/" | awk -F "|" '{ print $1"\t"$2}' > a
	n=`less a | awk 'NR==1{ print NF}'`
	m=$((n+2))
	~/_PAPD5/scripts/leftJoin.pl a expressed.id.gid 1 1 "1-$n,$m" NULL | \
		awk -F "\t" 'BEGIN{OFS="\t"}{if(NR==1){$NF="gid"; } else { $1=$1"|"$2;} print }' | \
		cut -f 1,3- > input_4_bayseq_DE.txt
	rm a

	#plot splom with the given likelihoodcut for down and upreguated
	likelihoodcut=0.95
	likelihoodcut=0.5
	R --no-save --args input_4_bayseq_DE.txt .pernucExp $likelihoodcut < ~/_CLIP/scripts/brandNewPrettyLookingScatter_4baySeq_wColumnSelection.R > /dev/null

fi

echo DONE
