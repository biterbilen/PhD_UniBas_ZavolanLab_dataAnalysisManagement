#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_very_long@@qc_nehalem
#$ -P project_zavolan
# $ -l h_vmem=60G
#$ -j y
#$ -cwd
#$ -o LOG._run_topHat.chr.rep_A

PATH=$HOME/bin:$PATH

bti_dir=~/DATA/Bowtie_DB_INDEX
#db_name=hg18_chr19 #hg_TR #TODO set here
db_name=hg18 #hg_TR #TODO set here
tag=A
outdir=rep_$tag.$db_name

inp=inp.fa
thread_num=8

mkdir -p $outdir; cd $outdir;

misosrcdir=~bilebi00/_DOWNLOAD/_FOR_TRANSCRIPT_REGULATION/MixtureIsoforms
#Convert TopHat sam file to sorted bam file
#1.
#mkdir -p $misodir/sam-output
#samtools view -Sbh $tophatdir/accepted_hits.sam > $misodir/sam-output/accepted_hits.bam
#samtools sort $misodir/sam-output/accepted_hits.bam $misodir/sam-output/accepted_hits.sorted
#samtools index $misodir/sam-output/accepted_hits.bam $misodir/sam-output/accepted_hits.sorted.bam
#2.
python $misosrcdir/sam_to_bam.py --convert tophat/accepted_hits.sam miso/sam-output --ref $bti_dir/$db_name.fa.fai


## pickle alternative events from cufflinks
##TODO gtf2gff conversion
exit;
python $misosrcdir/index_gff.py --index cufflinks/transcripts.gtf miso/indexed
## Run MISO on a pair of paired-end sample (with insert length distribution with mean 250,
## standard deviation 15) using the mouse genome skipped exon annotations using the
## the cluster

# Compute Psi values for control sample
python run_events_analysis.py --compute-events-psi control control.counts --event-type SE --read-len 35 --overhang-len 4 --output-dir control
python $misosrcdir/run_events_analysis.py --compute-genes-psi mm9/pickled/SE data/control.bam --output-dir SE/control/ --read-len 37 --paired-end 250 15 --use-cluster

# Compute Psi values for knockdown sample
python run_events_analysis.py --compute-genes-psi mm9/pickled/SE data/knockdown.bam --output-dir SE/knockdown/ --read-len 35 --paired-end 250 15 --use-cluster


## Summarize the output (only run this once --compute-genes-psi finished!)
## This will create a "summary" directory in SE/control/ and in SE/knockdown/
python run_miso.py --summarize-samples SE/control/ SE/control/
python run_miso.py --summarize-samples SE/knockdown/ SE/knockdown/

## Detect differentially expressed isoforms between "control" and "knockdown"
## This will compute Bayes factors and delta Psi values between the samples
## and place the results in the directory SE/comparisons/control_vs_knockdown
python run_miso.py --compare-samples SE/control/ SE/knockdown/ SE/comparisons/

less ~bilebi00/_SHORT_READ_ALIGNMENT/data/repA_unmapped.seq |perl -e ' $_=<>;  while(<>) { @t=split; print ">sequd_$t[0]\t$t[1]\n"; }' | sort -u | perl -e 'while(<>){@t=split; print $t[0], "\n", $t[1], "\n"; }' >> $inp

