#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
# $ -q fs_long@@qc_nehalem
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
# $ -l h_vmem=60G
#$ -j y
#$ -cwd
# $ -o LOG._run_new_Cufflinks.chr19.rep_A_upperquannorm_$TASK_ID
# $ -o LOG._run_new_Cufflinks.chr19.rep_A_bias_$TASK_ID
# $ -o LOG._run_new_Cufflinks.chr19.rep_A_mask_$TASK_ID
# $ -o LOG._run_new_Cufflinks.chr19.rep_A_mask2_$TASK_ID
# $ -o LOG._run_new_Cufflinks.chr19.rep_A_final_$TASK_ID
#$ -o LOG._run_new_Cufflinks.chr19.rep_A_final2_$TASK_ID
#$ -t 4-4
# $ -t 1-5

PATH=$HOME/bin:$PATH

bti_dir=~/DATA/Bowtie_DB_INDEX
db_name=hg18_chr19 #hg_TR #TODO set here
tag=A
outdir=rep_$tag.$db_name

thread_num=8

mkdir -p $outdir; cd $outdir;

#-N TRUE doesn't work

mfpt=(2 4 6 8 10)
i=$(($SGE_TASK_ID - 1));
m=${mfpt[$i]};

#gtff=~bilebi00/_SHORT_READ_ALIGNMENT/data/refseq_genes.gtf
gtff=~bilebi00/_SHORT_READ_ALIGNMENT/data/GMAP_EXONS.gtf
maskfile=maskfile

~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rich' ~/DATA/hg18_ucsc_tracks/rmsk.gtf.gz > $maskfile
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rRNA,Mitochondrial' ~/DATA/hg18_ucsc_tracks/rRNAGene.gtf.gz >> $maskfile

#maskf=~bilebi00/_SHORT_READ_ALIGNMENT/data/mask_example.gtf
#maskf2=~bilebi00/_SHORT_READ_ALIGNMENT/data/mask_example2.gtf

#guide GTF
#mask file -M mask.gtf 

#seq bias correction -b genome.fa OK
#multi-read correct -u OK
#upper quantile norm -N OK
#-j --pre-mrna-fraction <0.0-1.0>
#Warning: Using default Gaussian distribution due to insufficient paired-end reads in open ranges.  It is recommended that correct paramaters (--frag-len-mean and --frag-len-std-dev) be provided.

#tag=guide;cmd="cufflinks -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"
#echo $cmd > $tag.log; nohup $cmd 2>> $tag.log &

#tag=multireadcorrect;cmd="cufflinks -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"
#echo $cmd > $tag.log; nohup $cmd 2>> $tag.log &

#tag=upperquannorm;cmd="cufflinks -N -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"
#echo $cmd > $tag.log; nohup $cmd 2>> $tag.log &

#tag=bias;cmd="cufflinks -b $bti_dir/$db_name.fa -N -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"
#echo $cmd > $tag.log; nohup $cmd 2>> $tag.log &

#tag=mask;cmd="cufflinks --no-update-check -q -M $maskf  -b $bti_dir/$db_name.fa -N -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"
#echo $cmd > $tag.log; nohup $cmd 2>> $tag.log &

#tag=mask2;cmd="cufflinks --no-update-check -q -M $maskf2  -b $bti_dir/$db_name.fa -N -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"

#tag=final;cmd="cufflinks --no-update-check -q -M $maskfile  -b $bti_dir/$db_name.fa -N -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"

tag=final2;cmd="cufflinks -s 36 -s 1 --no-update-check -q -M $maskfile  -b $bti_dir/$db_name.fa -N -u -p $thread_num -o cufflinks_GTF_GMAP_mfpt${m}_$tag -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 tophat/accepted_hits.bam"
echo $cmd
$cmd

