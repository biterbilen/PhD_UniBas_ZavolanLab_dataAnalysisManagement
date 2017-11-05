#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -t 1-runtimes
#$ -o log_file_$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> log_file_$SGE_TASK_ID
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

outdir=output_directory
samples=samples_file
guide_gtff=guide_gtff
indir_tag="indir_tag"

mkdir -p $outdir/4GB_$indir_tag;
pushd $outdir/4GB_$indir_tag;

ft=`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`

if [ "` echo $ft | grep $indir_tag `" == "" ]; then exit; fi #the links should not exist here

ln -sf ../TopHat_$ft/accepted_hits.bam $ft.bam;
ln -sf ../TopHat_$ft/accepted_hits.bam.bai $ft.bam.bai;
ln -sf ../Cufflinks_$ft/transcripts.gtf $ft.trx.gtf

sleep $SGE_TASK_ID;
ln -sf ../Cuffcompare_$indir_tag/Cuffcompare_$indir_tag.combined.gtf .
ln -sf ../Cuffcompare_$indir_tag/Cuffcompare_$indir_tag.trusted.gtf .
ln -sf $guide_gtff .
