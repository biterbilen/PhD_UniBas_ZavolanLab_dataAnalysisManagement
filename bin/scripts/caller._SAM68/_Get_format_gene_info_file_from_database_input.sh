#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_format_gene_info_file_from_database_input

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> LOG._Get_format_gene_info_file_from_database_input
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV


dbf=~bilebi00/_SAM68/data/tab_for_biter_mm9
mcoordsf=~bilebi00/_SAM68/data/mm_auxData/mRNA_CDS.COORDS
allannot_f=~bilebi00/_SAM68/data/mm_auxData/all_annotations.fa

out=~bilebi00/_SAM68/data/mm_TR.info.all

~bilebi00/CDS_CLIP_ElMMo/scripts/generate_gene_info_file.pl $dbf $mcoordsf $allannot_f > $out 2> $out.stat
