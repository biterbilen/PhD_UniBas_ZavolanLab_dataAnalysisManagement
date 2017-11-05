#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -j y
#$ -cwd
#$ -o log_file

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> log_file
} 
trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PERL5LIB=$HOME/lib:$PERL5LIB
hierClust=/import/bc2/home/nimwegen/GROUP/PhyloGibbs/Scripts/hier_clus_wm_simil

outdir=output_directory
tag=tag
infile=in_file
t2ccut=t2c_cut
libcut=lib_cut
ftag=ftag

mkdir -p $outdir/$tag; cd $outdir/$tag;

#Cluster motifs
less $ftag*/*wm > $ftag.all_wms;
rm -rf $ftag-run*;

sed -i 's/DB_/DB/g' $ftag.all_wms;

$hierClust $ftag.all_wms 4 0 0.2 $ftag.all_wms_report $ftag.allcls_wms > /dev/null; 
grep _ $ftag.allcls_wms; 
~bilebi00/PUM/bperlscripts/make_logo_from_hierClust_file.pl $ftag.allcls_wms .cls_wm > $ftag.cls_wms

~bilebi00/bsoft/weblogo/mylogo -f $ftag.cls_wms -F PNG -a -c -n -Y 2> /dev/null; 

echo "Create Highest Scoring Windows page"; 
~bilebi00/PUM/bperlscripts/order_IC_of_WM_windows.pl $ftag.cls_wms $ftag; 

~bilebi00/bsoft/weblogo/mylogo -f $ftag.cls_wms.best -F PNG -a -c -n -Y 2> /dev/null; 

