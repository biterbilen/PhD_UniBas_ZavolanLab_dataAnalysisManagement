#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_seed_counts_and_cuffdiffDE.sh 

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

PATH=$HOME/bin:$PATH
export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

#outdir=Guo; pushd $outdir;
#outdir=Zavolan; pushd $outdir;
outdir=Scheiffele; pushd $outdir;

#This is not used
#~bilebi00/_SHORT_READ_ALIGNMENT/scripts/gtf2clsformat.pl ~bilebi00/_SHORT_READ_ALIGNMENT/Guo/merged_asm/merged.gtf "seqname,start,end,strand,gene_name,oId,nearest_ref,exon_number,class_code" > merged_asm/merged.cls

#Fields
#3xN+10+14
nsamples=`ls | grep "^Cufflinks_" | wc -l`;
bi=$(($((nsamples * 3)) + 10 + 5));
ei=$((bi + 9));
for cat in genes isoforms; do 
	echo $cat
	catv=`perl -e '$a=shift; chop $a; print $a;' $cat`;
	~bilebi00/_PAPD5/scripts/innerJoin.pl Cuffdiff/$cat.fpkm_tracking Cuffdiff/${catv}_exp.diff 1 1 "2,3,5,10,$bi-$ei" > $cat.all  
	#DE is ok
	less $cat.all | perl -e 'my $cat=shift; while(<>) { my ($cl,$id,$gid,$exs,$g1,$g2,$des,$fpkm1,$fpkm2,$lr,$ts,$pv,$qv,$sig)=split; if ($exs eq "OK" and $des eq "OK" and (($cl eq "=" and $cat eq "isoforms") or ($cat eq "genes" and $gid ne "-"))) { print join("\t",$id,$gid,$g1,$g2,$fpkm1,$fpkm2,$lr,$ts,$pv,$qv), "\n";} }' $cat > $cat.selected_based_on_cuffdiff_stat;
	#EXP estimation is OK
	less $cat.all | perl -e 'my $cat=shift; while(<>) { my ($cl,$id,$gid,$exs,$g1,$g2,$des,$fpkm1,$fpkm2,$lr,$ts,$pv,$qv,$sig)=split; if ($exs eq "OK" and (($cl eq "=" and $cat eq "isoforms") or ($cat eq "genes" and $gid ne "-"))) { print join("\t",$id,$gid,$g1,$g2,$fpkm1,$fpkm2,$lr,$ts,$pv,$qv), "\n";} }' $cat > $cat.selected_based_on_cufflinks_stat;
done

exit;
