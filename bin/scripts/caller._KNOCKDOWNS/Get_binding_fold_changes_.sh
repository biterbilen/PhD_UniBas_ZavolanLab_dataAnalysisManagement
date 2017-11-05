#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -cwd
#$ -j y
#$ -o log_file
#$ -t 1-5

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

export PATH=~bilebi00/bin:/import/bc2/soft/bin:$PATH

outdir=output_directory
qnTable=qn_table

#TODO delete me later
qnTable=~bilebi00/_KNOCKDOWNS/Normalization/293wt-RNAseq_HEK293_Nanofectin_RNAseq_HEK293_siAUF1_RNAseq_HEK293_siGFP_RNAseq_HEK293_siTIA1_RNAseq_HEK293_sihnRNPC_RNAseq_mRNASeq_No_4SU_No_XL_rep_A_Clip13_mRNASeq_No_4SU_No_XL_rep_B_Clip13_si_GFP_mRNASEQ_si_HuR_mRNASEQ.qn

outdir=FoldChange
outdir=FoldChangeControls
#outdir=FoldChangeExpTrx
#outdir=FoldChangeControlsExpTrx

mkdir -p $outdir; cd $outdir

#4th run
#Controls vs Controls
#20110121
fgroups=(si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq 293wt.RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13 mRNASeq_No_4SU_No_XL_rep_B_Clip13);
bgroups="HEK293_siGFP_RNAseq";
fgroups=(HEK293_siGFP_RNAseq si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq 293wt.RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13);
bgroups="mRNASeq_No_4SU_No_XL_rep_B_Clip13";
fgroups=(mRNASeq_No_4SU_No_XL_rep_B_Clip13 HEK293_siGFP_RNAseq si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13);
bgroups="293wt.RNAseq";

#f=${fgroups[$(($SGE_TASK_ID-1))]};
#echo Getting Fold changes for $f;
#rm -rf $f*stats*;
#for b in $bgroups; do
#	echo Doing for background group $b;
#	awk '{ print $2; }'  ~bilebi00/_KNOCKDOWNS/number_of_2_times_clip_sites | while read i; do ct=`perl -e '$s=shift; $s=~/AREBP\/(.+)\//; print $1;' $i`; echo $i $ct; ~bilebi00/_KNOCKDOWNS/scripts/getFoldChange_forBinnedIdMatch.pl ~bilebi00/_KNOCKDOWNS/$i $ct $f $b $qnTable 1; done
#done
#echo Plotting $f
#~bilebi00/_KNOCKDOWNS/scripts/run_call_wrapper_barweb.sh /import/bc2/soft/app/matlab/current/Linux/  $f
#convert $f*.pdf all_$f.pdf
#exit;

#3rd run
#Controls vs KDs expressed transcripts
#20110122
fgroups=(si_HuR_mRNASEQ HEK293_sihnRNPC_RNAseq HEK293_siTIA1_RNAseq HEK293_siAUF1_RNAseq);
bgroups="si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq 293wt.RNAseq HEK293_siGFP_RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13 mRNASeq_No_4SU_No_XL_rep_B_Clip13";
clipprottags=(HuR hnRNPC TIA1 Auf1p);

#f=${fgroups[$(($SGE_TASK_ID-1))]};
#cpt=${clipprottags[$(($SGE_TASK_ID-1))]};
#rm -rf $f*stats*;
#for b in $bgroups; do
#	echo Doing for background group $b;
#	awk '{ print $2; }'  ~bilebi00/_KNOCKDOWNS/number_of_2_times_clip_sites | grep $cpt | while read i; do ct=`perl -e '$s=shift; $s=~/AREBP\/(.+)\//; print $1;' $i`; echo $i $ct; ~bilebi00/_KNOCKDOWNS/scripts/getFoldChange_forBinnedIdMatch.pl ~bilebi00/_KNOCKDOWNS/$i $ct $f $b $qnTable 1; done
#done
#echo Plotting $f
#~bilebi00/_KNOCKDOWNS/scripts/run_call_wrapper_barweb.sh /import/bc2/soft/app/matlab/current/Linux/  $f
#convert $f*.pdf all_$f.pdf
#exit;

#2nd run
#Controls vs Controls
#20110121
fgroups=(si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq 293wt.RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13 mRNASeq_No_4SU_No_XL_rep_B_Clip13);
bgroups="HEK293_siGFP_RNAseq";
fgroups=(HEK293_siGFP_RNAseq si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq 293wt.RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13);
bgroups="mRNASeq_No_4SU_No_XL_rep_B_Clip13";
#fgroups=(mRNASeq_No_4SU_No_XL_rep_B_Clip13 HEK293_siGFP_RNAseq si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13);
#bgroups="293wt.RNAseq";

f=${fgroups[$(($SGE_TASK_ID-1))]};
echo Getting Fold changes for $f;
rm -rf $f*stats*;
for b in $bgroups; do
	echo Doing for background group $b;
	awk '{ print $2; }'  ~bilebi00/_KNOCKDOWNS/number_of_2_times_clip_sites | while read i; do ct=`perl -e '$s=shift; $s=~/AREBP\/(.+)\//; print $1;' $i`; echo $i $ct; ~bilebi00/_KNOCKDOWNS/scripts/getFoldChange_forBinnedIdMatch.pl ~bilebi00/_KNOCKDOWNS/$i $ct $f $b $qnTable; done
done
echo Plotting $f
~bilebi00/_KNOCKDOWNS/scripts/run_call_wrapper_barweb.sh /import/bc2/soft/app/matlab/current/Linux/  $f
echo convert $f*.pdf all_$f.pdf
convert $f*.pdf all_$f.pdf
exit;

#1st run
#Controls vs KDs
#20110120
fgroups=(si_HuR_mRNASEQ HEK293_sihnRNPC_RNAseq HEK293_siTIA1_RNAseq HEK293_siAUF1_RNAseq);
bgroups="si_GFP_mRNASEQ HEK293_Nanofectin_RNAseq 293wt.RNAseq HEK293_siGFP_RNAseq mRNASeq_No_4SU_No_XL_rep_A_Clip13 mRNASeq_No_4SU_No_XL_rep_B_Clip13";
clipprottags=(HuR hnRNPC TIA1 Auf1p);

for c in `seq 0 $((${#fgroups[@]}-1))`; do 
	cpt=${clipprottags[$c]};
	f=${fgroups[$c]};
	rm -rf $f*stats*;
	for b in $bgroups; do
		echo Doing for background group $b;
		awk '{ print $2; }'  ~bilebi00/_KNOCKDOWNS/number_of_2_times_clip_sites | grep $cpt | while read i; do ct=`perl -e '$s=shift; $s=~/AREBP\/(.+)\//; print $1;' $i`; echo $i $ct; ~bilebi00/_KNOCKDOWNS/scripts/getFoldChange_forBinnedIdMatch.pl ~bilebi00/_KNOCKDOWNS/$i $ct $f $b $qnTable; done
	done
	echo Plotting $f
	matlab -nojvm -r "d=dir('$f*mean'); [r c] = size(d); for i=1:r, [p n e] = fileparts(d(i).name); fprintf('%s\n', n); addpath ~bilebi00/_KNOCKDOWNS/scripts; wrapper_barweb(n); end; exit;";
	convert $f*.pdf all_$f.pdf
done

