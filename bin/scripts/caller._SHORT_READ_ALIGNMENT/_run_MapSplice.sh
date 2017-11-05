#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
# $ -q fs_short@@qc_nehalem
#$ -q fs_long@@qc_nehalem 
#$ -P project_zavolan
# $ -l h_vmem=60G
#$ -j y
#$ -cwd
#$ -o LOG._run_MapSplice.chr19.rep_A

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

chr_dir=~/DATA/hg18
bti_dir=~/DATA/Bowtie_DB_INDEX
chr=chr19 #""
#db_name=_hg18_$chr #hg_TR #TODO set here
db_name=hg18_$chr #hg_TR #TODO set here
tag=A
outdir=rep_$tag.$db_name/MapSplice
len=37;
seglen=$(($len / 2));
t=8

mkdir -p $outdir; cd $outdir;

inp=inp.fa

less ~bilebi00/CLIP_samples/mRNASeq_No_4SU_No_XL_rep_$tag/PROTS/mRNASeq_No_4SU_No_XL_rep_$tag.bwa_oligoAln | perl -e '$l=shift; $c=shift; while(<>) {($id)=split;$chr=<>;chomp $chr;$_=<>;($x,$y,$z,$str)=split; $seq=<>; chomp $seq; $_=<>; $_=<>; $_=<>; if ($chr =~ /chr/) { $seq=~tr/ACGT/TGCA/ if ($str eq "-"); $seq =~ s/\-//g; print ">$id\t$seq\n" if (length($seq) == $l and ($c eq "" or $c eq $chr)); }} ' $len $chr > $inp
less ~bilebi00/_SHORT_READ_ALIGNMENT/data/rep_${tag}_unmapped.seq |perl -e '$l=shift; $_=<>;  while(<>) { @t=split; print ">sequd_$t[0]\t$t[1]\n" if (length($t[1]) == $l); }' $len >> $inp

less $inp | sort -u | perl -e 'while(<>){@t=split; print $t[0], "\n", $t[1], "\n"; }' > $inp.tmp
mv $inp.tmp $inp

echo "python ~bilebi00/bsoft/MapSplice_1.15.1/bin/mapsplice_segments.py -i 50 -x 500000 -X $t -Q fa -o ./ -w $len -c $chr_dir -u $inp -B $bti_dir/$db_name -L $seglen --non-canonical --fusion-non-canonical --not-rem-temp 2> log"
python ~bilebi00/bsoft/MapSplice_1.15.1/bin/mapsplice_segments.py -i 50 -x 500000 -X $t -Q fa -o ./ -w $len -c $chr_dir -u $inp -B $bti_dir/$db_name -L $seglen --non-canonical --fusion-non-canonical --not-rem-temp 2> log


#cp ~bilebi00/_SHORT_READ_ALIGNMENT/data/mapsplice_ex/1M_36bp_fastq.txt $inp
#echo "python ~bilebi00/_DOWNLOAD/_FOR_ALIGNMENT/MapSplice_1.15.1/bin/mapsplice_segments.py -i 50 -x 500000 -X $t -Q fq -o ./ -w 36 -c $chr_dir -u $inp -B $bti_dir/hg18 -L 18 --non-canonical --fusion-non-canonical 2> 36bp_time.log"
#python ~bilebi00/_DOWNLOAD/_FOR_ALIGNMENT/MapSplice_1.15.1/bin/mapsplice_segments.py -i 50 -x 500000 -X $t -Q fq -o ./ -w 36 -c $chr_dir -u $inp -B $bti_dir/hg18 -L 18 --non-canonical --fusion-non-canonical 2> 36bp_time.log

#python ~bilebi00/_DOWNLOAD/_FOR_ALIGNMENT/MapSplice_1.15.1/bin/mapsplice_segments.py -Q fq -o $outdir -w 36 -c ~/_SHORT_READ_ALIGNMENT/data/mapsplice_ex/ -u ~/_SHORT_READ_ALIGNMENT/data/mapsplice_ex/1M_36bp_fastq.txt -B ~/_SHORT_READ_ALIGNMENT/data/mapsplice_ex/index -L 18
