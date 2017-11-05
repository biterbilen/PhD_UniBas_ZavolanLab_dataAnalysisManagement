#!/bin/bash
#$ -S /bin/bash
# $ -q fs_long@@qc_nehalem
#$ -q fs_long@@high_mem_node
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.Get_peaks.sh$TASK_ID
#$ -o LOG.Get_peaks.sh.task4.n10000
#$ -t 2-2
# $ -t 1-7

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1"
}
trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=$HOME/bin:$PATH

#outdir=Analysis/ChIPseq
outdir=Analysis/ChIPseq2
mkdir -p $outdir; pushd $outdir

#TODO set
task=1
indir=~bilebi00/_EWSR1/data/clipz_ChIPseq_tags/

#1.
if [ $task == 1 ]; then
	####mkdir -p $indir; rsync mirz@web07:~/BITER/ChIP/scratch/bed_files/*bed $indir/.
	mkdir -p $indir; 
	rsync -t mirz@web08:~/BITER/ChIP/scratch/bed_files_hg19/*bed $indir/.
	echo DONE
	exit;
fi

#name=rmsk; annotateBed -i GSM693953_HEK293T.Input.U0-2.hg18.reads.bed -files ~/DATA/hg18_ucsc_tracks/$name.gtf.gz -names $name -s -counts > $name.annot
#cp ~/DATA/Shivendra/background/GSM693953_HEK293T.Input.U0-2.hg18.reads.bed.gz 
#cp ~/DATA/Shivendra/background/GSM693952_HEK293T.BACH1.U0-2.hg18.reads.bed.gz .
#gzip -d *gz
#cp ~/_DOWNLOAD/_FOR_ALIGNMENT/BEDTools-GIT-Version/bedtools/genomes/human.hg18.genome .

#2.
if [ $task == 2 ]; then
	#ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep ChIPseq | cut -f 2`);
	tag="`basename $outdir`";
	ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -w $tag | cut -f 2`);
	id=$indir/${ids[$((SGE_TASK_ID-1))]}
	echo Doing $id.bed
	#bed files are not weighted; multiple copy sequences are repeated in the file
	~bilebi00/_EWSR1/scripts/ChIPseq.fromNick/shifts_repeats.pl $id.bed
	echo DONE
	exit;
fi

#3.
if [ $task == 3 ]; then
#merge data and plot
	tag="`basename $outdir`";
	less ~/_EWSR1/data/protein_sampleid_list | grep -w $tag | cut -f 2 | while read i; do 
		n=`grep -w $i ~/_EWSR1/data/protein_sampleid_list | cut -f 3`; 
		less $i.reg_cor | perl -e '$id=shift; while(<>){ print join ("\t", $id, $_); }' $n; 
	done > all.reg_cor
	~bilebi00/bin/R --no-save --args all.reg.cor < ~bilebi00/_EWSR1/scripts/plot_shifts.R > /dev/null
	echo DONE
fi

#TODO set
n=15
n=10000

dist=BinomialDistr
#id2s=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep ChIPseq_PolIIa | cut -f 2`);
id2s=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep Input_DNA | cut -f 2`);
i2=`echo $((SGE_TASK_ID-1)) / 4 | bc`;
i2=0
id2=${id2s[$i2]}
#4
if [ $task == 4 ]; then
	ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep EWSR1_Abcam | cut -f 2`);
	i=`echo $SGE_TASK_ID % 4 | bc`;
	i=0
	id=${ids[$i]}
	echo Doing $id $id2 $n 
	out=${dist}_${id}_${id2}_n$n/

	#4.1
	ip_file=$indir/$id.bed
	wce_file=$indir/$id2.bed
	#batch 1
	ip_shift=80; #$indir/GSM693952_HEK293T.reg_cor #d should be something 80-160
	wce_shift=80; #$indir/GSM693953_HEK293T.reg_cor #d should be 80-160; check from the file
	#batch 2
	ip_shift=90; #$indir/GSM693952_HEK293T.reg_cor #d should be something 80-160
	wce_shift=65; #$indir/GSM693953_HEK293T.reg_cor #d should be 80-160; check from the file
	#takes long
	python ~bilebi00/_EWSR1/scripts/ChIPseq.fromNick/peakfinder_update.py -i $ip_file -c $wce_file -s $ip_shift -t $wce_shift -z 0 -w 500 -b 2000 -p 250 -n $n -o $out -g hg19 
	pushd $out;
	mv peaks.bed peaks.bed.back
	less peaks.bed.back | awk '$2>=0{ print }' > peaks.bed #eliminate negative start coordinates
	###python ~bilebi00/_EWSR1/scripts/ChIPseq.fromPiotr/peakfinder_by_silvia_2.py -i $ip_file -c $wce_file -s $ip_shift -t $wce_shift -z 0 -w 500 -b 2000 -p 250 -o $out -g hg19 
	~bilebi00/bsoft/bilebi00/scripts/reverse_cum_dist_of_density.pl peaks.bed 6 0 1 > background_tags_reverse_cum
	~bilebi00/bsoft/bilebi00/scripts/reverse_cum_dist_of_density.pl peaks.bed 4 0 1 > zscore_reverse_cum
	#XXX set
	mutationcount=8
	zscore=6; 
	exon_distance=10000

	#extract top CLIPped w mutationCount method -reproduced in both libaries
	zless ~bilebi00/_ARE/Analysis/Pooled_Reproducibility/EWSR1.gtf.gz | awk '$6 >= mutc && $14 == "\"2\";" { print; }' mutc=$mutationcount | \
		bedtools closest -s -d -t first -a stdin -b ~bilebi00/_ARE/Analysis/xlinkEnrichment/regions/hg19_GMAP_GENEexons_refseq.gtf.gz | \
		perl -e '$ed=shift; while(<>) { chomp; @t=split/\t/; next if ($t[17] eq "." or $t[18]>$ed); $other = "nearest_$t[17] exon_distance \"$t[18]\""; $other =~ s/([;]) /$1 nearest_/g; splice(@t,9); $t[8] = "$t[8] $other;"; print join("\t", @t), "\n"; }' $exon_distance  > EWSR1.annot.$mutationcount

	#get set of ChIPped
	less peaks.bed | awk '$5 > zscore { print; }' zscore=$zscore | \
		bedtools closest -d -t first -a stdin -b ~bilebi00/_ARE/Analysis/xlinkEnrichment/regions/hg19_GMAP_GENEexons_refseq.gtf.gz | \
		perl -e '$ed=shift; while(<>) { chomp; @t=split/\t/; next if ($t[15] eq "." or $t[16]>$ed); print join("\t",$t[0],$t[1],$t[2],$t[3],$t[4],$t[5],$t[6],$t[15]),"\n"; }' $exon_distance > \
	 	peaks.bed.zscore$zscore

	#get ChIPed and CLIPed:common
	bedtools closest -d -t first -a peaks.bed.zscore$zscore -b EWSR1.annot.$mutationcount | cut -f 1-7,9- >  peaks.bed.zscore$zscore.wCLIP

	~bilebi00/bsoft/bilebi00/scripts/reverse_cum_dist_of_density.pl peaks.bed.zscore$zscore.wCLIP 16 0 1 1 > clip_distance_cum

	#for the browser #TODO upload it to clipz
	flank=5000; background=$n; xlinkdistance=10000
	echo "min_ChIPseq_zscore=$zscore, max_Distance2xlinkSite=$xlinkdistance min_ChIPseq_background=$background" > top_peaks_wCLIP.gff3
	less peaks.bed.zscore$zscore.wCLIP | perl -e '$count=1; $flank=shift; $zscore=shift; $background=shift; $xlinkdistance=shift; while(<>){ chomp; @t=split/\t/; @clip=@t[0..6]; splice(@t,1,7); next if ($clip[4]<$zscore or $clip[6]>$background or $t[9]>$xlinkdistance); $t[1]="ChIPseq"; $t[2]="peak"; $t[3]=$clip[1]-$flank; $t[4]=$clip[2]+$flank; $t[5]=0; $t[7]="."; $t[8] = "ID=ChIPseq$count; $t[8] xlinkDistance=\"$t[9]\"; ChIP_zscore=$clip[4]; ChIP_foreground=$clip[5]; ChIP_background=$clip[6];"; pop @t; $t[8]=~s/ "/=/g; $t[8]=~s/"//g; print join("\t", @t), "\n"; $count++; }' $flank $zscore $background $xlinkdistance >> top_peaks_wCLIP.gff3

	b=`less top_peaks_wCLIP.gff3 | awk 'NR>1{print $16}' | sort | uniq | wc -l`;
	clipped=`less EWSR1.annot.$mutationcount  | awk '{ print $22;}' | sort | uniq | wc -l`;
	chipped=`zless peaks.bed.zscore$zscore | cut -f 8 | sort | uniq | wc -l`;
	all=`zless ~bilebi00/_ARE/Analysis/xlinkEnrichment/regions/hg19_GMAP_GENEexons_refseq.gtf.gz | cut -f 9 | sort | uniq | wc -l`;
	R --no-save --args $b $clipped $chipped $all < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R  | grep "^Fisher" > fisher.stat

	echo DONE
	exit;

	#4.2
	#pushd $out
	#mv peaks.bed peaks.bed.1
	#python ~bilebi00/_EWSR1/Analysis/ChIPseq/last_peakfinder_for_single_exp_nwk.py -i ip_shifted_reads.bed.sort -c wce_shifted_reads.bed.sort -z 0 -w 500 -b 2000 -p 250 -o ./
	#exit;
fi


#5.
if [ $task == 5 ] ; then 
	#if [ $SGE_TASK_ID == 1 ]; then
	for n in 10000; do
	#for n in 15 10000; do
		for id2 in ${id2s[@]}; do
			~bilebi00/bsoft/bilebi00/scripts/prep_reverse_cum_forplotting.pl "$dist*${id2}_n$n/background_tags_reverse_cum" ~bilebi00/_EWSR1/data/protein_sampleid_list ${dist}_ > $dist${id2}_background_tag_reverse_cum_n$n
			R --no-save --args $dist${id2}_background_tag_reverse_cum_n$n "$id2 tag count" "Count" < /import/bc2/home/zavolan/bilebi00/bsoft/bilebi00/scripts/plot_reverse_cum.R > /dev/null

			~bilebi00/bsoft/bilebi00/scripts/prep_reverse_cum_forplotting.pl "$dist*${id2}_n$n/zscore_reverse_cum" ~bilebi00/_EWSR1/data/protein_sampleid_list ${dist}_ > $dist${id2}_zscore_reverse_cum_n$n
			R --no-save --args $dist${id2}_zscore_reverse_cum_n$n "window zscore" "Count" < /import/bc2/home/zavolan/bilebi00/bsoft/bilebi00/scripts/plot_reverse_cum.R > /dev/null

			~bilebi00/bsoft/bilebi00/scripts/prep_reverse_cum_forplotting.pl "$dist*${id2}_n$n/clip_distance_cum" ~bilebi00/_EWSR1/data/protein_sampleid_list ${dist}_ > $dist${id2}_clip_distance_cum_n$n
			R --no-save --args $dist${id2}_clip_distance_cum_n$n "CLIP xlink distance" "Count" Cumulative < /import/bc2/home/zavolan/bilebi00/bsoft/bilebi00/scripts/plot_reverse_cum.R > /dev/null
			#fi 
			/usr/bin/convert  $dist${id2}*_n$n.pdf ${dist}_${id2}_n$n.pdf
		done
		/usr/bin/convert  ${dist}_*_n$n.pdf ${dist}_n$n.pdf
	done
fi


#TODO set
n=15
n=10000

dist=BinD
#4
if [ $task == 6 ]; then
	id2s=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep ChIPseq_PolIIa | cut -f 2`);
	ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep ChIPseq_EWSR1 | cut -f 2`);
#	id2s=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep Input_DNA_ChIP_HeLa| cut -f 2`);
#	ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep EWSR1_ChIP_Abcam_HeLa | cut -f 2`);
	out=${dist}_n$n/

	#4.1
	wce_file=$indir/679.bed,$indir/680.bed,$indir/711.bed
	ip_file=$indir/648.bed,$indir/649.bed,$indir/650.bed,$indir/655.bed
	wce_shift=150,150,150; #$indir/GSM693953_HEK293T.reg_cor #d should be 80-160; check from the file
	ip_shift=160,160,160,160; #$indir/GSM693952_HEK293T.reg_cor #d should be something 80-160
	echo Doing $ip_file $wce_file
	python ~bilebi00/_EWSR1/scripts/ChIPseq.fromNick/peakfinder_update.py -i $ip_file -c $wce_file -s $ip_shift -t $wce_shift -z 0 -w 500 -b 2000 -p 250 -n $n -o $out -g hg19 
	pushd $out;
	mv peaks.bed peaks.bed.back
	less peaks.bed.back | awk '$2>=0{ print }' > peaks.bed #eliminate negative start coordinates
	###python ~bilebi00/_EWSR1/scripts/ChIPseq.fromPiotr/peakfinder_by_silvia_2.py -i $ip_file -c $wce_file -s $ip_shift -t $wce_shift -z 0 -w 500 -b 2000 -p 250 -o $out -g hg19 
	~bilebi00/bsoft/bilebi00/scripts/reverse_cum_dist_of_density.pl peaks.bed 6 0 1 > background_tags_reverse_cum
	~bilebi00/bsoft/bilebi00/scripts/reverse_cum_dist_of_density.pl peaks.bed 4 0 1 > zscore_reverse_cum
	less peaks.bed | awk '$5 > 1 { print; }' > peaks.bed.zscore1
	zless ~bilebi00/_ARE/Analysis/Pooled_Reproducibility/EWSR1.gtf.gz | awk '$6 >= 32 { print; }' | closestBed -s -d -t first -a stdin -b ~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/hg19_GMAP_EXONS.refseq.gtf | perl -e 'while(<>) { chomp; @t=split/\t/; $other = "nearest_$t[17] exon_distance \"$t[18]\""; $other =~ s/([;]) /$1 nearest_/g; splice(@t,9); $t[8] = "$t[8] $other;"; print join("\t", @t), "\n"; }'  > EWSR1.annot.32
	closestBed -d -t first -a peaks.bed.zscore1 -b EWSR1.annot.32 | grep -v "\-1$"  > peaks.bed.zscore1.wCLIP
	~bilebi00/bsoft/bilebi00/scripts/reverse_cum_dist_of_density.pl peaks.bed.zscore1.wCLIP 16 0 1 1 > clip_distance_cum

	flank=5000; zscore=2; background=15; xlinkdistance=10000;
	echo "min_ChIPseq_zscore=$zscore, max_Distance2xlinkSite=$xlinkdistance min_ChIPseq_background=$background" > top_peaks_wCLIP.gff3
	less peaks.bed.zscore1.wCLIP | perl -e '$count=1; $flank=shift; $zscore=shift; $background=shift; $xlinkdistance=shift; while(<>){ chomp; @t=split/\t/; @clip=@t[0..6]; splice(@t,1,7); next if ($clip[4]<$zscore or $clip[6]>$background or $t[9]>$xlinkdistance); $t[1]="ChIPseq"; $t[2]="peak"; $t[3]=$clip[1]-$flank; $t[4]=$clip[2]+$flank; $t[7]="."; $t[8] = "ID=ChIPseq$count; $t[8] xlinkDistance \"$t[9]\"; ChIP_zscore=$clip[4]; ChIP_foreground=$clip[5]; ChIP_background=$clip[6];"; pop @t; $t[8]=~s/ "/=/g; $t[8]=~s/"//g; print join("\t", @t), "\n"; $count++; }' $flank $zscore $background $xlinkdistance >> top_peaks_wCLIP.gff3

	echo DONE
	exit;

	#4.2
	#pushd $out
	#mv peaks.bed peaks.bed.1
	#python ~bilebi00/_EWSR1/Analysis/ChIPseq/last_peakfinder_for_single_exp_nwk.py -i ip_shifted_reads.bed.sort -c wce_shifted_reads.bed.sort -z 0 -w 500 -b 2000 -p 250 -o ./
	#exit;
fi


