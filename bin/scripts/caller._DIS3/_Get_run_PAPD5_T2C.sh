#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@high_mem_node
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG.normalization 

#TODO write the cluster submission version; not ready yet! todo this:
#TODO FIXME Cannot plot in nodes
#TODO some parts should be done manually like the selection of cut-offs
#TODO separate into several files (i.e. write the caller script for this pipeline

tag=MC;
tag=PC;
tag=PAPD5;
tag=PAPD5_25nc;
tag=PAPD5_T2C;
libs=`less ~bilebi00/_DIS3/data/protein_sampleid_list_$tag | cut -f 2`;
#TODO change this as output_directory

outdir=Normalization_$tag #output_directory
mkdir -p $outdir; cd $outdir;

#TODO set
dbpref=DB_copies_
dbpref=DB_

#set pathnames
number_samples=0
copies="";
for i in $libs; do
	ln -sf ~bilebi00/_DIS3/$dbpref$i/reverse*100 raw_$number_samples;
	copies="$copies /import/bc2/home/zavolan/bilebi00/_DIS3/$dbpref$i/copy*"
	number_samples=$((number_samples+1));
done	

#TODO compile and generalize the function plot_reverse_cum_2
#plot revcum of raw values
matlab -nojvm -r "addpath ~bilebi00/_DIS3/scripts; plot_reverse_cum('', 'raw*', 100); exit;";
matlab -nojvm -r "addpath ~bilebi00/_PAPD5/scripts; plot_reverse_cum('', 'raw*', 100); exit;";
i=`ls reverse*pdf`; mv $i prenorm_$i;

#----------
#FIXME table joining takes long time; around 1-2 hours
#TODO generalize script to to get library
nohup ~bilebi00/_PAPD5/scripts/tablize.pl 5 table.inp 1 $copies &> ti.log &
nohup ~bilebi00/_PAPD5/scripts/tablize.pl 5 tableall.inp 0 $copies &> tai.log &

#----------
#FIXME table trimming takes 15-20 minutes; TODO merge table summing and trimming
#prepare non overlapping windows; 
cut=0.01 #TODO decide cut-off value here 0.01
nohup ~bilebi00/_PAPD5/scripts/trim_min_table_values_for_normalization.pl table.inp $cut > table 2> tt.log &
nohup ~bilebi00/_PAPD5/scripts/trim_min_table_values_for_normalization.pl tableall.inp $cut > tableall 2> tat.log &

#----------
nohup less xaa | perl -e '$ns=shift; $ns--; while(<>) { chomp; @t=split; for(0..$ns) { $s{$_} += $t[$_+3]; }} for(0..$ns) { print $s{$_}, "\n"; } ;' $number_samples > xaa.totals 2> tt2.xaa.log &
nohup less table | perl -e '$ns=shift; $ns--; while(<>) { chomp; @t=split; for(0..$ns) { $s{$_} += $t[$_+3]; }} for(0..$ns) { print $s{$_}, "\n"; } ;' $number_samples > table.totals 2> tt2.log &
nohup less tableall | perl -e '$ns=shift; $ns--; while(<>) { chomp; @t=split; for(0..$ns) { $s{$_} += $t[$_+3]; }} for(0..$ns) { print $s{$_}, "\n"; } ;' $number_samples > tableall.totals 2> tat2.log & 

#----------
#scale normalization T=1M
n=1000000;
nohup less table | perl -e '$ns=shift; $ns--; $tf=shift; $T=shift; open F, $tf or die; @tots=<F>; chomp @tots; close F; while(<>) { chomp; @t=split; @k=@t; for(0..$ns)  { splice(@k,$_+3,0,$t[$_+3]/$tots[$_]*$T); } print join("\t", @k), "\n"; };' $number_samples table.totals $n > table.norm_and_raw 2> s.log &
nohup less tableall | perl -e '$ns=shift; $ns--; $tf=shift; $T=shift; open F, $tf or die; @tots=<F>; chomp @tots; close F; while(<>) { chomp; @t=split; @k=@t; for(0..$ns)  { splice(@k,$_+3,0,$t[$_+3]/$tots[$_]*$T); } print join("\t", @k), "\n"; };' $number_samples tableall.totals $n > tableall.norm_and_raw 2> sall.log &

##power law normalization
##transform data to a dist with T=1M, m=1.25
##TODO set cut-offs
lc=1
uc=10000
lc=.1
uc=10000
lc=.1
uc=100
nohup ~bilebi00/_PAPD5/scripts/piotrNormalizationScripts/regress_raw.pl $number_samples $lc $uc > fits 2> fits.log &
#----------
#scale normalization T=1M
#get normalized values using overlapping and non-overllapping tables
nohup ~bilebi00/_PAPD5/scripts/piotrNormalizationScripts/normalize_level1.pl table $number_samples &> norm.log &
nohup ~bilebi00/_PAPD5/scripts/piotrNormalizationScripts/normalize_level1.pl tableall $number_samples &> normtal.log &

#prepare rev_cum post norm
#FIXME This part is slow too

endi=$((number_samples-1));
i=0; #TODO set this to 0 normally
for index in `seq $i $endi`; do
	echo $index;
	less table.norm_and_raw | perl -e '$index = shift; while(<>) { ($c, $s, $b, @t) = split; $e=99+$b; print join("\t", "CC", $c, $s, $b, $e, $t[$index]), "\n";  }' $index > $index.inp ;
	sort -k 6,6gr $index.inp > $index.sorted;
	~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $index.sorted 5 > norm_$index;
done

#plot rev_cum norm values #TODO use the generalized plot_reverse_cum_2
matlab -nojvm -r "addpath ~bilebi00/_DIS3/scripts; plot_reverse_cum('', 'norm*', 100); exit;";
matlab -nojvm -r "addpath ~bilebi00/_PAPD5/scripts; plot_reverse_cum('', 'norm*', 100); exit;";
i=`ls reverse*pdf`; mv $i postnorm_$i;

#Error model 
#estimate sigmas from non-overlapping windows
#TODO copy fileOperations.pm and fit_sigmas.pl in scripts directory
#TODO run it from the scripts directory
#TODO manually set the hash content in the script for nsamp and ncols and push the comparison pairs if you want to compare all
#Fits sigmas using the non-overlapping windows
#takes ~1-2 hour
lc2=2.1
uc2=10000
nohup ~bilebi00/_DIS3/scripts/fit_sigmas.pl $tag $lc2 $uc2 > ~bilebi00/_DIS3/$outdir/noise_sigma.$lc2.$uc2 2> ~bilebi00/_DIS3/$outdir/noise_sigma.$lc2.$uc2.log & 

#TODO select sigma^2 from noise_sigma file
#variation is higher than PAPD5 for DIS3L2
s2=0.289307 #MC_DIS3L2
s2=0.313721 #MC_DIS3L2 scalae normalization variance
s2=0.166626 #PAPD5
s2=0.167847 #PAPD5_25nc
s2=0.0241852 #TIA1 T2C
s2=0.0631714 #PAPD5_T2C 2.5 cutoff
#run z-scores for all the comparison pairs
#takes ~ 1 hour in high_mem_node
endi=$((number_samples-1));
for i in `seq 0 $endi`; do 
	k=$((i+1));
	for j in `seq $k $endi`; do 
		echo $i$j;
		nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "$i $j" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.$i${j}zscore 2> z$i$j.log &
	done
done
#----------
#TODO set the indices by ls data which is the table order
#MC
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 1" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.01zscore 2> z01.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 2" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.02zscore 2> z02.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 3" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.03zscore 2> z03.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 4" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.04zscore 2> z04.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "1 2" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.12zscore 2> z12.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "1 3" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.13zscore 2> z13.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "1 4" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.14zscore 2> z14.log &
#PC
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 1" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.01zscore 2> z01.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 4" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.04zscore 2> z04.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 5" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.05zscore 2> z05.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "0 6" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.06zscore 2> z06.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "1 4" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.14zscore 2> z14.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "1 5" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.15zscore 2> z15.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "1 6" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.16zscore 2> z16.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "2 4" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.24zscore 2> z24.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "2 5" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.25zscore 2> z25.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "2 6" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.26zscore 2> z26.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "3 4" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.34zscore 2> z34.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "3 5" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.35zscore 2> z35.log &
nohup ~bilebi00/_PAPD5/scripts/zScore_multnoise.pl "3 6" tableall.norm_and_raw $number_samples 3 $s2 > tableall.norm_and_raw.36zscore 2> z36.log &
#PAPD5
#runned the loop and deleted unrelate comparisons

#----------
#MC
for i in 01 02 03 04 12 13 14; do nohup less tableall.norm_and_raw.${i}zscore | grep -v NA > tableall.norm_and_raw.${i}zscore.woNA 2> yy${i}woNA.log & done
#PC
for i in 01 04 05 06 14 15 16 24 25 26 34 35 36; do nohup less tableall.norm_and_raw.${i}zscore | grep -v NA > tableall.norm_and_raw.${i}zscore.woNA 2> yy${i}woNA.log & done
#PAPD5
for i in 01 23 04 14 02 12 03 13; do nohup less tableall.norm_and_raw.${i}zscore | grep -v NA > tableall.norm_and_raw.${i}zscore.woNA 2> yy${i}woNA.log & done

#plot zscore histograms
#MC
matlab -nojvm -r "addpath ~bilebi00/_DIS3/scripts; plot_hist_zscore_multnoise('./'); exit;";
#PC
matlab -nojvm -r "addpath ~bilebi00/_DIS3/scripts; plot_hist_zscore_multnoise_2('./'); exit;";
#PAPD5
matlab -nojvm -r "addpath ~bilebi00/_PAPD5/scripts; plot_hist_zscore_multnoise('./'); exit;";
matlab -nojvm -r "addpath /import/bc2/home/zavolan/bilebi00/_PAPD5/scripts; plot_hist_zscore_multnoise('./'); exit;";

#TODO inspect z-scores from the plot and get cutoff values
zcut=3;
zcut=2; #PAPD5
zcut=0; #PAPD5 T2C
zcut=3; #PAPD5
zcut=0; #PAPD5 T2C
for i in tableall.norm_and_raw.*zscore.woNA; do
	nohup less $i |perl -e '$zcut=shift; while(<>){@t=split; if ($t[5] > $zcut or $t[5] < -$zcut ) { print $_;} }' $zcut > $i.filtered  2> $i.log &
done

#DIS3L2 MC
#pairs
chrindex=1
iis=(12 12 13)
jjs=(13 14 14)
for i in `seq 0 $((${#iis[@]}-1))`; do 
	nohup ~bilebi00/_PAPD5/scripts/wrapper_for_innerJoin_with_chrs.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA tableall.norm_and_raw.${jjs[$i]}zscore.woNA 1 1 '1,2,3,6,8,9,12' $chrindex > _${iis[$i]}_${jjs[$i]}.notfiltered 2> zz_${iis[$i]}_${jjs[$i]}.log &
	#nohup innerJoin.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA.filtered tableall.norm_and_raw.${jjs[$i]}zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _${iis[$i]}_${jjs[$i]}.filtered 2> zz_${iis[$i]}_${jjs[$i]}.log &
done

#DIS3L2 PC
#pairs
chrindex=1
iis=(01 04 04 05 14 14 15 24 24 25 34 34 35)
jjs=(01 05 06 06 15 16 16 25 26 26 35 36 36)
for i in `seq 0 $((${#iis[@]}-1))`; do 
	nohup ~bilebi00/_PAPD5/scripts/wrapper_for_innerJoin_with_chrs.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA tableall.norm_and_raw.${jjs[$i]}zscore.woNA 1 1 '1,2,3,6,8,9,12' $chrindex > _${iis[$i]}_${jjs[$i]}.notfiltered 2> zz_${iis[$i]}_${jjs[$i]}.log &
	#nohup innerJoin.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA.filtered tableall.norm_and_raw.${jjs[$i]}zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _${iis[$i]}_${jjs[$i]}.filtered 2> zz_${iis[$i]}_${jjs[$i]}.log &
done

#PAPD5
#add other interesting comparisons with TIA1-40
#pairs
iis=(01 04 02 03)
jjs=(01 14 12 13)
for i in `seq 0 $((${#iis[@]}-1))`; do 
	nohup ~bilebi00/_PAPD5/scripts/wrapper_for_innerJoin_with_chrs.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA tableall.norm_and_raw.${jjs[$i]}zscore.woNA 1 1 '1,2,3,6,8,9,12' $chrindex > _${iis[$i]}_${jjs[$i]}.notfiltered 2> zz_${iis[$i]}_${jjs[$i]}.log &
	#nohup innerJoin.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA.filtered tableall.norm_and_raw.${jjs[$i]}zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _${iis[$i]}_${jjs[$i]}.filtered 2> zz_${iis[$i]}_${jjs[$i]}.log &
done

#PAPD5 copy and t2c
chrindex=1
iis=(04 14)
jjs=(04 14)
for i in `seq 0 $((${#iis[@]}-1))`; do 
	nohup ~bilebi00/_PAPD5/scripts/wrapper_for_innerJoin_with_chrs.pl tableall.norm_and_raw.${iis[$i]}zscore.woNA.copy tableall.norm_and_raw.${jjs[$i]}zscore.woNA.t2c 1 1 '1,2,3,6,8,9,12' $chrindex > _${iis[$i]}_${jjs[$i]}.copy_t2c.notfiltered 2> zz_${iis[$i]}_${jjs[$i]}.log &
done

#MC
innerJoin.pl tableall.norm_and_raw.02zscore.woNA.filtered tableall.norm_and_raw.03zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _02_03.filtered 2> zz_02_03.log &
innerJoin.pl tableall.norm_and_raw.02zscore.woNA.filtered tableall.norm_and_raw.04zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _02_04.filtered 2> zz_02_04.log &
innerJoin.pl tableall.norm_and_raw.03zscore.woNA.filtered tableall.norm_and_raw.04zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _03_04.filtered 2> zz_03_04.log &
innerJoin.pl tableall.norm_and_raw.12zscore.woNA.filtered tableall.norm_and_raw.13zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _12_13.filtered 2> zz_12_13.log &
innerJoin.pl tableall.norm_and_raw.12zscore.woNA.filtered tableall.norm_and_raw.14zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _12_14.filtered 2> zz_12_14.log &
innerJoin.pl tableall.norm_and_raw.13zscore.woNA.filtered tableall.norm_and_raw.14zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _13_14.filtered 2> zz_13_14.log &

#PC
innerJoin.pl tableall.norm_and_raw.01zscore.woNA.filtered tableall.norm_and_raw.01zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _01_01.filtered 2> zz_01_01.log &
innerJoin.pl tableall.norm_and_raw.04zscore.woNA.filtered tableall.norm_and_raw.05zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _04_05.filtered 2> zz_04_05.log &
innerJoin.pl tableall.norm_and_raw.04zscore.woNA.filtered tableall.norm_and_raw.06zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _04_06.filtered 2> zz_04_06.log &
innerJoin.pl tableall.norm_and_raw.05zscore.woNA.filtered tableall.norm_and_raw.06zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _05_06.filtered 2> zz_05_06.log &
innerJoin.pl tableall.norm_and_raw.14zscore.woNA.filtered tableall.norm_and_raw.15zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _14_15.filtered 2> zz_14_15.log &
innerJoin.pl tableall.norm_and_raw.14zscore.woNA.filtered tableall.norm_and_raw.16zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _14_16.filtered 2> zz_14_16.log &
innerJoin.pl tableall.norm_and_raw.15zscore.woNA.filtered tableall.norm_and_raw.16zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _15_16.filtered 2> zz_15_16.log &
innerJoin.pl tableall.norm_and_raw.24zscore.woNA.filtered tableall.norm_and_raw.25zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _24_25.filtered 2> zz_24_25.log &
innerJoin.pl tableall.norm_and_raw.24zscore.woNA.filtered tableall.norm_and_raw.26zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _24_26.filtered 2> zz_24_26.log &
innerJoin.pl tableall.norm_and_raw.25zscore.woNA.filtered tableall.norm_and_raw.26zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _25_26.filtered 2> zz_25_26.log &
innerJoin.pl tableall.norm_and_raw.34zscore.woNA.filtered tableall.norm_and_raw.35zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _34_35.filtered 2> zz_34_35.log &
innerJoin.pl tableall.norm_and_raw.34zscore.woNA.filtered tableall.norm_and_raw.36zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _34_36.filtered 2> zz_34_36.log &
innerJoin.pl tableall.norm_and_raw.35zscore.woNA.filtered tableall.norm_and_raw.36zscore.woNA.filtered 1 1 '1,2,3,6,8,9,12' > _35_36.filtered 2> zz_35_36.log &

##Do not run in clusters
##TODO write the looped version; alternative inner Join for big tables
#~bilebi00/_PAPD5/scripts/wrapper_for_innerJoin_with_chrs.pl tableall.norm_and_raw.04zscore.woNA tableall.norm_and_raw.14zscore.woNA 1 1 '1,2,3,6,8,9,12' > _04_14 &

#rest is done by 
# run_webprep_for_filtered.pl
# _Get_last4nc_from3Prime.sh
# _Get_copy_t2c_intersection.sh
# _Get_annotations_stats.sh
# _Get_all_annotations.sh
# scripts/fisher.ex.t.R
# scripts/calcNucPosFreq.pl for AUaddition
# TABLE for manuscript
c=5; c2=5; d=-5; d2=-5;
annot=_04_14.csv
outfile=TableS1.csv
less $annot | head -n 1 > $outfile
less $annot | awk 'NR>1{ if ($5>c && $6>c2 || $5<d && $6<d2) print $0; }' c=$c c2=$c2 d=$d d2=$d2 | sort -k 5,6gr >> $outfile

#Tables for copy and T2C analysis from individual replicates
c=4; c2=1.5; d=-4; d2=-1.5;
annot=../Normalization_PAPD5_copyAndT2C/_04_04.csv
outfile=04_copy_t2c_TableS1.csv
less $annot | head -n 1 > $outfile
less $annot | awk 'NR>1{ if ($5>c && $6>c2 || $5<d && $6<d2) print $0; }' c=$c c2=$c2 d=$d d2=$d2 | sort -k 5,6gr >> $outfile
#---
annot=../Normalization_PAPD5_copyAndT2C/_14_14.csv
outfile=14_copy_t2c_TableS1.csv
less $annot | head -n 1 > $outfile
less $annot | awk 'NR>1{ if ($5>c && $6>c2 || $5<d && $6<d2) print $0; }' c=$c c2=$c2 d=$d d2=$d2 | sort -k 5,6gr >> $outfile

#-----------------

#ATTENTION this part is wrong doesnot consider the second z-score
#Bed format change to check in the browser
less _04_14.filtered | perl -e 'while(<>) { $_=~/(chr\w+)_([\-+])_(\d+)\t([\-\.\d]*)/;@t=split; print join ("\t",$1,$3,$3+99,"PAPD5_IGF",$4,$2), "\n" if($4>0);}' > _04_14.bed &
less _24_34.filtered | perl -e 'while(<>) { $_=~/(chr\w+)_([\-+])_(\d+)\t([\-\.\d]*)/;@t=split; print join ("\t",$1,$3,$3+99,"TIA_IGF",$4,$2), "\n" if($4>0);}' > _24_34.bed &
less _03_13.filtered | perl -e 'while(<>) { $_=~/(chr\w+)_([\-+])_(\d+)\t([\-\.\d]*)/;@t=split; print join ("\t",$1,$3,$3+99,"PAPD5_TIA1",$4,$2), "\n" if($4>0);}' > _03_13.bed &

#Bed format change to check in the browser for depleted set
less _04_14.filtered | perl -e 'while(<>) { $_=~/(chr\w+)_([\-+])_(\d+)\t([\-\.\d]*)/;@t=split; print join ("\t",$1,$3,$3+99,"PAPD5_IGF",$4,$2), "\n" if($4<0);}' > _40_41.bed &

#Annotate
cat _03_13.filtered | awk '{if($2 > 0) print $0}' > enriched_replicates_vs_TIA1.z5
cat enriched_replicates_vs_TIA1.z5 | perl -e 'while(<>) {@s = split(/\s+/, $_); @c = split(/\_/, $s[0]); print "$s[0]\t$c[0]\t$c[1]\t$c[2]\t", $c[2]+99, "\t$s[1]\n";}' > enriched_replicates_vs_TIA1.z5.summary
~/../GROUP/miRNA/CentralProg/match_repeats enriched_replicates_vs_TIA1.z5.summary /import/bc2/home/zavolan/rodak/auxDataNewest/1/auxData/repeatMasks _rmsk.filtered 20 > enriched_replicates_vs_TIA1.z5.repoverlap20
grep rRNA enriched_replicates_vs_TIA1.z5.repoverlap20 | awk '{print $1}' | sort | uniq -c | wc -l

cat _04_14.filtered | awk '{if($2 < 0) print $0}' > enriched_replicates_vs_TIA1.z-5
cat enriched_replicates_vs_TIA1.z-5 | perl -e 'while(<>) {@s = split(/\s+/, $_); @c = split(/\_/, $s[0]); print "$s[0]\t$c[0]\t$c[1]\t$c[2]\t", $c[2]+99, "\t$s[1]\n";}' > enriched_replicates_vs_TIA1.z-5.summary
~/../GROUP/miRNA/CentralProg/match_repeats enriched_replicates_vs_TIA1.z-5.summary /import/bc2/home/zavolan/rodak/auxDataNewest/1/auxData/repeatMasks _rmsk.filtered 20 > enriched_replicates_vs_TIA1.z-5.repoverlap20
grep rRNA enriched_replicates_vs_TIA1.z-5.repoverlap20 | awk '{print $1}' | sort | uniq -c | wc -l

#vs IGF
cat _04_14.filtered | awk 'BEGIN{OFS="\t";}{if($4 > 0) print $1,$4}' > enriched_replicates_vs_IGF.z6
cat enriched_replicates_vs_IGF.z6 | perl -e 'while(<>) {@s = split(/\s+/, $_); @c = split(/\_/, $s[0]); print "$s[0]\t$c[0]\t$c[1]\t$c[2]\t", $c[2]+99, "\t$s[1]\n";}' > enriched_replicates_vs_IGF.z6.summary
~/../GROUP/miRNA/CentralProg/match_repeats enriched_replicates_vs_IGF.z6.summary /import/bc2/home/zavolan/rodak/auxDataNewest/1/auxData/repeatMasks _rmsk.filtered 20 > enriched_replicates_vs_IGF.z6.repoverlap20
grep rRNA enriched_replicates_vs_IGF.z6.repoverlap20 | awk '{print $1}' | sort | uniq -c | wc -l

cat _04_14.filtered | awk 'BEGIN{OFS="\t";}{if($4 < 0) print $1,$4}' > enriched_replicates_vs_IGF.z-6
cat enriched_replicates_vs_IGF.z-6 | perl -e 'while(<>) {@s = split(/\s+/, $_); @c = split(/\_/, $s[0]); print "$s[0]\t$c[0]\t$c[1]\t$c[2]\t", $c[2]+99, "\t$s[1]\n";}' > enriched_replicates_vs_IGF.z-6.summary
~/../GROUP/miRNA/CentralProg/match_repeats enriched_replicates_vs_IGF.z-6.summary /import/bc2/home/zavolan/rodak/auxDataNewest/1/auxData/repeatMasks _rmsk.filtered 20 > enriched_replicates_vs_IGF.z-6.repoverlap20
grep rRNA enriched_replicates_vs_IGF.z-6.repoverlap20 | awk '{print $1}' | sort | uniq -c | wc -l


#significance with Fisher's Exact Test
#> a <- matrix(c(2,29,30,71), nr=2);
#> fisher.test(a, 'greater')
#
#Fisher's Exact Test for Count Data
#
#data:  a 
#p-value = 0.00779
#alternative hypothesis: true odds ratio is not equal to 1 
#95 percent confidence interval:
#0.01796179 0.72391587 
#sample estimates:
#odds ratio 
#0.1648749 

#> a
#[,1] [,2]
#[1,]    2   30
#[2,]   29   71

