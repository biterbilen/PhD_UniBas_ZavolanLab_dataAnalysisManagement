#!/bin/bash

index=$1; #starts from 0; select from libraries

#TODO set here the input dir name for T2C peaks
profile_region=TIA_ENHANCED_CE;
#profile_region=POLYA_217_single_polyA;
#profile_region=POLYA_217;
#profile_region=POLYA;
#profile_region=AGO;

t2c_peak_indir=(__69_nhood10000 __70_nhood10000 DB_67-68-212_nhood10000);
#if [ "$profile_region" == "AGO" ] && [ $index == 3 ] ; then exit; fi
t2c_peak_file_nametag=${t2c_peak_indir[$index]};

top_t2c_density=10000;
profile=SS;
prof_dist_file_pat="%s%d%s%s%c%d%d%s%s%s%c%d%d%f"; #TODO could be different in different types of data
prof_index=8; #TODO could be different in different types of data
dist_index=2; #TODO could be different in different types of data
t2c_index=14; #TODO could be different in different types of data
pwd=`pwd`;
rand=$RANDOM;
inp=T2C_region.inp;
cls=$inp.cls
nhood=10000;
prof_file=$cls.nhood$nhood.$profile_region.profile;
prof_dist_file=${prof_file}_distribution

#TODO set here the names of the region files and 
if [ $profile_region == POLYA ]; then
	region_file1=PolyA-DB-mammal-human-hg18
elif [ $profile_region == AGO ]; then
	region_file1=DB_67-68-212_nhood10000
elif [ $profile_region == "POLYA_217" ]; then
	region_file1=PolyA_217
elif [ $profile_region == "POLYA_217_single_polyA" ]; then
	region_file1=PolyA_217_single_polyA
elif [ $profile_region == "TIA_ENHANCED_CE" ]; then
	region_file1=TIA_ENHANCED_CE
fi
#region_file2=PolyA-DB-mammal-human-hg18 #TODO change this if there are multiple regions in the input file; can be maximum 2 for the moment

t2c_peak_infile="$pwd/$t2c_peak_file_nametag/T2C_exon_cls.cls";
t2c_density_file="T2C_density";

#TODO select here the T2Cs peak files and region files and set the $outdir name according to that
#region_files="$region_file1 region_file2"
region_files="$region_file1";
outdir=rp_`echo $region_files | sed 's/ //g'`T2C_$t2c_peak_file_nametag;

mkdir -p $outdir; cd $outdir;

#TODO set here the region file content
#remark: start coordinates in UCSC is 0 based
if [ $profile_region == "POLYA_217" ] || [ $profile_region == "POLYA_217_single_polyA" ]; then
	if [ $profile_region == "POLYA_217_single_polyA" ]; then
		less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/polyAsites | grep -v NULL | perl -e '$pr=shift; while(<>){ chomp; @t = split; if ($t[7] == 1) { $t[6]=$pr; print join("\t", @t), "\n"; } }' $profile_region > $region_file1;
	else
		less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/polyAsites | grep -v NULL > $region_file1;
	fi
	prof_dist_file_pat="%s%d%s%s%c%d%d%f%s%d%s%s%c%d%d%f"; #TODO could be different in different types of data
	prof_index=9; #TODO could be different in different types of data
	t2c_index=16; #TODO could be different in different types of data
elif [ $profile_region == POLYA ]; then
	less ~bilebi00/DATA/PolyA_DB_mammal_human_hg18.gz | awk '{OFS="\t";} { if (NR > 1 ) { print $5, $2, $7, $3+1, $4, "POLYA"; } }' > $region_file1;
elif [ $profile_region == AGO ]; then
	less $pwd/$region_file1/T2C_exon_cls.cls | grep $profile | sort -gr -k 6,6 | head -n $top_t2c_density | perl -e ' $profile_reg=shift; while(<>){ $b=2; $e=9; @t=split; if($t[2] == "-") { $b=9; $e=2;} $t[3]-=$b; $t[4]+=$e; $t[0]=$t[5]; $t[5]=$profile_reg; print join("\t", @t, "\n"); } ' $profile_region > $region_file1;
elsif [ $profile_region == TIA_ENHANCED_CE ]; then
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/TIA_ENHANCED_CE > $region_file1;
fi

#TODO set this if there are multiple type of regions (eg. ALTERNATIVE, CONSTITUTIVE)
#less ~bilebi00/DATA/PolyA_DB_mammal_human_hg18.gz | awk '{OFS="\t";} { if (NR > 1 ) { print "POLYA", $2, $7, $3+1, $4, $6; } }' > $regionfile2; 

less $t2c_peak_infile | grep $profile | sort -gr -k 6,6 | head -n $top_t2c_density | perl -e '@a=(); $k=0; while(<>){ @t=split; print join("\t", @t, "\n"); push(@a, $t[$#t]); } print STDERR $a[int($#a/2)],"\n";' > $inp 2> $t2c_density_file

#clustering input file
cat $region_files >> $inp

#cluster region and T2C peaks
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls;

#get profile of T2C around region
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nhood $cls $profile_region $profile 1 2> $rand;
grep "^Distance" $rand > $prof_dist_file;
rm $rand;

#plot now
cd $pwd;
echo "addpath('~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts'); plot_profile_error('$outdir',$nhood,$top_t2c_density,'$outdir/$prof_dist_file','$outdir/$t2c_density_file','$prof_dist_file_pat','$profile_region','$profile_region','$profile_region',1, $prof_index, $dist_index, $t2c_index);exit;";
matlab -nojvm -r "addpath('~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts'); plot_profile_error('$outdir',$nhood,'$outdir/$prof_dist_file','$outdir/$t2c_density_file','$prof_dist_file_pat','$profile_region','$profile_region','$profile_region', 1, $prof_index, $dist_index, $t2c_index);exit;";

c=""; for i in 10 100 1000 10000; do c="$c $outdir/*nhood-$i[._]*.pdf"; done; echo $c;
/usr/bin/convert $c ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/OUT/`basename $outdir`.pdf;
