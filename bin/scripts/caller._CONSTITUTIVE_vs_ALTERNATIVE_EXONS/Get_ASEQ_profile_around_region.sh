#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-libCount

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=output_directory
ptag=tag
libs=(libs)
libCuts=(libCuts)
profile=profile;
profile_type=profileType
top_event_density=topEventdensity
nhood=nhood

mkdir -p $outdir/$ptag; cd $outdir/$ptag;

dist_index=2; 
prof_index=8; 
t2c_index=17; 

i=$(($SGE_TASK_ID - 1));
lib=${libs[$i]};
libCut=${libCuts[$i]};

inp=$lib.inp;
stat=$lib.stat
cls=$lib.cls
bed=$lib.bed

#The format of the files should be like this
# $idl,$c,$s,$b,$e,$name,$dummy,$dummy2,$dumm3
if [ $profile_type == POLYA ]; then
	zless ~bilebi00/DATA/PolyA_DB_mammal_human_hg18.gz | awk '{OFS="\t";} { if (NR > 1 ) { print $5, $2, $7, $3+1, $4,"POLYA","dummy1","dummy2","dummy3";; } }' > $inp
elif [[ $profile_type =~ POLYA_ ]]; then
	if [[ $profile_type =~ _SINGLE ]]; then
		file=`echo $profile_type | sed -e 's/_SINGLE//'`;
		less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/$file | grep -v NULL | perl -e '$pr=shift; while(<>){ chomp; @t = split; if ($t[7] == 1) { splice(@t,5,0,$pr); splice(@t,6,1,$t[0]); $t[0]=$t[1].$t[2].$t[3]; print join("\t", @t), "\n"; } }' $profile_type > $inp;
	else
		less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/$profile_type | grep -v NULL | perl -e '$pr=shift; while(<>){ chomp; @t = split; splice(@t,5,0,$pr); splice(@t,6,1,$t[0]); $t[0]=$t[1].$t[2].$t[3]; print join("\t", @t), "\n"; }' $profile_type > $inp;
	fi
elif [ $profile_type == "TIA_ENHANCED_CE" ]; then
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/konig_wang_2010_iclip_papers/splicing-change_TIA_hnRNPC.unix | perl -e '$pr=shift; while(<>){chomp; @t=split/\t/; @a=split/[:\-]/,$t[4];  if ($t[15]>1 and $t[13] eq "CE") { print join "\t", ($t[3]."_".$t[12],$a[0],$t[6], $a[1], $a[2], $pr, $t[15], "dummy2", "dummy3"),"\n";  }}' $profile_type > $inp;
elif [ $profile_type == "TIA_SILENCED_CE" ]; then
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/konig_wang_2010_iclip_papers/splicing-change_TIA_hnRNPC.unix | perl -e '$pr=shift; while(<>){chomp; @t=split/\t/; @a=split/[:\-]/,$t[4];  if ($t[15]<-1 and $t[13] eq "CE") { print join "\t", ($t[3]."_".$t[12],$a[0],$t[6], $a[1], $a[2], $pr, $t[15], "dummy2", "dummy3"),"\n";  }}' $profile_type > $inp;
elif [ $profile_type == "hnRNPC_ENHANCED_CE" ]; then
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/konig_wang_2010_iclip_papers/splicing-change_TIA_hnRNPC.unix | perl -e '$pr=shift; while(<>){chomp; @t=split/\t/; @a=split/[:\-]/,$t[4];  if ($t[14]>1 and $t[13] eq "CE") { print join "\t", ($t[3]."_".$t[12],$a[0],$t[6], $a[1], $a[2], $pr, $t[14], "dummy2", "dummy3"),"\n";  }}' $profile_type > $inp;
elif [ $profile_type == "hnRNPC_SILENCED_CE" ]; then
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/konig_wang_2010_iclip_papers/splicing-change_TIA_hnRNPC.unix | perl -e '$pr=shift; while(<>){chomp; @t=split/\t/; @a=split/[:\-]/,$t[4];  if ($t[14]<-1 and $t[13] eq "CE") { print join "\t", ($t[3]."_".$t[12],$a[0],$t[6], $a[1], $a[2], $pr, $t[14], "dummy2", "dummy3"),"\n";  }}' $profile_type > $inp;
elif [ $profile_type == "TSR" ]; then 
	less ~bilebi00/_AGO1/TSR/human_TSR_v1.gff.cls_format.w_gene_symbol.w_event_count > $inp 
elif [ $profile_type == "TSR_SINGLE" ]; then 
	less ~bilebi00/_AGO1/TSR/human_TSR_v1.gff.cls_format.w_gene_symbol.w_event_count | perl -e '$pr=shift;while(<>){chomp;@t=split; if($t[8]==1){$t[5]=$pr; print join("\t",@t),"\n"; } }' $profile_type > $inp 
elif [[ $profile_type =~ STRICT ]]; then
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/EXONS_ExpressedStrictBoundary/exons.strict_w_genesymbols.expressed_genes | grep $profile_type | awk 'BEGIN{OFS="\t"}{ print $0, "dummy2", "dummy3"}' > $inp
elif [[ $profile_type =~ AGO ]]; then
	lc=`grep $profile_type ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/protein_file_libcut | cut -f 3`;
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/Reproducibility_v05/$profile_type | perl -e ' $pr=shift; $lc=shift; while(<>){ @t=split; next if($t[6]<$lc); $b=2; $e=9; if($t[1] == "-") { $b=9; $e=2;} $t[3]-=$b; $t[4]+=$e; $t[5]=$pr; splice(@t,0,0,$t[7]); splice(@t,3,1); $t[6] = $t[0]; $t[0] = $t[1].$t[2].$t[3]; splice(@t, 7); print join("\t", @t, "dummy2", "dummy3", "\n"); } ' $profile_type $lc > $inp
else #other protein profiles should exist in the protein_file_libcut
	lc=`grep $profile_type ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/protein_file_libcut | cut -f 3`;
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/Reproducibility_v05/$profile_type | perl -e ' $pr=shift; $lc=shift; while(<>){ @t=split; next if($t[6]<$lc); $b=0; $e=0; if($t[1] == "-") { $b=0; $e=0;} $t[3]-=$b; $t[4]+=$e; $t[5]=$pr; splice(@t,0,0,$t[7]); splice(@t,3,1); $t[6] = $t[0]; $t[0] = $t[1].$t[2].$t[3]; splice(@t, 7); print join("\t", @t, "dummy2", "dummy3", "\n"); } ' $profile_type $lc > $inp
fi

#add densities to the input file for clustering  and some statistics in the stat file
echo -e "totalNumberOfRegions\t"`wc -l $inp | awk '{ print $1}'` > $stat;
echo -e "nhood\t$nhood" >> $stat

#GAPDH chr12 + 6517785 6517785 76871.9526  POLYA_217 1
#NULL  chrM  - 228 228 44532.0702  POLYA_217 11606
less ../$lib | sort -k 6,6gr | head -n $top_event_density | perl -e '$p=shift; $libc=shift; @a=(); while(<>){@t=split;  push(@a, $t[5]); $t[8] = $t[0]; $t[0] = $p; print join ("\t",@t)."\n"; } print STDERR "medianEventSum\t$a[int($#a/2)]\n", "totalNumberOfEvents\t",scalar @a,"\n"; ' $profile $libCut >> $inp 2>> $stat;

prof_file=$lib.profile;
prof_dist_file=${prof_file}_distribution

#------------------------------------
#main job; clustering and checking the boundary

#cluster region and Event peaks
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls;

#get profile of event around region
echo "~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nhood $cls $profile_type $profile 1 0 $prof_file 2> $prof_dist_file.tmp";
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nhood $cls $profile_type $profile 1 0 $prof_file 2> $prof_dist_file.tmp;
grep "^Distance" $prof_dist_file.tmp > $prof_dist_file;

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/profile2bed.pl $prof_dist_file > $bed;

#TODO remove $prof_dist_file as well when done with the bug in profile2bed.pl script
echo -e "totalNumberofHitRegions\t"`wc -l $bed | awk '{print $1; }'` >> $stat;

rm $prof_dist_file.tmp $inp $cls;
exit;

