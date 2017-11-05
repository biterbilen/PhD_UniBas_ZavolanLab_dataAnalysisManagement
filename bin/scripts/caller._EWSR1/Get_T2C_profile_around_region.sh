#!/bin/bash
#$ -S /bin/bash
#$ -q queue_type@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-libCount

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

outdir=output_directory
ptag=tag
libs=(libs)
libCuts=(libCuts)
gblibs=(gblibs)
profile=profile;
profile_type=profileType
top_T2C_density=topT2Cdensity
nhood=nhood
prot_libs_libcuts_file=protein_libs_libcuts_file

pushd $outdir;

mkdir -p $ptag; pushd $ptag;

dist_index=2; 
prof_index=8; 
t2c_index=17; 

i=$(($SGE_TASK_ID - 1));
lib=${libs[$i]};
gblib=${gblibs[$i]};
libCut=${libCuts[$i]};

inp=$lib.inp;
stat=$lib.stat
cls=$lib.cls
gtf=$lib.gtf

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
	lc=`grep $profile_type $prot_libs_libcuts_file | cut -f 3`;
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/Reproducibility_v05/$profile_type | perl -e ' $pr=shift; $lc=shift; while(<>){ @t=split; next if($t[6]<$lc); $b=2; $e=9; if($t[1] == "-") { $b=9; $e=2;} $t[3]-=$b; $t[4]+=$e; $t[5]=$pr; splice(@t,0,0,$t[7]); splice(@t,3,1); $t[6] = $t[0]; $t[0] = $t[1].$t[2].$t[3]; splice(@t, 7); print join("\t", @t, "dummy2", "dummy3", "\n"); } ' $profile_type $lc > $inp
else #other protein profiles should exist in the protein_file_libcut
	lc=`grep $profile_type $prot_libs_libcuts_file | cut -f 3`;
	less ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/Reproducibility_v05/$profile_type | perl -e ' $pr=shift; $lc=shift; while(<>){ @t=split; next if($t[6]<$lc); $b=0; $e=0; if($t[1] == "-") { $b=0; $e=0;} $t[3]-=$b; $t[4]+=$e; $t[5]=$pr; splice(@t,0,0,$t[7]); splice(@t,3,1); $t[6] = $t[0]; $t[0] = $t[1].$t[2].$t[3]; splice(@t, 7); print join("\t", @t, "dummy2", "dummy3", "\n"); } ' $profile_type $lc > $inp
fi

#add densities to the input file for clustering  and some statistics in the stat file
echo -e "totalNumberOfRegions\t"`wc -l $inp | awk '{ print $1}'` > $stat;
echo -e "nhood\t$nhood" >> $stat
less ../$lib | sort -k 8,8gr | head -n $top_T2C_density | perl -e '$p=shift; $libc=shift; @a=(); while(<>){@t=split; if ($t[6] >= $libc) { push(@a, $t[$#t]); $pos=$t[3]+($t[4]-$t[3])/2; print join ("\t",$p,$t[0],$t[1],$pos,$pos,$t[7])."\n";}} print STDERR "medianEventSum\t$a[int($#a/2)]\n", "totalNumberOfEvents\t",scalar @a,"\n"; ' $profile $libCut >> $inp 2>> $stat;

prof_file=$lib.profile;
prof_dist_file=${prof_file}_distribution

#cluster region and T2C peaks
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls;

#get profile of T2C around region
echo "~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nhood $cls $profile_type $profile 1 0 $prof_file 2> $prof_dist_file.tmp";
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nhood $cls $profile_type $profile 1 0 $prof_file 2> $prof_dist_file.tmp;
grep "^Distance" $prof_dist_file.tmp > $prof_dist_file;

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/profile2gtf.pl $prof_dist_file $nhood "$gblib" "$ptag $lib" "$ptag $lib" > $gtf

echo -e "totalNumberofHitRegions\t"`grep "^chr" $gtf | wc -l | awk '{print $1; }'` >> $stat;


rm $prof_dist_file* $inp $cls;

echo DONE

exit;

