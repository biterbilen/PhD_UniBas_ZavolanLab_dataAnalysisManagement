#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@high_mem_node
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file

outdir=output_directory
nh=nhood

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

mkdir -p $outdir; cd $outdir;

cls=T2C_exon_cls.cls
strict=T2C_exon_cls.strict

echo Total Number of Strict Constitutive Exon Intervals: `less $strict | grep STRICT_CONSTITUTIVE_EXON | wc -l`;
echo Total Number of Strict Alternative Exon Intervals: `less $strict | grep STRICT_ALTERNATIVE_EXON | wc -l`;
echo Total Number of Strict Other Exon Intervals: `less $strict | grep STRICT_OTHER_EXON | wc -l`;
echo exon unification: `less $strict | grep EXON | wc -l` exons exist in $strict;

#Don't split this
#select T2C peaks around EXONS
echo ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $strict STRICT_ALTERNATIVE_EXON SS 1;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $strict STRICT_ALTERNATIVE_EXON SS 1 0;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $strict STRICT_CONSTITUTIVE_EXON SS 1 0;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $strict STRICT_OTHER_EXON SS 0 0;

##select T2C peaks around EXONS
#~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $cls 'AK294343.*ALTERNATIVE_EXON\(2' SS 0;
#~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $cls ALTERNATIVE_EXON SS 0;
#~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $cls CONSTITUTIVE_EXON SS 0;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $cls SINGLE_EXON SS 0 0;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $cls FIRST_EXON SS 0 0;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/check_profile.pl $nh $cls LAST_EXON SS 0 0;
