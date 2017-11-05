#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file_$TASK_ID
#$ -t 1-3

#!/bin/bash

#TODO get this from Chris for the new server 
outdir=SP

mkdir -p $outdir; pushd $outdir

annot=mRNA; #ALL; #mRNA
d=50 
N=10000
feat=T2C
l1s=(288 288 289)
l2s=(288 289 289)
ps=(EWSR1_1 EWSR1_12 EWSR1_2)

l1=${l1s[$((SGE_TASK_ID-1))]}
l2=${l2s[$((SGE_TASK_ID-1))]}
p=${ps[$((SGE_TASK_ID-1))]}

of=${p}_${annot}_$feat.$d.top$N.sp
#TODO ths doesn't work; ask chris later
~/progs/clipz/analysis/cluster/superCluster.sh $annot $d $l1:$feat:$N,$l2:$feat:$N $of

exit;
#~bilebi00/bin/java -jar ~rodak/topClusters.jar $annot $d $l1:$feat:$N,$l2:$feat:$N $of
~bilebi00/_KNOCKDOWNS/scripts/getCommonGroupId_fromSP.pl $of > ${of}_cls

source=superClusters
name=${p}_${annot}_$feat.$d.top$N.sp_cls
less $name | cut -f 11 | sort | uniq > $name.genenames;
less $name | perl -e '$src=shift; $name=shift; $_=<>; while(<>){@t=split; next if ($t[11] eq "NULL"); my (@s)=split/,/,$t[12]; print join("\t", $t[1],$src,$name,$t[3],$t[4],@s,$t[2],".",".",".",".",".",".",".",$t[10]),"\n"; }' $source $name | uniq > $name.data

rm $of
