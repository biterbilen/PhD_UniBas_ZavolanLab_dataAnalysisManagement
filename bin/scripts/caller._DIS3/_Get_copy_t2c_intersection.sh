#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j n
#$ -cwd
# $ -e LOG._Get_copy_t2c_intersection_PAPD5.sh 
#$ -e LOG._Get_copy_copy_intersection_PAPD5.sh 
# $ -o PAPD5.copy_t2c.intersection.txt
#$ -o PAPD5.copy_copy.intersection.txt

export PERL5LIB=$HOME/lib:$PERL5LIB

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> LOG._Get_copy_t2c_intersection.sh 
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

#!/bin/bash

#cat copy_t2c.txt | grep "^\-*3 \-*1" 

rand=$RANDOM
#1. copy t2c
copydir=~bilebi00/_DIS3/Normalization_PAPD5/
t2cdir=~bilebi00/_DIS3/Normalization_PAPD5_T2C/

t2c_zvalues=(1  1.5 1.5 2 2 1.5);
copy_zvalues=(3 4   5   5 6 6);

#2. copy copy 
copy_zvalues=(3 4 5 6);
for i in ${copy_zvalues[@]}; do
	less $copydir/tableall.norm_and_raw.04zscore.woNA | awk '{ if ($6>cut) {print $0;} }' cut=$i > 04.$i.copy$rand
	less $copydir/tableall.norm_and_raw.14zscore.woNA | awk '{ if ($6>cut) {print $0;} }' cut=$i > 14.$i.copy$rand
	less $copydir/_04_14.notfiltered | awk '{ if ($4>cut && $7>cut) {print $0;} }' cut=$i > 04_14.$i.copy$rand
	less $copydir/tableall.norm_and_raw.04zscore.woNA | awk '{ if ($6<cut) {print $0;} }' cut=-$i > 04.-$i.copy$rand
	less $copydir/tableall.norm_and_raw.14zscore.woNA | awk '{ if ($6<cut) {print $0;} }' cut=-$i > 14.-$i.copy$rand
	less $copydir/_04_14.notfiltered | awk '{ if ($4<cut && $7<cut) {print $0;} }' cut=-$i > 04_14.-$i.copy$rand
done
for i in ${t2c_zvalues[@]}; do
	less $t2cdir/tableall.norm_and_raw.04zscore.woNA | awk '{ if ($6>cut) {print $0;} }' cut=$i > 04.$i.t2c$rand
	less $t2cdir/tableall.norm_and_raw.14zscore.woNA | awk '{ if ($6>cut) {print $0;} }' cut=$i > 14.$i.t2c$rand
	less $t2cdir/_04_14.notfiltered | awk '{ if ($4>cut && $7>cut) {print $0;} }' cut=$i > 04_14.$i.t2c$rand
	less $t2cdir/tableall.norm_and_raw.04zscore.woNA | awk '{ if ($6<cut) {print $0;} }' cut=-$i > 04.-$i.t2c$rand
	less $t2cdir/tableall.norm_and_raw.14zscore.woNA | awk '{ if ($6<cut) {print $0;} }' cut=-$i > 14.-$i.t2c$rand
	less $t2cdir/_04_14.notfiltered | awk '{ if ($4<cut && $7<cut) {print $0;} }' cut=-$i > 04_14.-$i.t2c$rand
done

#2.
echo minzScore.copy minzScore.copy lib count.inters count.copy count.copy ...
for i in `seq 0 $((${#copy_zvalues[@]} - 1))`; do
	ccopy=${copy_zvalues[$i]};
	l=04
	l2=14
	echo -en $ccopy $ccopy"\t";
	echo -en "${l}_${l2}:" `~bilebi00/_PAPD5/scripts/innerJoin.pl $l.$ccopy.copy$rand $l2.$ccopy.copy$rand 1 1 '' | wc -l` `less $l.$ccopy.copy$rand | wc -l` `less $l2.$ccopy.copy$rand | wc -l` "\t" 
	echo 

	echo -en -$ccopy -$ccopy"\t";
	echo -en "${l}_${l2}:" `~bilebi00/_PAPD5/scripts/innerJoin.pl $l.-$ccopy.copy$rand $l2.-$ccopy.copy$rand 1 1 '' | wc -l` `less $l.-$ccopy.copy$rand | wc -l` `less $l2.-$ccopy.copy$rand | wc -l` "\t" 
	echo 
done
rm *copy$rand *t2c$rand;

exit;

#1.
echo minzScore.copy minzScore.t2c lib count.inters count.copy count.t2c
for i in `seq 0 $((${#copy_zvalues[@]} - 1))`; do
	ct2c=${t2c_zvalues[$i]};
	ccopy=${copy_zvalues[$i]};
	echo -en $ccopy $ct2c "\t";
	for l in 04 14 04_14; do 
		echo -en "$l:" `~bilebi00/_PAPD5/scripts/innerJoin.pl $l.$ccopy.copy$rand $l.$ct2c.t2c$rand 1 1 '' | wc -l` `less $l.$ccopy.copy$rand | wc -l` `less $l.$ct2c.t2c$rand | wc -l` "\t" 
	done
	echo 
done

for i in `seq 0 $((${#copy_zvalues[@]} - 1))`; do
	ct2c=-${t2c_zvalues[$i]};
	ccopy=-${copy_zvalues[$i]};
	echo -en $ccopy $ct2c "\t";
	for l in 04 14 04_14; do 
		echo -en "$l:" `~bilebi00/_PAPD5/scripts/innerJoin.pl $l.$ccopy.copy$rand $l.$ct2c.t2c$rand 1 1 '' | wc -l` `less $l.$ccopy.copy$rand | wc -l` `less $l.$ct2c.t2c$rand | wc -l` "\t" 
	done
	echo 
done

rm *copy$rand *t2c$rand;

