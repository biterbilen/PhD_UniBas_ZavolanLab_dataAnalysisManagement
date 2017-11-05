#!/bin/bash
#$ -S /bin/bash
#$ -q long_m8
#$ -pe MT8 4
#$ -e stderr_$TASK_ID
#$ -o stdout_$TASK_ID
#$ -j n
#$ -cwd
#$ -l mem_total=6000M
#$ -t 1-10

#run it max 3 at a time; main script submits to clusters
libs=(984 1341 1342 240 288 289 314 315)
SGE_TASK_ID=$1
seqLenCut=17

if [ $SGE_TASK_ID -gt ${#libs[@]} ]; then
	echo "$SGE_TASK_ID > ${#libs[@]}"
	exit;
fi

#make directory
soutdir=UniqueRawData_mRNA_repeat_none_seqLenCut$seqLenCut
mkdir -p $soutdir; pushd $soutdir;

#clusters or the 
SGE_TASK_ID=$(($SGE_TASK_ID - 1));
id=${libs[$SGE_TASK_ID]}
echo $SGE_TASK_ID $id;

dir=/import/wnz/home/mirz/progs/clipz/administration/data/density
odir=DB_all_$id
#annot=mRNA
annot=mRNA,repeat,none
uniq=1
#size=`cat ~mirz/clipzServer/data/database/$id/sample_size`;
if [ -e $odir ]; then
	echo $odir exists;
	exit;
fi
#-----------------------------------------------
#density
mapped_sequences=~mirz/clipzServer/data/database/$id/t_mapped_sequence.gz
#2012.08.08 always get it from the server
mapped_sequences=t_mapped_sequence$id
if [ ! -e $mapped_sequences ]; then
	#get unique
	if [ ! -e $mapped_sequences ]; then
		echo Database dump for t_mapped_sequence_$id
		/import/bc2/soft/bin/mysql -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped_sequences <<END
        select * from t_mapped_sequence_$id where genome_count_total=1 and seq_len>$seqLenCut;
        quit
	exit;
END
		less $mapped_sequences | grep -v "^id" | sed 's/NULL/\\N/g' > tmp$mapped_sequences
		mv tmp$mapped_sequences $mapped_sequences
	else
		echo Database dump done before
	fi
fi

version=0
~/BITER/scripts/selectAnnot.pl $mapped_sequences ~mirz/clipzServer/data/cluster/$id/genome.sorted $annot $uniq $version > $id.genome.sorted 2> $id.effective_sample_size 
#~/BITER/scripts/selectAnnot.pl $mapped_sequences ~mirz/clipzServer/data/cluster/$id/genome.sorted $annot $uniq $version 1 > $id.other.genome.sorted 2>> $id.effective_sample_size 

#size=`cat $id.sample_size`; #normalized copy numbers #size for filtered
size=1000000 #to get the raw copy numbers
$dir/density.sh $mapped_sequences $id.genome.sorted $size $odir
#$dir/density.sh $mapped_sequences $id.other.genome.sorted $size other$odir

tar -zcvf $odir.tar.gz $odir
#tar -zcvf other$odir.tar.gz other$odir

rm -rf $odir # $id.genome.sorted
rm -rf other$odir # $id.genome.sorted

echo DONE

