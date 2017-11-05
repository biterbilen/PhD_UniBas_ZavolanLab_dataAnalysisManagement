#!/bin/bash
#$ -S /bin/bash
#$ -q long_m8
#$ -pe MT8 4
#$ -e stderr_$TASK_ID
#$ -o stdout_$TASK_ID
#$ -j n
#$ -cwd
#$ -l mem_total=6000M
#$ -t 1-2

soutdir=UniqueRawData_ALL_woContamination
mkdir -p $soutdir; pushd $soutdir;

#run it max 3 at a time; main script submits to clusters
#TODO later on all not only mRNA
libs=(1341 1342 240 306 490 403 492 498 499 314 315 1037 1039 984 308 309 299 300 211 213 445 446 312 313 288 289 271 272 502 503 238 239 294 295 486 478 479 482 483 480 481 484 485)
SGE_TASK_ID=$(($SGE_TASK_ID - 1));
id=${libs[$SGE_TASK_ID]}
echo $SGE_TASK_ID $id;

dir=/import/wnz/home/mirz/progs/clipz/administration/data/density
odir=DB_all_$id
#annot=mRNA
annot=ALL
uniq=1
#size=`cat ~mirz/clipzServer/data/database/$id/sample_size`;
if [ -e $odir ]; then
	echo $odir exists;
	exit;
fi
#-----------------------------------------------
#density
mapped_sequences=~mirz/clipzServer/data/database/$id/t_mapped_sequence.gz
if [ ! -e $mapped_sequences ]; then
	mapped_sequences=t_mapped_sequence$id
	if [ ! -e $mapped_sequences ]; then
		echo Database dump for t_mapped_sequence_$id
		/import/bc2/soft/bin/mysql -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped_sequences <<END
        select * from t_mapped_sequence_$id where genome_count_total=1;
        quit
END
		less $mapped_sequences | grep -v "^id" | sed 's/NULL/\\N/g' > tmp$mapped_sequences
		mv tmp$mapped_sequences $mapped_sequences
	else
		echo Database dump done before
	fi
fi


~/BITER/scripts/selectAnnot.pl $mapped_sequences ~mirz/clipzServer/data/cluster/$id/genome.sorted $annot $uniq > $id.genome.sorted 2> $id.effective_sample_size 

#size=`cat $id.sample_size`; #normalized copy numbers #size for filtered
size=1000000 #to get the raw copy numbers
$dir/density.sh $mapped_sequences $id.genome.sorted $size $odir

tar -zcvf $odir.tar.gz $odir
#-----------------------------------------------
#clean

rm -rf $odir $id.genome.sorted
echo $odir DONE

