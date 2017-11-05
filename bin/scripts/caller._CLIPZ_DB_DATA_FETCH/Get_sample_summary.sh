#!/bin/bash

soutdir=UniqueRawData_ALL_woContamination
soutdir=Summary
mkdir -p $soutdir; pushd $soutdir;

#run it max 3 at a time; main script submits to clusters
libs=(1617 1618 1619 1620 445 446 308 309 1344 1218 1219 1341 1342 240 288 289 984 238 239 502 503 482 483 484 481 478 479 480 485 294 295 240 1037 1039 271 272 647 651)
libs=(1341 1342 238 239 240 288 289 308 309)
libs=(1055)
libs=(271 272 651 647 1037 1039)

SGE_TASK_ID=$(($1 - 1));
id=${libs[$SGE_TASK_ID]}
echo $SGE_TASK_ID $id;

dir=/import/wnz/home/mirz/progs/clipz/administration/data/density
#size=`cat ~mirz/clipzServer/data/database/$id/sample_size`;
#-----------------------------------------------
#extraction and mapping statistics
mysql --quick -u clipzRead -h wnzmsql01 -D clipz --password=readClipz -A > $id.general_summary <<END
	select reads_c,extracted_c,mapped_c,mapping_rate from t_annotation_sample_stats where t_sample_id=$id;
	quit
END
exit;
#-----------------------------------------------
#annotation summary for the selected -non contamination and unique mappers
mysql -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $id.annotation_summary <<END
	select annotation,sum(copies),count(*) from t_sequence_genome_alignment_$id,t_mapped_sequence_$id where t_mapped_sequence_$id.id=t_sequence_genome_alignment_$id.t_sequence_id && genome_count_total=1 && annotation!="bacterial" && annotation!="fungus" && annotation!="vector" && annotation!="viral" && annotation!="Markers_Adaptors" GROUP BY annotation;
	quit
END
#-----------------------------------------------
