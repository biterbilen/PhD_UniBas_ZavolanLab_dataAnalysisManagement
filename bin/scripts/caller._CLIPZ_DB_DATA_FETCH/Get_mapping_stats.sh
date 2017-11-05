#!/bin/bash

soutdir=Summary
mkdir -p $soutdir; pushd $soutdir;

libs=(288 289 299 300 211 213 445 446 312 313 271 272 502 503 238 239 294 295 486 478 479 482 483 480 481 484 485)
libs=(1341 1342 238 239 240 288 289 308 309)
libs=(1055)

for id in ${libs[@]}; do 
	echo $id;
	mysql -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $id.summary <<END
	select * from t_sample_statistics where t_sample_id=$id;
	select reads_c,extracted_c,extracted_p,filtered_c,filtered_p,unique_c,unique_p,mapped_c,mapped_p from t_annotation_sample_stats where t_sample_id=$id;
	quit
END
done

#mapped_c corresponds to sum(copies) from t_mapped_sequence_$id
#select count(*),sum(copies) from t_mapped_sequence_288 where genome_count_total IS NOT NULL;

