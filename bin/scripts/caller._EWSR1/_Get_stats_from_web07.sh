#!/bin/bash

outdir=Analysis

pushd $outdir

soutdir=Summary
mkdir $soutdir; pushd $soutdir;

#expressed

#clipped by different regions
#scp mirz@bc2-web07.bc2.unibas.ch:~/BITER/28*.summary .

for id in 288 289; do
	typename=mutation
	type=7
	total=`less $id.summary | awk 'BEGIN{c=0;} {if ($3 == g) { c = c+$5; print c; }}' g=$type |tail -n 1`
	less $id.summary | awk 'BEGIN{OFS="\t"; print "lib","type","value";} { if($3==g) print $2, $4, 100*$5/t; }' t=$total g=$type > $id.$typename

	typename=annotation
	type=1
	less $id.summary | awk 'BEGIN{OFS="\t"; print "lib","type","value";} { if($3==g) print $2, $4, 100*$5; }' g=$type | sort -k 3,3r > $id.$typename

	typename=mRNA
	type=2
	total=`less $id.summary | awk 'BEGIN{c=0;} {if ($3 == g) { c = c+$5; print c; }}' g=$type |tail -n 1`
	less $id.summary | awk 'BEGIN{OFS="\t"; print "lib","type","value";} { if($3==g) print $2, $4, 100*$5/t; }' t=$total g=$type > $id.$typename
done

for i in mutation annotation mRNA; do
	cat 288.$i > all.$i ; awk 'NR>1{print;}' 289.$i >> all.$i;
	R --no-save --args all.$i < ~bilebi00/_EWSR1/scripts/pie.R > /dev/null
done



#
#RSA keys are set for login wo password after ssh-keygen
#scp  ~/.ssh/id_rsa.pub mirz@bc2-web07.bc2.unibas.ch:~/.ssh/authorized_keys4
#
#pushd BITER
#for id in 288 289; do 
#	./stats.sh $id > $id.summary 
#done
#	
#stats.sh in web07
##!/bin/bash
#
#id=$1;
##id=288;
#mysql -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A <<!!
# select * from t_sample_statistics where t_sample_id=$id;
#	quit
#	!!

exit;

echo DONE
