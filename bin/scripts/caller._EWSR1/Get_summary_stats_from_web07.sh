#!/bin/bash

function prep4R {
	id=$1
	topN=$2
	typename=$3
	type=$4
	idname=$5;
	total=`less $id.summary | awk -F "\t" 'BEGIN{c=0;} {if ($3 == g) { c = c+$5; print c; }}' g=$type | tail -n 1`
	less $id.summary | sort -k5,5gr | \
		awk -F "\t" 'BEGIN{OFS="\t";ft=0;i=0;} { if($3==g) {i=i+1; if(i<topN) {f=100*$5/t;ft=ft+f;print n,$4,f;} else if(i==topN) { print n,"other",100-ft;}}}' \
	 	t=$total g=$type topN=$topN n="$idname"
	return 0;
}
outdir=Analysis
pushd $outdir

soutdir=Summary
mkdir -p $soutdir; pushd $soutdir;

ids=(288 289)
idnames=("EWSR1 (80kDa)" "EWSR1 (120kDa)")
#expressed

#clipped by different regions
rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/Summary/28*.summary .

topN=100; #a big number not to filter any
typenames=(dummy annotation mRNA pre-mRNA dummy dummy dummy mutation)
for type in `eval echo {0..$((${#typenames[@]}-1))}`; do 
	typename=${typenames[$type]};
	if [ $typename == "dummy" ]; then continue; fi
	echo $typename
	rm -rf all.$typename*
	echo "lib type value" | awk 'BEGIN{OFS="\t";}{ print $1, $2, $3; }' > all.$typename
	for i in `eval echo {0..$((${#ids[@]}-1))}`; do
		id=${ids[$i]}
		idn=${idnames[$i]}
		prep4R $id $topN $typename $type "$idn" >> all.$typename;
		#prep4R $id $topN $typename $type "$idn" | perl -e 'while(<>){@t=split/\t/; print join("\t",$t[0], sprintf("%s (%.2f \%)", $t[1],$t[2]), $t[2]), "\n"; }' >> all.$typename
	done
	#R --no-save --args all.$typename < ~bilebi00/_EWSR1/scripts/pie.R > /dev/null
	R --no-save --args all.$typename < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
	pwd=`pwd`;
	pushd ~bilebi00/www/Summary/EWSR1
	ln -sf $pwd/all.$typename.barchart.pdf .
	popd
done
exit;



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
