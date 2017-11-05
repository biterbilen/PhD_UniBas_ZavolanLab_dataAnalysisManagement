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
sid_prot_f=~bilebi00/_DIS3L2/data/protein_sampleid_list
outdir=Analysis
pushd $outdir

soutdir=Summary
mkdir -p $soutdir; pushd $soutdir;

libs=(`less $sid_prot_f | cut -f 2`)

for lib in ${libs[@]}; do 
	scp mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/$lib.summary .
done

#expressed

#clipped by different regions
#scp mirz@bc2-web07.bc2.unibas.ch:~/BITER/28*.summary .

topN=100; #a big number not to filter any
typenames=(dummy annotation mRNA pre-mRNA dummy dummy dummy mutation)
type=1
for type in `eval echo {0..$((${#typenames[@]}-1))}`; do 
	typename=${typenames[$type]};
	if [ $typename == "dummy" ]; then continue; fi
	echo $typename
	rm -rf all.$typename*
	echo "lib type value" | awk 'BEGIN{OFS="\t";}{ print $1, $2, $3; }' > all.$typename
	for i in `eval echo {0..$((${#libs[@]}-1))}`; do
		id=${libs[$i]}
		idn=`grep -w $id $sid_prot_f | cut -f 3`;
		prep4R $id $topN $typename $type "$idn" >> all.$typename;
	done
	R --no-save --args all.$typename < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
	pwd=`pwd`;
	pushd ~bilebi00/www/Summary/ARE
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

exit;

echo DONE
