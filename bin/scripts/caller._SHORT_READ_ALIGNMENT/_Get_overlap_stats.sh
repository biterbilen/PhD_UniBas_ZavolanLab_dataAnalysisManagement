#!/bin/bash

function run {
tag=$1
c1=$2
c2=$3
t1=$4
t2=$5
s1=$6
s2=$7
tr1=$8
tr2=$9
~bilebi00/_PAPD5/scripts/innerJoin.pl ${t1}RNAseq-$c1-0VS${t1}RNAseq-$c1-${tr1}Dis3L2.gene_exp.diff.$s1 ${t2}RNAseq-$c2-0VS${t2}RNAseq-$c2-${tr2}Dis3L2.gene_exp.diff.$s2 1 1 3 | sort | uniq > $tag-$c1-$c2-$t1-$t2-$s1-$s2-$tr1-$tr2; 
echo $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 `less ${t1}RNAseq-$c1-0VS${t1}RNAseq-$c1-${tr1}Dis3L2.gene_exp.diff.$s1 | cut -f 3 | sort | uniq | wc -l` `less ${t2}RNAseq-$c2-0VS${t2}RNAseq-$c2-${tr2}Dis3L2.gene_exp.diff.$s2 | cut -f 3 | sort | uniq | wc -l` `less $tag-$c1-$c2-$t1-$t2-$s1-$s2-$tr1-$tr2 | wc -l`;

return 0	
}


#TODO ln *diff from CLIP_comparison; run CLIP_comparison again


for i in *diff; do less $i | awk '$8>$9 && $14~/yes/{ print; }' > $i.pos; less $i | awk '$8<$9 && $14~/yes/{ print; }' > $i.neg;  echo $i `less $i | grep "yes$" | wc -l` `less $i.pos | wc -l` `less $i.neg | wc -l`; done
#mRNAseq-HEK-0VSmRNAseq-HEK-oeDis3L2.gene_exp.diff 1276 741 535
#mRNAseq-HEK-0VSmRNAseq-HEK-siDis3L2.gene_exp.diff 1183 528 655
#totRNAseq-HEK-0VStotRNAseq-HEK-oeDis3L2.gene_exp.diff 2121 1047 1074
#totRNAseq-HEK-0VStotRNAseq-HEK-siDis3L2.gene_exp.diff 2373 1149 1224
#totRNAseq-Hela-0VStotRNAseq-Hela-siDis3L2.gene_exp.diff 1930 1267 663


echo "Destabilized genes"

c1=HEK; c2=Hela; t1=tot; t2=tot; s1=neg; s2=neg; tr1=si; tr2=si;
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=neg; s2=neg; tr1=si; tr2=si;
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=pos; s2=pos; tr1=oe; tr2=oe;
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=tot; t2=tot; s1=neg; s2=pos; tr1=si; tr2=oe;
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=m; t2=m; s1=neg; s2=pos; tr1=si; tr2=oe;
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

#HEK Hela tot tot neg neg si si 413 277 34
#HEK HEK tot m neg neg si si 413 231 23
#HEK HEK tot m pos pos oe oe 369 309 27
#HEK HEK tot tot neg pos si oe 413 369 48
#HEK HEK m m neg pos si oe 231 309 18

less destabilized* | sort | uniq > all_destabilized


#-----------------------
echo "Stabilized genes"

c1=HEK; c2=Hela; t1=tot; t2=tot; s1=pos; s2=pos; tr1=si; tr2=si;
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=pos; s2=pos; tr1=si; tr2=si;
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=neg; s2=neg; tr1=oe; tr2=oe;
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=tot; t2=tot; s1=pos; s2=neg; tr1=si; tr2=oe;
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

c1=HEK; c2=HEK; t1=m; t2=m; s1=pos; s2=neg; tr1=si; tr2=oe;
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2

#HEK Hela tot tot pos pos si si 358 403 38
#HEK HEK tot m pos pos si si 358 207 14
#HEK HEK tot m neg neg oe oe 350 222 18
#HEK HEK tot tot pos neg si oe 358 350 29
#HEK HEK m m pos neg si oe 207 222 25

less stabilized* | sort | uniq > all_stabilized

#get maximum clip score
#FIXME very very slow
for sig_gene_file in all_destabilized all_stabilized; do
	echo $sig_gene_file;
	for clipf in ~bilebi00/_DIS3/4WEB_DIS3L2/PC*\.data; do
		outf=`basename $clipf .data`.$sig_gene_file;
		echo $outf.$sig_gene_file;
		less $sig_gene_file | while read i; do if [ $i == "-" ]; then continue; fi; less $clipf | grep -w $i | grep -w EXON | sort -k 7,7gr | head -n 1 ; done > $outf
	done
done

header="chr source name beg end score score2 str frm foreground foreground2 background background2 annot_type annot_ids annot_id"

echo $header | sed 's/ /\t/g' > PC.all_stabilized; 
cat PC_DIS3L2_*all_stabilized >> PC.all_stabilized

echo $header | sed 's/ /\t/g' > PC.all_destabilized; 
cat PC_DIS3L2_*all_destabilized >> PC.all_destabilized

#run a.R in the same folder

R --no-save --args PC.all_stabilized  < a.R
R --no-save --args PC.all_destabilized  < a.R
