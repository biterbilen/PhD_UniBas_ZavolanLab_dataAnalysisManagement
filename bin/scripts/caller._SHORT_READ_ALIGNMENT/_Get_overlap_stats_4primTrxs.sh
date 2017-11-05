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
two=${10}
two2=${11}
~bilebi00/_PAPD5/scripts/innerJoin.pl ${t1}RNAseq-$c1-0${two}VS${t1}RNAseq-$c1-${tr1}Dis3L2$two.gene_exp.diff.$s1 ${t2}RNAseq-$c2-0${two2}VS${t2}RNAseq-$c2-${tr2}Dis3L2$two2.gene_exp.diff.$s2 1 1 3 | sort | uniq > $tag-$c1-$c2-$t1-$t2-$s1-$s2-$tr1-$tr2-$two-$two2; 
echo "#" $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 rep2:$two rep2:$two2: `less ${t1}RNAseq-$c1-0${two}VS${t1}RNAseq-$c1-${tr1}Dis3L2$two.gene_exp.diff.$s1 | cut -f 3 | sort | uniq | wc -l` `less ${t2}RNAseq-$c2-0${two2}VS${t2}RNAseq-$c2-${tr2}Dis3L2$two2.gene_exp.diff.$s2 | cut -f 3 | sort | uniq | wc -l` `less $tag-$c1-$c2-$t1-$t2-$s1-$s2-$tr1-$tr2-$two-$two2 | wc -l`;

return 0	
}

trusted=_trusted; #"";
outdir=Stepanka/OVERLAP_STATS$trusted; 
mkdir $outdir; pushd $outdir;

for d in ../CLIP_comparison$trusted/*diff; do 
	i=`basename $d`; 
	less $d | awk '$8>$9 && $14~/yes/{ print; }' > $i.pos; 
	less $d | awk '$8<$9 && $14~/yes/{ print; }' > $i.neg;  
	echo "#" $i `less $d | awk '$14~/yes/{ print; }' | wc -l` `less $i.pos | wc -l` `less $i.neg | wc -l`; 
done
#mRNAseq-HEK-0VSmRNAseq-HEK-oeDis3L2.gene_exp.diff 126 77 49
#mRNAseq-HEK-0VSmRNAseq-HEK-siDis3L2.gene_exp.diff 94 44 50
#totRNAseq-HEK-0_2VStotRNAseq-HEK-0.gene_exp.diff 146 69 77
#totRNAseq-HEK-0_2VStotRNAseq-HEK-oeDis3L2_2.gene_exp.diff 186 109 77
#totRNAseq-HEK-0_2VStotRNAseq-HEK-siDis3L2_2.gene_exp.diff 150 89 61
#totRNAseq-HEK-0VStotRNAseq-HEK-oeDis3L2.gene_exp.diff 165 91 74
#totRNAseq-HEK-0VStotRNAseq-HEK-siDis3L2.gene_exp.diff 172 98 74
#totRNAseq-HEK-oeDis3L2_2VStotRNAseq-HEK-oeDis3L2.gene_exp.diff 158 80 78
#totRNAseq-HEK-siDis3L2_2VStotRNAseq-HEK-siDis3L2.gene_exp.diff 155 88 67
#totRNAseq-Hela-0VStotRNAseq-Hela-siDis3L2.gene_exp.diff 126 74 52

echo "# Destabilized genes"
#TODO there should be more of these comparisons
c1=HEK; c2=Hela; t1=tot; t2=tot; s1=neg; s2=neg; tr1=si; tr2=si; two="_2"; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=Hela; t1=tot; t2=tot; s1=neg; s2=neg; tr1=si; tr2=si; two=""; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=neg; s2=neg; tr1=si; tr2=si; two="_2"; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=neg; s2=neg; tr1=si; tr2=si; two=""; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=pos; s2=pos; tr1=oe; tr2=oe; two=""; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=pos; s2=pos; tr1=oe; tr2=oe; two="_2"; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=tot; s1=neg; s2=pos; tr1=si; tr2=oe; two=""; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=tot; s1=neg; s2=pos; tr1=si; tr2=oe; two="_2"; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=m; t2=m; s1=neg; s2=pos; tr1=si; tr2=oe; two=""; two2="";
run destabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

less destabilized* | sort | uniq > all_destabilized

# Destabilized genes
# HEK Hela tot tot neg neg si si rep2:_2 rep2:: 61 52 4
# HEK Hela tot tot neg neg si si rep2: rep2:: 74 52 9
# HEK HEK tot m neg neg si si rep2:_2 rep2:: 61 50 2
# HEK HEK tot m neg neg si si rep2: rep2:: 74 50 5
# HEK HEK tot m pos pos oe oe rep2: rep2:: 91 77 7
# HEK HEK tot m pos pos oe oe rep2:_2 rep2:: 109 77 11
# HEK HEK tot tot neg pos si oe rep2: rep2:: 74 91 6
# HEK HEK tot tot neg pos si oe rep2:_2 rep2:: 61 91 4
# HEK HEK m m neg pos si oe rep2: rep2:: 50 77 5

#-----------------------
echo "# Stabilized genes"

c1=HEK; c2=Hela; t1=tot; t2=tot; s1=pos; s2=pos; tr1=si; tr2=si; two=""; two2="";
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=pos; s2=pos; tr1=si; tr2=si; two=""; two2="";
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=m; s1=neg; s2=neg; tr1=oe; tr2=oe; two=""; two2="";
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=tot; t2=tot; s1=pos; s2=neg; tr1=si; tr2=oe; two=""; two2="";
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

c1=HEK; c2=HEK; t1=m; t2=m; s1=pos; s2=neg; tr1=si; tr2=oe; two=""; two2="";
run stabilized $c1 $c2 $t1 $t2 $s1 $s2 $tr1 $tr2 $two $two2

# Stabilized genes
# HEK Hela tot tot pos pos si si 98 73 9
# HEK HEK tot m pos pos si si 98 44 5
# HEK HEK tot m neg neg oe oe 74 49 5
# HEK HEK tot tot pos neg si oe 98 74 5
# HEK HEK m m pos neg si oe 44 49 2

#TODO add rep 2 for HEK

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
