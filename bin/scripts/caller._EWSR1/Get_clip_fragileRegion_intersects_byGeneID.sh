#!/bin/bash

rd=~bilebi00/_EWSR1/data/ReplicationDomain/;
rybas=(${rd}GE_BG01_ESC_All.txt ${rd}GE_BG01_NPC_All.txt ${rd}GE_BG02_ESC_All.txt ${rd}GE_H7_ESC_All.txt)

chimerdiff=~bilebi00/_EWSR1/data/ChimerDB/mRNA_Information_Diff.txt
chimersame=~bilebi00/_EWSR1/data/ChimerDB/mRNA_Information_Same.txt
gif=~bilebi00/_EWSR1/data/hg19_TR.info.all
expressed=~bilebi00/_EWSR1/data/clipz_RNAseq/RNAseq.raw.qn.w_geneSymbol
outdir=Analysis

pushd $outdir

soutdir=FragileRegion
mkdir -p $soutdir; pushd $soutdir;

#expressed
#~bilebi00/_DIS3/scripts/assign_gene_id.pl $expressed $gif 0 NA | awk ' BEGIN{OFS="\t";} { if (NR>1) { print $NF; } }' > expressed
less $expressed | awk ' BEGIN{OFS="\t";} { if (NR>1) { print $NF; } }' > expressed

#clipped
#pools gene_ids mapping to the same region
less ../EWSR1.annot | awk '$1 == "mRNA" && $11>2 && $9 > 0.5 && $10==2{ print;}' | cut -f 3 | perl -e 'while(<>){chomp; @t=split/\|/; for(@t){print $_,"\n";}}' | sort | uniq | grep -v -w NA > clipped 

#Replication domain
#Ryba T, Hiratani I, Lu J, Itoh M et al. Evolutionarily conserved replication timing profiles predict long-range chromatin interactions and distinguish closely related cell types. Genome Res 2010 Jun;20(6):761-70.
echo exp type pval > ryba.test_results 
for ryba in ${rybas[@]}; do
	N=`less $ryba | awk -F "\t" 'NF==7 && $7==1 {print $0;}' | cut -f 2 | sort | uniq | wc -l`;
	frac=`echo "$N / 20 " | bc`; #5% of the top
	of=ryba_`basename $ryba _All.txt`;

	#early
	f=$of.early;
	less $ryba | awk -F "\t" 'NF==7 && $7==1 {print $0;}' | sort -k 6,6gr | cut -f 2 | head -n $frac | sort | uniq > $f;
	grep -w -f $f expressed | sort | uniq > $f.expressed
	f=$f.expressed
	~bilebi00/_PAPD5/scripts/innerJoin.pl clipped.expressed $f 1 1 1 > $f-clipped
#	echo -n early $f-clipped clipped.expressed $f expressed;
	counts=`wc -l $f-clipped clipped.expressed $f expressed | awk 'NR<5{print $1}'`;
	p=`~bilebi00/bin/R --no-save --args $counts < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^TEST.RESULTS" -A 5 | tail -n 1` 
	echo $of EARLY $p

	#late
	f=$of.late; 
	less $ryba | awk -F "\t" 'NF==7 && $7==1 {print $0;}' | sort -k 6,6g | cut -f 2 | head -n $frac | sort | uniq > $f;
	grep -w -f $f expressed | sort | uniq > $f.expressed
	f=$f.expressed
	~bilebi00/_PAPD5/scripts/innerJoin.pl clipped.expressed $f 1 1 1 > $f-clipped
#	echo -n late clipped-x clipped.expressed x expressed;
	counts=`wc -l $f-clipped clipped.expressed $f expressed | awk 'NR<5{print $1}'`;
	p=`~bilebi00/bin/R --no-save --args $counts < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^TEST.RESULTS" -A 5 | tail -n 1` 
	echo $of LATE $p

	#random
	for i in `seq 1 1000`; do 
		f=$of.random
		less $ryba | perl -e 'while(<>){ chomp; @t=split/\t/; next if (@t !=7 or $t[6] != 1); $r=rand(); print $t[1],"\n" if ($r>0.5); }' | head -n $frac | sort | uniq > $f;
		grep -w -f $f expressed | sort | uniq > $f.expressed
		f=$f.expressed
		~bilebi00/_PAPD5/scripts/innerJoin.pl clipped.expressed $f 1 1 1 > $f-clipped
#		echo -n random clipped-x clipped.expressed x expressed;
		counts=`wc -l clipped-$f clipped.expressed $f expressed | awk 'NR<5{print $1}'`;
		p=`~bilebi00/bin/R --no-save --args $counts < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^TEST.RESULTS" -A 5 | tail -n 1` 
		echo $of RANDOM $p
	done
	rm $of*
done >> ryba.test_results

exit;
#ChimerDB Diff
less $chimerdiff | awk -F "\t" 'NR>2 && $3 == "A" { print $4; }' | sort | uniq > chimerdiff.head
less $chimerdiff | awk -F "\t" 'NR>2 && $3 == "A" { print $9; }' | sort | uniq > chimerdiff.tail
#ChimerDB Same
less $chimersame | awk -F "\t" 'NR>2 && $3 == "A" { print $4; }' | sort | uniq > chimersame.head
less $chimersame | awk -F "\t" 'NR>2 && $3 == "A" { print $9; }' | sort | uniq > chimersame.tail

#Durkin nd Glover 2007 fragile sites
less ~bilebi00/_EWSR1/data/FragileSites/Table2_Durkin_and_Glover_2007_Annu_Rev_Genet  | cut -f 3 | perl -e '$_=<>;while(<>){@t=split/,/; print join("\n", @t); }' | grep -v "Not identified" > durkin_and_glover

for f in ${orybas[@]} durkin_and_glover clipped chimersame.head chimersame.tail chimerdiff.head chimerdiff.tail; do
	echo $f
	grep -w -f $f expressed | sort | uniq > $f.expressed
done

#ChimerDB All 
cat chimer*expressed | sort | uniq > chimerall.expressed

for f in ${orybas[@]} durkin_and_glover chimerall chimersame.head chimersame.tail chimerdiff.head chimerdiff.tail; do
	f=$f.expressed
	~bilebi00/_PAPD5/scripts/innerJoin.pl clipped.expressed $f 1 1 1 > clipped-$f
	echo clipped-$f clipped.expressed $f expressed;
	counts=`wc -l clipped-$f clipped.expressed $f expressed | awk 'NR<5{print $1}'`;
	~bilebi00/bin/R --no-save --args $counts < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep TEST.RESULT -A 3 | grep "^\["
done > test_results



