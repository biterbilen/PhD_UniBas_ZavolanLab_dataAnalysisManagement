#!/bin/bash

rd=~bilebi00/_EWSR1/data/ReplicationDomain/;
rybas=(${rd}GE_BG01_ESC_All.txt ${rd}GE_BG01_NPC_All.txt ${rd}GE_BG02_ESC_All.txt ${rd}GE_H7_ESC_All.txt)

tag=EWSR1
subregion=intron
chimerdiff=~bilebi00/_EWSR1/data/ChimerDB/mRNA_Information_Diff.txt
chimersame=~bilebi00/_EWSR1/data/ChimerDB/mRNA_Information_Same.txt
expf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/annotation_xlinkEnrichment_topN10000_cellHEK293/exp_genes
clipf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/annotation_xlinkEnrichment_topN10000_cellHEK293/$subregion/common_subregion_genes

soutdir=Analysis/ChimerDB
mkdir -p $soutdir; pushd $soutdir;

#expressed
ln -sf $expf expressed
ln -sf $clipf clipped

#Replication domain
#Ryba T, Hiratani I, Lu J, Itoh M et al. Evolutionarily conserved replication timing profiles predict long-range chromatin interactions and distinguish closely related cell types. Genome Res 2010 Jun;20(6):761-70.
orybas=();
for ryba in ${rybas[@]}; do
	N=`less $ryba | awk -F "\t" 'NF==7 && $7==1 {print $0;}' | cut -f 2 | sort | uniq | wc -l`;
	f=`echo "$N / 20 " | bc`; #5% of the top
	of=ryba_`basename $ryba _All.txt`;
	less $ryba | awk -F "\t" 'NF==7 && $7==1 {print $0;}' | sort -k 6,6gr | cut -f 2 | head -n $f | sort | uniq > $of.early;
	less $ryba | awk -F "\t" 'NF==7 && $7==1 {print $0;}' | sort -k 6,6g | cut -f 2 | head -n $f | sort | uniq > $of.late;
	#less $ryba | perl -e 'while(<>){ chomp; @t=split/\t/; next if (@t !=7 or $t[6] != 1); $r=rand(); print $t[1],"\n" if ($r>0.5);  }' | head -n $f | sort | uniq > $of.random;
	orybas=(${orybas[@]} $of.early $of.late);
done

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
	~bilebi00/_PAPD5/scripts/innerJoin.pl clipped.expressed $f 1 1 1 > $f-clipped.expressed
	echo -en $f-clipped.expressed clipped.expressed $f expressed '\t';
	counts=`wc -l $f-clipped.expressed clipped.expressed $f expressed | awk 'NR<5{print $1}'`;
	echo `~bilebi00/bin/R --no-save --args $counts < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" | tail -n 5`
done > test_results



