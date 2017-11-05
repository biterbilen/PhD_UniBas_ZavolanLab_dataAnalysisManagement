#!/bin/bash

nucs="A C G T"

indir=~grubera/snoRNA-paper/rRNA
cs=(LSU SSU 58S)
fas=(bwa_files/NR_003287_LSU_28S.fa bwa_files/NR_003286_SSU_18S.fa bwa_files/NR_003285_5-8S.fa)
#TODO set here
SGE_TASK_ID=2; # 1 2 3
id=mRNASeq_4SU_365_XL_rep_A
id=mRNASeq_4SU_365_XL_rep_B
id=NOP58_repA_Clip22
id=NOP58_repB_Clip22

s="+"
fa=${fas[$((SGE_TASK_ID-1))]}
c=${cs[$((SGE_TASK_ID-1))]}
chr=`perl -e '$a=shift; $chr= ($a =~ /(NR_\d+)/); print $1;' $fa`;
r=$id.$c.$chr.$s
genomefa=$r.fa

#XXX hard coded for NR accessions
less $indir/$fa | perl -e ' $_=<>; ($id) = ($_ =~ /(NR\w+)/); $_=<>;  print ">$id\n$_"; ' > $genomefa
#1 nuc version
less $genomefa | perl -e ' $src = shift; $fea=shift; $_=<>; ($id) = ($_ =~ />(\w+)/); $_=<>; chomp; @nucs= split//,$_; for (my $i=0; $i<@nucs; $i++) { print join("\t", $id, $src, $fea, $i+1, $i+1, ".", "+", ".", "gene_id \"$id\"; transcript_id \"$id\";", $nucs[$i], "\n");} ' GENBANK rRNA | gzip -c > $r.nucbed.gz
#get coords and sequences
less $genomefa | awk 'BEGIN{OFS="\t";} { if (NR>1){ print chr, src, fea, 1, length($0),".", "+", ".", "gene_id \"" chr "\"; transcript_id \"" chr "\";" ; }}' chr=$chr src=GENBANK fea=rRNA | gzip -c > $r.gtf.gz

#get copies coverage
if [ ! -e $id.$c.copies_count.gtf.gz ]; then
	cat $indir/$id.$c | \
		awk 'BEGIN{OFS="\t"; } { print chr,$2-1,$3,$1,$4,strand; }' chr=$chr strand=$s | \
		~bilebi00/www/MARA/scripts/extend_bed_fromScoreField.pl 1 | \
		bedtools coverage -d -s -a stdin -b $r.gtf.gz | \
		awk -F "\t" 'BEGIN{OFS="\t";} { $4=($10-$4+1); $5=$4; $10=1; print $0; }' | \
		gzip -c > $id.$c.copies_count.gtf.gz
fi
for nuc in $nucs; do
	if [ ! -e $r.copies${nuc}_count.gtf.gz ]; then
			zless $id.$c.copies_count.gtf.gz | ~bilebi00/_EWSR1/scripts/get_nuc_coverage.pl  $r.nucbed.gz $nuc 0 copies_count 1 | \
			gzip -c > $r.copies${nuc}_count.gtf.gz
	fi			
done

echo coverage for copies finished 


#get $nucs mutations coverage
for nuc in $nucs; do
	if [ ! -e $r.mut${nuc}_count.gtf.gz ]; then
		cat $indir/$id.$c | \
			perl -e '$chr=shift; $s=shift; $nuc=shift; $strand=shift; while(<>){ @t=split; @ps=split/;/,$t[5]; @ms=split/;/,$t[6]; for (0..$#ms) { if ($ms[$_]=~/M$nuc/) { print join("\t", $chr, $ps[$_]-1, $ps[$_], $t[0], $t[3], $strand ),"\n";  } }}' $chr $s $nuc $s | \
			~bilebi00/www/MARA/scripts/extend_bed_fromScoreField.pl 1 | \
			bedtools coverage -d -s -a stdin -b $r.gtf.gz | \
			awk -F "\t" 'BEGIN{OFS="\t";} { $4=($10-$4+1); $5=$4; $10=1; print $0; }' | \
			~bilebi00/_EWSR1/scripts/get_nuc_coverage.pl  $r.nucbed.gz $nuc 0 mut_count 1 | \
			gzip -c > $r.mut${nuc}_count.gtf.gz
	fi
done

for nuc in $nucs; do
	plot=0
	if [ $nuc == "T" ] && [ "`echo $id | grep -i RNAseq`" == "" ] ; then
		plot=1
	fi
	for mu in 0.02 0.05 0.08 0.1; do
		~bilebi00/bin/R --no-save --args $r.copies${nuc}_count.gtf.gz $r.mut${nuc}_count.gtf.gz $mu $plot < ~bilebi00/_EWSR1/scripts/plot_region_fromGtfPlus.R > /dev/null;
		#~bilebi00/bin/R --no-save --args $id.$c.copies${nuc}_count.gtf.gz $id.$c.mut${nuc}_count.gtf.gz $mu 1 < ~bilebi00/_EWSR1/scripts/plot_region_fromGtfPlus.R > /dev/null;
	done
done

echo coverage for mut finished 

less $id.$c.*freqstat | sort -r | uniq > $id.$c.stat
rm $id.$c.*freqstat

exit;
#
#
#cat *stat | sort -r | uniq > all.stat
#R --no-save --args all.stat < ~bilebi00/_EWSR1/scripts/CI.R

~bilebi00/_PAPD5/scripts/outerJoin.pl NOP58_repA_Clip22.58S*mutT_mutrate0.02.txt NOP58_repB_Clip22.58S*mutT_mutrate0.02.txt \
	  1,2,3,4,5,6,7,8,9 1,2,3,4,5,6,7,8,9 '1-9,12,10-11,13-14,24-25,27-28' NA > _NOP58_repA_repB.58S.mutT_mutrate0.02.txt
~bilebi00/bin/R --no-save --args  _NOP58_repA_repB.58S.mutT_mutrate0.02.txt < ~bilebi00/_EWSR1/scripts/get_truested_pbeta.R > /dev/null

~bilebi00/_PAPD5/scripts/outerJoin.pl NOP58_repA_Clip22.LSU*mutT_mutrate0.02.txt NOP58_repB_Clip22.LSU*mutT_mutrate0.02.txt \
	  1,2,3,4,5,6,7,8,9 1,2,3,4,5,6,7,8,9 '1-9,12,10-11,13-14,24-25,27-28' NA > _NOP58_repA_repB.LSU.mutT_mutrate0.02.txt
~bilebi00/bin/R --no-save --args  _NOP58_repA_repB.LSU.mutT_mutrate0.02.txt < ~bilebi00/_EWSR1/scripts/get_truested_pbeta.R > /dev/null

~bilebi00/_PAPD5/scripts/outerJoin.pl NOP58_repA_Clip22.SSU*mutT_mutrate0.02.txt NOP58_repB_Clip22.SSU*mutT_mutrate0.02.txt \
	  1,2,3,4,5,6,7,8,9 1,2,3,4,5,6,7,8,9 '1-9,12,10-11,13-14,24-25,27-28' NA > _NOP58_repA_repB.SSU.mutT_mutrate0.02.txt
~bilebi00/bin/R --no-save --args  _NOP58_repA_repB.SSU.mutT_mutrate0.02.txt < ~bilebi00/_EWSR1/scripts/get_truested_pbeta.R > /dev/null
#---------------
#0.0256788
#is the average mutT rate in mRNAseq
#cat all.stat | grep mutT | grep mRNASeq | less| awk 'BEGIN{c=0; k=0;}{c=c+($2*$4);k=k+$4; print c/k;}' | tail -n 1

