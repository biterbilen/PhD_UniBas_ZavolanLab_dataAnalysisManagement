#!/bin/bash

s=68;
c=ALL;
f=T2C;

#index starts from 0

~/CDS_CLIP_ElMMo/scripts/generateGenomeAlignmets.pl $s 1 + $f $c > $s.$c
~bilebi00/CDS_CLIP_ElMMo/scripts/generateCumulativeScores.pl $s.$c | awk 'BEGIN{OFS="\t";}{print "chr14",$0;}' | sort -k 4,4gr > $s.$c.old

~/GEORGES_POLYA/scripts/copy_allFeaturesToWiggle.pl $s 1 + $s$c $c 
#~/GEORGES_POLYA/scripts/allFeaturesToWiggle.pl $s 1 + $s$c $c 
#less $s$c/$f | awk '{if ($4 != 0) { print $0;} }' | sort -k 4,4gr > $s$c.$f
less $s$c/$f | sort -k 4,4gr > $s$c.$f

wc -l $s.$c.old $s$c.$f
head $s.$c.old $s$c.$f
tail $s.$c.old $s$c.$f


exit;
#('SELECT * FROM (SELECT alignment, sequence, chr_start, chr_end, 
#			copies/genome_count_total*(1000000/((SELECT sum(copies) FROM t_mapped_sequence_'.$ARGV[0].' ) + 
#							    (SELECT sum(copies) FROM t_unmapped_sequence_'.$ARGV[0].' ) )) as copies 
#			from (SELECT alignment, chr_start, chr_end, t_sequence_id FROM t_sequence_genome_alignment_'.$ARGV[0].' 
#								WHERE t_chromosome_contig_id = ? AND strand = ? '.$annot.') AS AL, 
#			t_mapped_sequence_'.$ARGV[0].' AS M WHERE M.id = t_sequence_id) AS K ORDER BY chr_start');
#
#('SELECT alignment, sequence, chr_start, chr_end, 
#	copies/genome_count_total*(1000000/((SELECT sum(copies) FROM t_mapped_sequence_'.$ARGV[0].' ) + 
#					    (SELECT sum(copies) FROM t_unmapped_sequence_'.$ARGV[0].' ) )) as copies 
#	from t_sequence_genome_alignment_'.$ARGV[0].' AS AL, t_mapped_sequence_'.$ARGV[0].' AS M 
#					WHERE t_chromosome_contig_id = ? AND strand = ? AND M.id = t_sequence_id '.$annot.' ORDER BY chr_start')
#
