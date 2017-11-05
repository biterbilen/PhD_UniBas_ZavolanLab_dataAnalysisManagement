#Number of DIS3L2 and AGO2 joint targets for 72 nc (as well as all tRNAs) is significant (check presentation -- p-value = 3.8e-40, odd-ratio=12)

#RUN MIRZA; mirza
#pushd ./launch_MIRZA ./file_samples/miRNA_expression ./file_samples/miRNA_sequence ./file_samples/mRNA_sequence noupdate

odir=Analysis/xlinkEnrichment/prep4Mirza
mkdir -p $odir; pushd $odir

#statistical significance of the joint 72nc-long AGO2 and DIS3L2 tRNAs targets
for i in ../xlinkPositionEnrichment_tRNA/*/*nucbed*; do echo $i; zless $i | wc -l; done
for i in ../xlinkPositionEnrichment_tRNA/*/*nucbed*; do zless $i ; done | sort | uniq -D | uniq | wc -l
R --no-save --args 50 68 76 209 < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" > fisher.stats 

#prep guides
sslen=21
for i in ../xlinkPositionEnrichment_tRNA/*/*nucbed*; do zless $i ; done | \
	sort | uniq -D | uniq | \
	~bilebi00/_DIS3L2/scripts/split_seq_from_nucbed.pl $sslen 1 gtf 9 \
 	> tiRNA_sequence
# 	> tiRNA_sequence.w.originalIds
##patch for mirza getting coredump for long names
#less tiRNA_sequence.w.originalIds | grep -v "^>" | awk 'BEGIN{c=0;}{c=c+1; print ">"c"\n"$1; }' > tiRNA_sequence
#less tiRNA_sequence.w.originalIds | grep "^>" | awk 'BEGIN{c=0;}{c=c+1; print $1"\t"c; }' > tiRNA_sequence.id.givenId

less tiRNA_sequence | awk '$1 ~ /^>/ { print $1"\t"expression_prior;} ' expression_prior=1 | cut --bytes 2- \
	> tiRNA_expression

#intersection of mRNA (exon) targets of DIS3L2 and AGO2 is significant
#TODO N is approximate
Napproximate=10000
cut=1 #added for practical reasons; to patch segmentation fault giving mirza
less ../hyperRegionStats/AGO2/pernuc_merged_XL.genes_exon | awk '$1>cut{print}' cut=$cut > AGO2_exon
less ../hyperRegionStats/DIS3L2/pernuc_merged_XL.genes_exon | awk '$1>cut{print}' cut=$cut > DIS3L2_exon
cat *exon | awk '{ print $3;}' | sort | uniq -D | uniq | sed 's/";*//g' > common_gene_names
wc -l common_gene_names
wc -l *exon
R --no-save --args 215 534 1731 $Napproximate < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" >> fisher.stats
R --no-save --args 49 782 141 $Napproximate < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" >> fisher.stats

#prep targets
#get representative transcripts to get the sequences
~bilebi00/_xlinkEnrichment/scripts/grep_f_w_like_w_column.pl common_gene_names ~/_DIS3L2/Analysis/clipz_RNAseq/FoldChange_wSN_woMM_w_clip_w_chip/_RNAseq_common_scale.raw.w_geneSymbol.w_clip.w_chip 8 1 | grep "^NM_" | cut -f 1,9 > trx.common_gene_names

ln -sf /import/bc2/home/zavolan/bilebi00/../GROUP/miRNA/RefSeq/DB_11-01-18/hg/hg_3UTR.fa .

sslen2=51
shiftlen=$(($sslen2-$sslen))
~bilebi00/_xlinkEnrichment/scripts/grep_fa_given_id.pl trx.common_gene_names  hg_3UTR.fa | \
	~bilebi00/_DIS3L2/scripts/split_seq_from_nucbed.pl $sslen2 $shiftlen dummy 1 0 2 \
	> mRNA_sequence_$cut


popd $odir
