for i in Cuffcompare_*/Cuffcompare*.trusted.gtf*; do 
	echo $i; 
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/get_stats_from_gtf.pl $i tss_id transcript_id gene_id gene_name class_code
	echo; 
done | less > CUFFCOMPARE.stats
