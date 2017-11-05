for t in gene_exp isoform_exp promoters splicing tss_group_exp; do
	i=Cuffdiff_trusted_RNAseq/$t.diff;
	echo $i; 
	awk ' NR>1 && $NF == "yes" { print $0; }' $i | cut -f 5-7,14 | sort | uniq -c | sort -k 1,1gr;
	echo; 
done | less > CUFFDIFF.stats
exit;

for t in gene_exp isoform_exp promoters splicing tss_group_exp; do
	i=w_old_trusted_Cuffdiff_trusted_RNAseq/$t.diff;
	echo $i; 
	awk ' NR>1 && $NF == "yes" { print $0; }' $i | cut -f 5-7,14 | sort | uniq -c | sort -k 1,1gr;
	echo; 
done | less > CUFFDIFF_w_old_trusted.stats
exit;
