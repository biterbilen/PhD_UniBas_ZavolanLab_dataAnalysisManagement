getDataFromGtf <- function (filename) {
	aa <- read.table(filename, header=F);

	im <- 1:8
	io <- 9:dim(aa)[2]

	ns <- aa[1,io[io %% 3==0]]
	iv <- io[io %% 3==1]

	a <- aa[,c(im,iv)];
	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)));
	return(a)
}

#library(plotrix)
#library(hexbin)
#library(reshape)

cargs <- commandArgs();
#cargs <- c(1,1,1,"~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage.meanCoverage/271.introns.gtf", "~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage.meanCoverage.meanCoverage/271.exons.gtf","mRNAseq_EWSR1_siCTRL_1");
x <- getDataFromGtf(cargs[4]);
y <- getDataFromGtf(cargs[5]);

#get exon information
yy <- unique(y[,c("gene_id", names(y)[grep("exon",names(y))])]);
a <- merge(x, yy, by=c("gene_id"))

a$sample <- cargs[6];

#densities are in log scale
a$gene_exon_density   <- log2((a$all_exons_RNAseq_coverage+1) / a$all_exons_length)
a$gene_intron_density <- log2((a$all_introns_RNAseq_coverage+1) / a$all_introns_length)
a$intron_length  <- a$end - a$start + 1
a$intron_density <- log2((a$RNAseq_coverage + 1) / a$intron_length)

head(a,2)
a$intron_exp_over_exons   <- a$intron_density - a$gene_exon_density;
a$intron_exp_over_introns <- a$intron_density - a$gene_intron_density;
a$intron_score <- a$intron_exp_over_exons + a$intron_exp_over_introns;
gexp <- quantile(a$gene_exon_density);

dim(a)
#restrictions on: 
#expression, gene count in locus, gene is xlinked, region length
b <- a[a$gene_exon_density>gexp[2] & a$score==1 & a$xlink_sites>-1 & a$intron_length>1000,]
dim(b)
a <- b


#TODO uncomment these lines later
write.table(a[order(a$intron_score,decreasing=T),c("sample", "gene_id","seqname", "start", "end", "strand", "region_xlink_sites", "intron_length", "gene_exon_density", "intron_exp_over_exons", "intron_exp_over_introns", "intron_score")],
	file=paste(cargs[4],'.all.txt',sep=''),row.names=F,quote=F,sep="\t");

library(latticeExtra)
library(hexbin) 

p1 <- xyplot(log2(intron_length) ~ intron_score, data= a, type=c("p","g"), panel=panel.hexbinplot); #short introns excluded
p2 <- xyplot(intron_score ~ gene_exon_density, data= a, type=c("p","g"), panel=panel.hexbinplot); #gene exp is not correlated to intron score
p3 <- xyplot(log2(intron_length) ~ log2(region_xlink_sites+1), data= a, type=c("p","g"), panel=panel.hexbinplot); # 
p4 <- xyplot(intron_exp_over_exons ~ intron_exp_over_introns, data=a, type=c("p","g"), panel=panel.hexbinplot);

plot.new(); pdf(paste(cargs[4], "_params.pdf",sep=''), h=8, w=8);
plot(p1, split = c(1, 1, 2, 2)) 
plot(p2, split = c(1, 2, 2, 2), newpage = FALSE) 
plot(p3, split = c(2, 1, 2, 2), newpage = FALSE) 
plot(p4, split = c(2, 2, 2, 2), newpage = FALSE) 
dev.off();


#update(c(p1,p2,p3,p4), x.same=F, y.same=F, merge.legends=F );


quit();


