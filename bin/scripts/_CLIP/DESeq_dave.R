library(DESeq)
library("RColorBrewer")
library("gplots")

#Uses DEseq for calling differentially expressed genes
#Expression values are normalized using geometric mean of the genes as a reference sample
#New variance stabilizing transformation feature is used
#Variance in designs wo/replicates is fit by ignoring the sample labels and outliers are not

cargs <- c("input_4_deseq","selected","expressed.id.gid", "siCTRL", "0.2")
cargs <- c("input_4_deseq","selected","expressed.id.gid", "siUntreated", "0.2");
cargs <- c("input_4_deseq.anneke.1","selected.anneke.1","expressed.id.gid.anneke","Ago", "0.2");
cargs <- c("input_4_bayseq.1","selected.EWSR1.1","expressed.id.gid.EWSR1","siCTRL", "0.2");
cargs <- c("input_4_deseq","selected","expressed.id.gid","siCTRL", "0.1");
cargs <- c("input_4_deseq","selected","null","cont","0.1")
cargs <- commandArgs(trailingOnly=T)

gid <- NULL
cond1 <- "siUntreated"
padj <- 0.2
fitType <- "parametric"
a   <- read.delim(cargs[1], header=T, sep="\t")
design <- read.delim(cargs[2], header=T, sep="\t")
if (length(cargs) > 2 && !(cargs[3]=="null" || cargs[3]=="NULL")) 
	gid <- read.table(cargs[3],header=F,col.names=c("ids","gid"))
if (length(cargs) > 3) 
	cond1 <- cargs[4]
if (length(cargs) > 4) 
	padj <- as.numeric(cargs[5])
if (length(cargs) > 5) 
	fitType <- cargs[6]

print(head(gid,1))
print(dim(a))
print(fitType)

#filter out not expressed genes
a <- a[rowSums(a[grep("exp",names(a))])>0,];

#ids vector
ids  <- a[,1]

#length matrix 
#lens <- as.matrix(a[,grep("len",names(a))])

#read count matrix
all  <- a[,grep("exp",names(a))]

#dimensions of data
ngenes   <- dim(all)[1]
nsamples <- dim(all)[2]

#total expression
sums     <- colSums(all)

#remove tag from names
names(all) <- gsub(".exp","",names(all))

rownames(all) <- ids
print(dim(all))
print(head(all,2))
print(design)

cds <- newCountDataSet(all,design$condition)
cds <- estimateSizeFactors(cds)
#This is the published method which uses local dispersion and estimates the variance of each gene ignoring the labels
  
print(length(unique(design$condition)))
print(length(design$condition));

if (length(unique(design$condition)) == length(design$condition)) {
	if (fitType == "local")
		#local dispersion fit as in the paper
		cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local") 
	else
		#parametric dispersion fit: new method
		cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only")
} else {
	if (fitType == "local")
	#local dispersion fit as in the paper
		cds <- estimateDispersions(cds,fitType="local")
	else
	#parametric dispersion fit: new method
		cds <- estimateDispersions(cds)
}

print(str(fitInfo(cds)))

plot.new(); pdf(paste(cargs[1],"_",cond1,".pdf", sep=""),h=10,w=10) 
plotDispEsts(cds, main="Dispersion estimates")
#sizeFactors(cds)
for (cond2 in unique(design$condition)) {
	if (cond1 == cond2) {
		all_norm <- counts(cds, normalized=TRUE) 
		hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#		heatmap( as.matrix( dist(t(all)) ), symm=TRUE, main="raw" )
		heatmap.2( as.matrix( dist(t(all_norm)) ), main="normalized", col = hmcol, trace="none", margin=c(10,6))
	} else {
		#vsd <- getVarianceStabilizedData( cds )

		res <- nbinomTest (cds, cond2, cond1)
		plotMA(res, main=paste(cond2,"vs",cond1))
		hist(res$pval, breaks=100, col="skyblue", border="slateblue", main=paste(cond2,"vs",cond1))
		hist(res$padj, breaks=100, col="skyblue", border="slateblue", main=paste(cond2,"vs",cond1))

		if (is.null(gid)) {
			deseq_DE <- res;
		} else {
			deseq_DE <- merge(res, gid, by.x="id", by.y="ids")
		}
		ordered_deseq_DE <- deseq_DE[order(deseq_DE$padj),]
		write.table(ordered_deseq_DE, file=paste(cargs[1],"_", cond1,"_", cond2,"_DE.txt",sep=""),sep="\t",quote=F,row.names=F)
		if (!is.null(gid)) {
			write.table(ordered_deseq_DE[ordered_deseq_DE$padj < padj, "gid"], file=paste(cargs[1],"_", cond1,"_", cond2,".genenames",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
		}
	}
}
dev.off()

