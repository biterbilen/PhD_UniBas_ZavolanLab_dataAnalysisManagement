library(latticeExtra)
#library(plotrix)
#library(hexbin)
#library(reshape)

cargs <- commandArgs();
#cargs <- c(1,1,1,'regionFinal.plot');
#cargs <- c(1,1,1,"~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/sRNAseq.joined");
aa <- read.table(cargs[4], header=T);
head(aa,2)

#a <- aa[aa$intronLength>0,]
a <- aa;
a$intronDensity <- (a$regionCoverage+1)/a$regionLength;
a$exonDensity <- (a$pexonCoverage+1 + a$fexonCoverage+1) / (a$pexonLength + a$fexonLength);
a$intronDensity.1 <- (a$regionCoverage.1+1)/a$regionLength.1;
a$exonDensity.1 <- (a$pexonCoverage.1+1 + a$fexonCoverage.1+1) / (a$pexonLength.1 + a$fexonLength.1);

#a$geneCLIPstat <- cut(a$xlinkSites, breaks=c(-1,0,max(a$xlinkSites)), labels=c("NOTCLIPPEDGENE","CLIPPEDGENE"))
a$regionCLIPstat <- cut(a$regionXlinkSites, breaks=c(-1,0,max(a$regionXlinkSites)), labels=c("NOTCLIPPEDREGION","CLIPPEDREGION"))
a$batch <- "BATCH0"; 
a$batch[grep("_1$", a$name.1, perl=T)] <- "BATCH1"; 
a$batch[grep("_2$", a$name.1, perl=T)] <- "BATCH2"; 
print(unique(a$batch));
dim(a)

#b <- a[a$exonDensity>1,]
b <- a[a$xlinkSites>0,] #xlinked genes
a <- b;
dim(a)

samples <- levels(unique(a$name.1))

alt <- "less";
#alt <- "gr";
for (s in samples) {
	xi <- a$name.1==s & a$regionCLIPstat=="CLIPPEDREGION";
	x <- log2((a[xi,"intronDensity"] / a[xi,"exonDensity"]) / (a[xi,"intronDensity.1"] / a[xi,"exonDensity.1"]));
	yi <- a$name.1==s & a$regionCLIPstat=="NOTCLIPPEDREGION";
	y <- log2((a[yi,"intronDensity"] / a[yi,"exonDensity"]) / (a[yi,"intronDensity.1"] / a[yi,"exonDensity.1"]));
	ks <- ks.test(x,y,alternative=alt)
	write.table(data.frame(list(method="Two sample KS test for CLIPPEDREGION and NOTCLIPPEDREGION RNAseq expression", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F);
}

#dim(a[a$clipstat=="CLIPPED",])
plot.new(); pdf(paste(cargs[4], '_ecdf.pdf',sep=''), h=8, w=8);
plt <- ecdfplot(~ log2((intronDensity / exonDensity) / (intronDensity.1 / exonDensity.1)) | factor(batch), groups=paste(name.1, regionCLIPstat), data=a,
	subset=locusGeneCount==1, horizontal=F,
	main=paste("EWSR1 sRNAseq tag density over mRNAseq&&GROseq tag density", "for regions of xlinked multiexon genes v CLIP status",sep="\n"),
	scales = list(x = list(rot = 90)),
 	auto.key=list(space = "top", title = "samples") )
update(plt,
	panel = function(...) {
		panel.grid(h = -1, v = -1)
		panel.ecdfplot(...)
	})

dev.off()
quit();


