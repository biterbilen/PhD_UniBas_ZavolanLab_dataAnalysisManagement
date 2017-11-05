library(latticeExtra)
#library(plotrix)
#library(hexbin)
#library(reshape)

#TODO some parts added for tagged with high
quit();

cargs <- commandArgs();
#cargs <- c(1,1,1,'regionFinal.plot');
#cargs <- c(1,1,1,'a.plot');
aa <- read.table(cargs[4], header=T);
head(aa,2)

#a <- aa[aa$intronLength>0,]
a <- aa;
a$intronDensity <- (a$regionCoverage+1)/a$regionLength;
a$exonDensity <- (a$pexonCoverage+1 + a$fexonCoverage+1) / (a$pexonLength + a$fexonLength);

#high
dim(a)
qnt <- quantile(a$intronDensity);
b <- a[a$intronDensity > qnt[4],];
a <- b;
dim(a)

#a$geneCLIPstat <- cut(a$xlinkSites, breaks=c(-1,0,max(a$xlinkSites)), labels=c("NOTCLIPPEDGENE","CLIPPEDGENE"))
a$regionCLIPstat <- cut(a$regionXlinkSites, breaks=c(-1,0,max(a$regionXlinkSites)), labels=c("NOTCLIPPEDREGION","CLIPPEDREGION"))
a$batch <- "BATCH0"; 
a$batch[grep("_1$", a$name, perl=T)] <- "BATCH1"; 
a$batch[grep("_2$", a$name, perl=T)] <- "BATCH2"; 
print(mean(a$intronDensity));

#b <- a[a$exonDensity>1,]
dim(a)
b <- a[a$xlinkSites>0,] #xlinked genes
a <- b;
dim(a)

samples <- levels(unique(a$name))

alt <- "less";
#alt <- "gr";
for (s in samples) {
	xi <- a$name==s & a$regionCLIPstat=="CLIPPEDREGION";
	x <- log2(a[xi,"intronDensity"] / a[xi,"exonDensity"]);
	yi <- a$name==s & a$regionCLIPstat=="NOTCLIPPEDREGION";
	y <- log2(a[yi,"intronDensity"] / a[yi,"exonDensity"]);
	ks <- ks.test(x,y,alternative=alt)
	write.table(data.frame(list(method="Two sample KS test for CLIPPEDREGION and NOTCLIPPEDREGION RNAseq expression", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F, quote=F);
}

#dim(a[a$clipstat=="CLIPPED",])
#high
#plot.new(); pdf(paste(cargs[4], '_ecdf.pdf',sep=''), h=8, w=8);
plot.new(); pdf(paste(cargs[4], "_ecdf_highly_exp.pdf",sep=''), h=8, w=8);
plt <- ecdfplot(~ log2(intronDensity / exonDensity) | factor(batch), groups=paste(name, regionCLIPstat), data=a,
	subset=locusGeneCount==1 , horizontal=F,
	#high
	main="EWSR1 RNAseq tag density for highly expressed regions \nof xlinked multiexon genes v CLIP status",
	#main="EWSR1 RNAseq tag density for regions of xlinked multiexon genes v CLIP status",
	scales = list(x = list(rot = 90)),
 	auto.key=list(space = "top", title = "samples") )
update(plt,
	panel = function(...) {
		panel.grid(h = -1, v = -1)
		panel.ecdfplot(...)
	})

dev.off()
quit();

