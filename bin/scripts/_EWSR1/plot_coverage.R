library(latticeExtra)
#library(plotrix)
#library(hexbin)
#library(reshape)

cargs <- c(1,1,1,'final.plot');
aa <- read.table(cargs[4], header=T);
a <- aa[aa$intronLength>0,]
a$exonDensity <- (a$exonRNAseqCoverage+1)/a$exonLength;
a$intronDensity <- (a$intronRNAseqCoverage+1)/a$intronLength;
a$clipstat <- cut(a$xlinkSites, breaks=c(-1,0,max(a$xlinkSites)), labels=c("NOTCLIPPED","CLIPPED"))
print(mean(a$exonDensity));
print(mean(a$intronDensity));
dim(a)
b <- a[a$exonDensity>1,]
a <- b;
dim(a)
dim(a[a$clipstat=="CLIPPED",])
plot.new(); pdf(paste(cargs[4], '.pdf',sep=''), h=8, w=8);
#xyplot(exonDensity ~ intronDensity | name, data=a,
plt <- bwplot(log2(intronDensity / exonDensity) ~ paste(name,clipstat,sep="\n"), data=a,
#barchart(mean(log2(intronDensity / exonDensity)) ~ paste(name,clipstat,sep="\n"), data=a,
#	panel=function(...) { panel.grid(h = -1, v = -1); panel.bwplot(...) },
	subset=locusGeneCount==1, horizontal=F,
	main="EWSR1 RNAseq tag density v CLIP status for expressed multiexon genes",
	scales = list(x = list(rot = 90)),
 	auto.key=list(space = "top", title = "samples") )
update(plt,
	panel = function(...) {
		panel.grid(h = -1, v = -1)
		panel.bwplot(...)
	})

dev.off()
quit();

