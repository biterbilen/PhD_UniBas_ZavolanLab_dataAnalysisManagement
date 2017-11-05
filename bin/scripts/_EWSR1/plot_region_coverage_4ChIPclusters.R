cargs <- commandArgs();

#cargs <- c(1,1,1,"/scratch/ENCODE23/chip.CTCF.gtfplus",'CTCF',1);

a <- read.table(cargs[4],header=F);
head(a,1)
names(a) <- sub("V6", "locusGeneCount", names(a))
names(a) <- sub("V16", "xlinkSites", names(a))
names(a) <- sub("V19", "regionXlinkSites", names(a))
names(a) <- sub("V25", "fraction", names(a))
head(a,1)

dim(a)
b <- a[a$xlinkSites>0,] #xlinked genes
a <- b;
dim(a)

a$regionCLIPstat <- cut(a$regionXlinkSites, breaks=c(-1,0,max(a$regionXlinkSites)), labels=c("NOTCLIPPEDREGION","CLIPPEDREGION"));

head(a,1)
alt <- "less"
ps <- 0.0001;
x <- log2(a[a$locusGeneCount==1&a$regionCLIPstat=="CLIPPEDREGION",c("fraction")]+ps);
y <- log2(a[a$locusGeneCount==1&a$regionCLIPstat=="NOTCLIPPEDREGION",c("fraction")]+ps);
ks <- ks.test(x,y, alternative=alt)
ks
ks$p.value
ks$statistic
write.table(data.frame(list(pvalue=ks$p.value, statistic=ks$statistic)),quote=F, sep="\t",row.names=F);
write.table(data.frame(list(method="Two sample KS test for CLIPPEDREGION and NOTCLIPPEDREGION ChIPseq coverage fraction", alternative=alt, chipSeqSampleName=cargs[5], pvalue=ks$p.value, statistic=ks$statistic)),quote=F, sep="\t",row.names=F);


if (length(cargs)<6 | cargs[6] == 0) {
	quit();
}

library(latticeExtra)

plot.new(); pdf(paste(cargs[4], "_ecdf_",cargs[5],"_coverage_fraction.pdf",sep=''), h=8, w=8);
plt <- ecdfplot(~ log2(fraction + ps), groups=paste(regionCLIPstat), data=a,
	subset=locusGeneCount==1 , horizontal=F,
	main=paste(cargs[5], " coverage fraction for regions of\nEWSR1 xlinked multiexon genes", sep=''),
	scales = list(x = list(rot = 90)),
	auto.key=list(space = "top", title = "CLIP status") )
update(plt,
	panel = function(...) {
		panel.grid(h = -1, v = -1)
		panel.ecdfplot(...)
	})
dev.off()
