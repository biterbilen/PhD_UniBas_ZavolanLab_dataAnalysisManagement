library(lattice)
library(reshape)

cargs <- commandArgs()
#cargs <- c(1,1,1,"../Analysis/Reproducibility/TIA1.kmer.stat")
header <- c("count", "protein", "fea", "k", "lib", "scr", "kmer")
a <- read.table(cargs[4],header=F)
names(a) <- header;
protein <- as.character(a$protein[1])

#cargs <- c(1,1,1,"../Analysis/Reproducibility/kmer.stats")
#a <- read.table(cargs[4],header=T)

scrs <- sort(unique(a$scr), decreasing=T);
feas <- levels(sort(unique(a$fea), decreasing=T));
ks <-   sort(unique(a$k), decreasing=T);
libs <-   sort(unique(a$libs), decreasing=T);

newheader <- c("protein", "scr", "fea", "k", "corr")
m <- matrix(0, nrow=length(scrs)*length(feas)*length(ks), ncol=length(newheader));

i <- 1;
for (scr in scrs) {
	for (fea in feas) {
		for (k in ks) {
			suba <- a[a$protein==protein & a$scr>=scr & a$fea==fea & a$k==k,c("kmer","lib","count")];
			ma <- reshape(suba, direction="wide", timevar="lib", idvar=c("kmer"))
			corr <- NA
			if (dim(ma)[2] == 3) {
				yy <- aggregate(ma[,-1], by=list(ma$kmer), sum);
				corr <- cor(yy[,2], yy[,3], use="pairwise.complete.obs")
			}
			m[i,] <- c(protein, scr, fea, k, corr)
			i <- i+1;
		}
	}
}
head(m)
aplot <- as.data.frame(m, col.names=newheader, stringsAsFactors=F)
names(aplot) <- newheader;
head(aplot)

#write.table(cbind(T2C,m),quote=F,file=paste(cargs[4],'.cumcor', sep=''),sep="\t",row.names=F)
#write.table(aplot,quote=F,file=paste(cargs[4],'.cum', sep=''),sep="\t",row.names=F);
ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)      ## set as default

#write.table(aplot,quote=F,file=paste(cargs[4],'.cum', sep=''),sep="\t",row.names=F);

#kmer reproducibility
plot.new(); pdf(paste(cargs[4], '.kmercor','.pdf',sep=''), width=8, height=8)
barchart( as.numeric(corr) ~ fea | factor(as.numeric(k)), 
	group = as.numeric(scr), data=aplot, origin = 0,
	auto.key = list(title = "Minimum Feature Frequency", space="bottom",columns=3, cex=1.5,font=2),
	par.strip.text = list(cex=1.5,font=2),layout=c(1,3),
	par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1.5,
	main=paste("Correlation of kmers around non-reproduced position -", protein, sep=''),
	ylab="Proportion (%)", scales=list(x=list(rot=90)))
dev.off();

