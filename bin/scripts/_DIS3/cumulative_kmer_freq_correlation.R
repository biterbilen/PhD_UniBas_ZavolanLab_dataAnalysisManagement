library(lattice)
cargs <- commandArgs()
#cargs <- c(1,1,1,"~bilebi00/_ARE/Analysis/Reproducibility/combined.freq7mer", "7mer");
a <- read.table(cargs[4], header=T)
head(a)
T2C <- sort(unique(a$minT2C), decreasing=T);

cc <- dim(a)[2];
n <- cc - 2;

libnames <- names(a[3:cc]);

pairnames <- 1:((n*n - n) / 2); 

k <- 1; 
for (i in 2:length(libnames)) {
	for (j in 1:(i-1)) {
		pairnames[k] <- paste(libnames[i], 'v', libnames[j], sep=' '); 	
		k <- k+1;
	}	
}


m <- matrix(0, nrow=length(T2C), ncol=length(pairnames));
rownames(m) <- T2C;
colnames(m) <- pairnames;

aplot <- NULL;
k <- 1;
t2c <- T2C[2];
for (t2c in T2C) {
	suba <- a[a$minT2C>=t2c,2:cc];
	yy <- aggregate(suba[,-1], by=list(suba$kmer), sum);
#	s <- matrix(rep(colSums(yy[,-1]), times=dim(yy)[1]), nrow=dim(yy)[1], byrow=T);
#	yy[,-1] <- yy[,-1]/s
	yy.na <- replace(yy,yy==0,NA);
#	corsuba <- cor(yy.na[,-1], method="spearman", use="pairwise.complete.obs")
	corsuba <- cor(yy.na[,-1], method="pearson", use="pairwise.complete.obs")
	m[k,] <- as.vector(corsuba[upper.tri(corsuba)]) 
	k <- k + 1;
	yy.na$kmer <- yy.na$Group.1;
	yy.na$Group.1 <- NULL;
	yy.na$T2C <- t2c;
	aplot <- rbind(aplot,yy.na)
}
write.table(cbind(T2C,m),quote=F,file=paste(cargs[4],'.cumcor', sep=''),sep="\t",row.names=F)
write.table(aplot,quote=F,file=paste(cargs[4],'.cum', sep=''),sep="\t",row.names=F);
#write.table(cbind(T2C,m),quote=F,file=paste(cargs[4],'.spearman.cumcor', sep=''),sep="\t",row.names=F)
#write.table(aplot,quote=F,file=paste(cargs[4],'.spearman.cum', sep=''),sep="\t",row.names=F);

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

#plot.new(); pdf(paste(cargs[4],'.spearman.cumcor.pdf',sep=''), h=8,w=8);
plot.new(); pdf(paste(cargs[4],'.cumcor.pdf',sep=''), h=8,w=8);
levelplot(m, ylab="sample pairs",xlab="T2C",
	par.strip.text = list(cex=1.5,font=2),
	par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1.5,
	scales=list(x=list(rot=90)),
	main=paste("Correlation of ", cargs[5], "-freqs around non-reproduced sites", sep=''))
dev.off();

