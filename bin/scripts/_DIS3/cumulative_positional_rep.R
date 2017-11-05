library(latticeExtra)
cargs <- commandArgs()
#cargs <- c(1,1,1,"~/_DIS3/Analysis/aa/posstat");
a <- read.table(cargs[4], header=T);

T2Cs <- sort(unique(a$T2C), decreasing=T);
proteinNames <- sort(unique(a$proteinName), decreasing=T);

aplot <- NULL; #initialize
k <- 1;
for (t2c in T2Cs) {
	suba <- a[a$T2C >= t2c,];
	yy <- aggregate(suba[,c(4:5)], by=list(suba$proteinName, suba$sampleId), sum);
	yy$T2C <- t2c;
	yy$proteinName <- yy$Group.1;
	yy$sampleId <- yy$Group.2;
	yy$Group.1 <- NULL;
	yy$Group.2 <- NULL;
	aplot <- rbind(aplot,yy)
}

write.table(aplot,quote=F,file=paste(cargs[4],'.cum', sep=''),sep="\t",row.names=F);


ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

p <- proteinNames[1]
for (p in proteinNames) {
	a <- aplot[aplot$proteinName==p,];
	print(dim(a))
	tit = paste("Fraction of reproduced T2C positions ", p, sep="-");
	ofilename <- paste(cargs[4], "-", p, ".cum.pdf", sep="")
	plot.new(); pdf(ofilename, width=8, height=8)
	plt <- levelplot(reproducedT2Cpositions / totalT2Cpositions ~ log2(T2C) + log2(totalT2Cpositions) 
		| sampleId, data=a, main=tit,
		auto.key = list(space = "top", title="Sample Id"),
		par.strip.text = list(cex=1.5,font=2),
		par.settings=list(axis.text=list(font=2,cex=1.5),
			par.ylab.text=list(font=2,cex=1.5),
			par.xlab.text=list(font=2,cex=1.5),
			par.main.text=list(font=2,cex=1.5)),
		lwd=2,cex=2.5,
		panel = panel.levelplot.points, type=c("g","p","l") )
	print(plt)
	dev.off()
}

#levelplot(reproducedT2Cpositions / totalT2Cpositions ~ log2(T2C) + log2(totalT2Cpositions) 
#	| paste(sampleId), data=a, main=tit, cex=2, groups=sampleId, allow.multiple=T,
#	panel = panel.levelplot.points, type=c("g","p","l"))

