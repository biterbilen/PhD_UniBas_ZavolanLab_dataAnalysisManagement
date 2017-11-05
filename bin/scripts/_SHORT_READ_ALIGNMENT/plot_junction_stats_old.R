cargs <- commandArgs(); 

#cargs <- c(1,1,1,'../aaa','Tophat');
#cargs <- c(1,1,1,'../JUNC.inp.stat','Tophat');
#cargs <- c(1,1,1,"../rep_A.hg18_chr19/MapSplice/best_junction.bed.GMAP.stat", 'MapSplice');
a <- read.table(cargs[4], header=T, check.names=F);


library("latticeExtra")
tit = cargs[5];
#tit = paste("Overlap of ", cargs[5], " and GMAP splice junctions", sep="")
ofilename <- paste(cargs[4],".pdf", sep="")
plot.new(); pdf(ofilename, width=8, height=8)
nop <- xyplot(hitJunctionCount ~ log2(coverage) | overlapfrac, data=a,
#nop <- xyplot(hitJunctionCount ~ coverage | overlapfrac, data=a,
	par.strip.text=list(cex=1.5), 
	par.settings=list(axis.text=list(cex=1.5)),
	xlab=list(cex=2), ylab=list(cex=2),
	main=list(label=tit, cex=2),
	cex=2,
	lcw=c(2),
	groups = overlapfrac, auto.key=list(title=tit), type=c("p"));
forp <- xyplot(hitJunctionCount/allJunctionCount ~ log2(coverage) | overlapfrac, data=a,
#forp <- xyplot(hitJunctionCount/allJunctionCount ~ coverage | overlapfrac, data=a,
	lcw=c(2),
	groups = overlapfrac, type=c("h"));
doubleYScale(forp, nop, style1 = 0, style2 = 3, add.ylab2 = T)
dev.off();


