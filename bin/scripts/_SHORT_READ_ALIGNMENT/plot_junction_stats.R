cargs <- commandArgs(); 

#cargs <- c(1,1,1,'../aaa','Tophat');
#cargs <- c(1,1,1,'../JUNC.inp.stat','Tophat');
#cargs <- c(1,1,1,"../rep_A.hg18_chr19/tophat/junctions.bed.GMAP.stat", 'TopHat');
a <- read.table(cargs[4], header=T, check.names=F);

library(lattice)
tit = cargs[5];
#tit = paste("GMAP SJs junctions hit by ", cargs[5], sep="")
ofilename <- paste(cargs[4],".pdf", sep="")
#a$coverageSubgroups <- cut(a$coverage, 2);
a$coverageSubgroups <- cut(a$coverage, breaks=c(0,1,max(a$coverage)), include.lowest=T, right=F);
b <- aggregate(a$hitJunctionCount, by=list(a$coverageSubgroups, a$overlapfrac), sum)
names(b) <- c("coverageSubgroups", "overlapfrac", "hitJunctionCount");
plot.new(); pdf(ofilename, width=8, height=8)
barchart(hitJunctionCount~coverageSubgroups|overlapfrac, data=b,
	par.strip.text = list(cex=1.0,abbreviate = TRUE, minlength = 5),
	stack=T,
	par.settings=list(axis.text=list(cex=1.0)),
	xlab=list(cex=2), ylab=list(cex=2),
	main=list(label=tit, cex=2),
	cex=1.0,
#	lcw=c(2),
#	layout=c(1,4),
	groups = overlapfrac, 
	auto.key=list(title="overlapfrac",space='top'), 
	type=c("p"),
	panel = function(...) {
		panel.grid(v = -1, h = -1)
		panel.barchart(...)
	})

dev.off();
 

