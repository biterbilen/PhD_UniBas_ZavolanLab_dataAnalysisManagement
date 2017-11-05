library(plotrix)

cargs <- commandArgs(trailingOnly=T);
#cargs <- c("out","pat","universe.txt","a", "b", "c");

otag <- cargs[1];
suff <- cargs[2];
universe <- unique(as.vector(t(read.table(cargs[3], header=F)$V1)))
index <- 4;

nsets <- length(cargs) - (index - 1);
nms <- NULL;

Counts <- matrix(0, nrow=length(universe), ncol=nsets)

for (i in index:length(cargs)) {
	nms <- c(nms,sub(suff, '', cargs[i]));
	set <- as.vector(t(read.table(cargs[i], header=F)$V1))
	ci <- i-index+1
	print(ci)
	print(nms)
	for (j in 1:length(universe))
	{
		Counts[j,ci] <- universe[j] %in% set;
	}
}
colnames(Counts) <- nms

Countslist <- makeIntersectList(Counts)

plot.new(); pdf(paste(otag,'.pdf',sep=''), h=24, w=24)
intersectDiagram(Countslist, show.nulls=T, pct=F, main=otag)
dev.off()
				      
write.table(data.frame(universe,Counts),file=paste(otag,'.tab', sep=""),row.names=F, quote=F, sep="\t");
