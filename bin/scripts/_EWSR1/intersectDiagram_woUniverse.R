library(plotrix)

cargs <- commandArgs(trailingOnly=T);
#cargs <- c(3, "a", "b", "c", "out","");

nsets <- as.numeric(cargs[1]);
universe <- NULL;
for (i in 1:nsets) {
	universe <- c(universe,as.vector(t(read.table(cargs[i+1], header=F)$V1)))
}
universe <- unique(universe)

otag <- cargs[1+nsets+1];
reppat <- "";
if (length(cargs) > 1+nsets+1) {
	reppat <- cargs[1+nsets+2];
}
nms <- NULL;

Counts <- matrix(0, nrow=length(universe), ncol=nsets)

for (i in 1:nsets) {
	nms <- c(nms,sub(reppat, '', cargs[i+1]));
	set <- as.vector(t(read.table(cargs[i+1], header=F)$V1))
	for (j in 1:length(universe))
	{
		Counts[j,i] <- universe[j] %in% set;
	}
}
colnames(Counts) <- nms

Countslist <- makeIntersectList(Counts)

plot.new(); pdf(paste(otag,'.pdf',sep=''), h=12, w=12)
intersectDiagram(Countslist, show.nulls=T, pct=F, main=otag)
dev.off()
				      
write.table(data.frame(universe,Counts),file=paste(otag,'.txt', sep=""),row.names=F, quote=F, sep="\t");
