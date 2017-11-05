require(limma)
#generalize the output names
Venn <- function(files, n, universe)
{
	stopifnot( n == 3 | n == 2)
	# Form universe as union of all three sets
	Counts <- matrix(0, nrow=length(universe), ncol=n)

	for (j in 1:n) {
		set <- as.vector(t(read.table(files[j], header=F)$V1))
		for (i in 1:length(universe)) {
			Counts[i,j] <- universe[i] %in% set
		}
	}
	return (Counts);
}


#cargs <- system(paste(commandArgs(trailingOnly=T), collapse=" "))
cargs <- commandArgs(trailingOnly=T)
#cargs <- commandArgs(trailingOnly=T);
#cargs <- c('ENCODEHekRep3.INTRON.xlink.gid,ENCODEHekRep3.INTRON.Rloop.gid','ENCODEHekRep3.INTRON.all.gid', "venn", 'EWSR1','Rloop');


#split names ny space
files <- unlist(strsplit(cargs[1],",")); 
n <- length(files)

print(files)

universe <- NULL
if (length(cargs) > 1 && cargs[2] != "") {
	universe <- as.vector(t(read.table(cargs[2], header=F)$V1));
} else {
	for (f in files) {
		universe <- unique(c(universe, as.vector(t(read.table(f, header=F)$V1))))
	}
}

otag <- "venn";
if (length(cargs) > 2 && cargs[3] != "") {
	otag <- cargs[3];
} 

nms <- sub(".genenames", '', files)
if (length(cargs) > 3) {
	nms <- c(cargs[4:(3+n)]);
}
print(nms)
print(nms[0])


print(c(files, n, nms, length(universe), otag));
Counts <- Venn(files, n, universe);
colnames(Counts) <- nms

plot.new(); pdf(paste(otag,'.pdf',sep=''), h=8, w=8)
#	vennDiagram( vennCounts(Counts, include="both"), cex=1.5, font=2, lwd=2, mar=c(0,0,0,0), fontfamily="sans", circle.col=c("red","green","blue") )
vennDiagram( vennCounts(Counts, include="both"), cex=1.5, font=2, lwd=2)
dev.off();

write.table(data.frame(universe,Counts),file=paste(otag,'.txt', sep=""),row.names=F, quote=F, sep="\t");

