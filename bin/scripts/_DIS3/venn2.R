require(limma)
#generalize the output names
Venn2 <- function(set1, set2, names, universe)
{
	stopifnot( length(names) == 2)
	# Form universe as union of all three sets
	Counts <- matrix(0, nrow=length(universe), ncol=2)
	colnames(Counts) <- names

	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
	}
	vennDiagram( vennCounts(Counts, include="both"), cex=2, font=2, lwd=2 )

	write.table(data.frame(universe,Counts),file='venn2.txt',row.names=F, quote=F, sep="\t");

}


cargs <- commandArgs();
#cargs <- c(1,1,1,"DIS3L2.genenames","DIS3.genenames",".genenames","clipz.genenames")

set1 <- as.vector(t(read.table(cargs[4], header=F)));
set2 <- as.vector(t(read.table(cargs[5], header=F)));
if (length(cargs) < 7) {
	universe <- unique( c(set1, set2) ) 
} else {
	universe <- as.vector(t(read.table(cargs[7], header=F)));
}

names <- sub(cargs[6], '', c(basename(cargs[4]),basename(cargs[5])))

plot.new(); pdf("venn2.pdf", h=8, w=8)
Venn2(set1, set2, names, universe);
dev.off();


