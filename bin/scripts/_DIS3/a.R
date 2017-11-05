require(limma)
#generalize the output names
Venn3 <- function(set1, set2, set3, names, universe, otag)
{
	stopifnot( length(names) == 3)
	# Form universe as union of all three sets
	Counts <- matrix(0, nrow=length(universe), ncol=3)
	colnames(Counts) <- names

	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
		Counts[i,3] <- universe[i] %in% set3
	}
	vennDiagram( vennCounts(Counts, include="both"), cex=1.5, font=2, lwd=2 )

	write.table(data.frame(universe,Counts),file=paste(otag,'.txt', sep=""),row.names=F, quote=F, sep="\t");

}


cargs <- commandArgs(trailingOnly=T);
#cargs <- c("DIS3L2.genenames","DIS3L.genenames","DIS3.genenames",".genenames","clipz.genenames","out","DIS3L2","DIS3L","DIS3")

if (length(cargs) > 6 && cargs[5] != "") {
#	universe <- as.vector(t(read.table(cargs[5], header=F)$V1));
} 
otag <- "venn3";
if (length(cargs) > 7 && cargs[6] != "") {
	otag <- cargs[6];
} 
names <- sub(cargs[4], '', c(basename(cargs[1]),basename(cargs[2]), basename(cargs[3])))
if (length(cargs) > 8) {
	names <- c(cargs[7:9]);
}
print(names)
print(names[1])




