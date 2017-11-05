cargs <- commandArgs()

cargs <- c(1,1,1,"~bilebi00/_DIS3/data/clipz_RNAseq/mRNAseq.raw", '.mclustexp');
ext <- cargs[5];
if (length(cargs) < 5) {
	ext <- '.mclustExp';
}

library(mclust);

getpc <- function(x,cs) { 
	colmin <- 1:cs * 0;
	for (c in 1:cs) {
		colmin[c] <- sort(x[x[,c]>0,c])[1] #get min that is greater than 0
	};
	return(colmin);
}; 

aa <- read.table(cargs[4], header=T); #for writing with row.names column header
a <- read.table(cargs[4], header=T, row.names="id");
rs <- dim(a)[1]
cs <- dim(a)[2] 
mins <- getpc(a,cs);
d <- a + matrix(rep(mins,rs),nrow=rs,byrow=T)

classification <- matrix(rep(0,rs*cs),nrow=rs,byrow=T);
i=4
for(i in 1:cs) { 
	paste(i)
	fgcoverage <- sample(d[,i], 10000);
	clustRes <- Mclust(log2(fgcoverage), G=2, modelNames="V");
	means <- clustRes$parameters$mean;
	vars <- clustRes$parameters$variance$sigmasq;
	ra <- range(log2(fgcoverage));
	indxfg <- ((dnorm(log2(d[,i]), means[2], sqrt(vars[2]))>dnorm(log2(d[,i]), means[1], sqrt(vars[1]))) & (log2(d[,i])>means[1]))
	classification[,i] <- as.integer(indxfg); 

	#plot the distributions
	n=1000;
	plot.new(); pdf(paste(cargs[4], names(a)[i],'.norm.and.classification.pdf',sep=''),h=8,w=8);
	par(mfrow=c(1,2));
	plot(density(log2(d[,i])), main=names(a)[i], ylim=c(0, 0.5)); 
	lines(density(rnorm(n, means[1], sqrt(vars[1]))), col="green");
	lines(density(rnorm(n, means[2], sqrt(vars[2]))), col="red");
	#plot the classification
	plot(density(log2(d[,i])), main=names(a)[i], ylim=c(0, 0.5)); 
	lines(density(log2(d[indxfg,i])), col="red");
	lines(density(log2(d[!indxfg,i])), col="green");
	dev.off();
}

indxExpressed <- apply(classification,1,sum)>0

#Classification=as.integer(indxfg);
outfile=paste(cargs[4], ext, sep="");
write.table(aa[indxExpressed,], file=outfile, sep="\t", row.names=F, quote=F)

