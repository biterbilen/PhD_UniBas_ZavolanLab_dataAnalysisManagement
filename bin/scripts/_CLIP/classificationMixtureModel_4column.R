cargs <- c("1344.exp", '.mclustexp',1);
cargs <- c("mRNAseq_siUntreatedHeLa.exp", '.mclustexp',1);
cargs <- commandArgs(trailingOnly=T)
ext <- '.mclustExp';
if (length(cargs) > 1) {
	ext <- cargs[2];
}

#fits
library(mclust);

header <- c("Transcript","Length","Counts")
#a <- read.table(cargs[1], header=F, col.names=header)
a <- read.table(cargs[1], header=F,sep="\t")
names(a) <- header;
a <- a[a$Counts>0,]
head(a)
#add uniform ps wrt length
a$Expression <- a$Counts / a$Length

#20120401 fitted gaussians were different each time so for consistency; we use whole data
#fgcoverage=sample(a$Expression, 5000);
fgcoverage=sample(a$Expression, dim(a)[1])
clustRes <- Mclust(log2(fgcoverage), G=2, modelNames="V");
summary(clustRes)
means <- clustRes$parameters$mean;
vars <- clustRes$parameters$variance$sigmasq;
ra <- range(log2(fgcoverage));
#20120401 smaller mean component had higher value members than bigger value component
#indxfg=((dnorm(log2(a$Expression), means[2], sqrt(vars[2]))>dnorm(log2(a$Expression), means[1], sqrt(vars[1]))) & (log2(a$Expression)>means[1]))
indxfg <- ((dnorm(log2(a$Expression), means[2], sqrt(vars[2]))>dnorm(log2(a$Expression), means[1], sqrt(vars[1]))) & (log2(a$Expression)>means[1]) | 
				   (dnorm(log2(a$Expression), means[2], sqrt(vars[2]))<dnorm(log2(a$Expression), means[1], sqrt(vars[1]))) & (log2(a$Expression)>means[2]))
#TODO solve the gaussian to set the classification line 
#plot the distributions
n <- 1000;

plot.new(); pdf(paste(cargs[1],'.classification.pdf',sep=''),h=8,w=8);
par(mfrow=c(1,2));
plot(density(log2(a$Expression)), main=basename(cargs[1]), ylim=c(0, 0.5)); 
lines(density(rnorm(n, means[1], sqrt(vars[1]))), col="green");
lines(density(rnorm(n, means[2], sqrt(vars[2]))), col="red");
#plot the classification
plot(density(log2(a$Expression)), main=basename(cargs[1]), ylim=c(0, 0.5)); 
lines(density(log2(a$Expression[indxfg])), col="red");
lines(density(log2(a$Expression[!indxfg])), col="green");
dev.off();
  
#Classification=as.integer(indxfg);
cat(sum(indxfg), "clusters in foreground", "\n");
outfile=paste(cargs[1], ext, sep="");
a$Expression <- NULL
write.table(a[indxfg,], file=outfile, sep="\t", row.names=F, quote=F)

