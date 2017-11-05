source(file.path("/import/bc2/home/zavolan/bilebi00/_EWSR1/scripts/","getDataFromGtfPlus.R"))

cargs <- commandArgs(trailingOnly = T)
#cargs <- c("6mers", "6mers.onames", "6mers.background")

plusnames <- NULL
if (length(cargs) > 1) {
	plusnames <- levels(unlist(read.table(cargs[2], header=F)))
}

#uniform pseudocount
aa <- 1
bb <- 4**(nchar(gsub("kmer.","",plusnames[1]))) - aa;

a <- getDataFromGtfPlus(cargs[1], plusnames);
N <- colSums(a["kmer.all"])


#uniform background for mu_ref 
D <- NULL
D$kmer <- names(a)[grep("kmer.[ACGT]",names(a), perl=T)];
D$freqs <- colSums(a[,D$kmer]);
D$N <- N
d <- data.frame(D)
D <- d

tag <- NULL
head (D)
if (length(cargs) > 2) {
	b <- getDataFromGtfPlus(cargs[3], plusnames);
	print(head (D))
	D$mu_ref <- colSums(b[,levels(unlist(D$kmer))]) / colSums(b["kmer.all"]);
	D$mu_ref[D$mu_ref == 0] <- min(D$mu_ref[D$mu_ref > 0]);
	head (D)
	tag <- ".backgroundFile_mu_ref";
} else {

	D$mu_ref <- aa/(bb+aa)
	D$logpbeta <- pbeta(D$mu_ref, D$freq + aa, D$N - D$freq + bb, log.p=T)
	tag <- ".uniform_mu_ref";
}

D$logpbeta <- pbeta(D$mu_ref, D$freq + aa, D$N - D$freq + bb, log.p=T)
D$regionCount <- dim(a)[1]
write.table(D,row.names=F,sep="\t",quote=F, file=paste(cargs[1],tag,".pbeta", sep=''))

