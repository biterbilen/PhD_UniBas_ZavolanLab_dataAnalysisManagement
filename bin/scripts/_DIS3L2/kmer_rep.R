library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("4mers288.backgroundFile_mu_ref.pbeta", "4mers984.backgroundFile_mu_ref.pbeta", "aa", "8", "288", "984");

topn <- 5;
if (length(cargs)>3) {
	topn = as.numeric(cargs[4]);
}
name1 <- "1"
if (length(cargs)>4) {
	name1 = cargs[5];
}
name2 <- "2"
if (length(cargs)>5) {
	name2 = cargs[6];
}

a <- read.table(cargs[1],header=T)
b <- read.table(cargs[2],header=T)

k <- merge(a,b,by="kmer",suffixes=c(".1",".2"));
dim(k)
k$lab1 <- "";
tx <- tail(head(sort(k$logpbeta.1, decreasing=F),topn),1)
ty <- tail(head(sort(k$logpbeta.2, decreasing=F),topn),1)
lis <- (k$logpbeta.1 <= tx) | (k$logpbeta.2 <= ty)
k$lab1[lis] <- gsub("kmer.","",k$kmer)[lis];

#plot.new(); pdf(paste(cargs[3],'.pbeta.pdf',sep=''),h=8,w=8);
p1 <- xyplot(-logpbeta.2 ~ -logpbeta.1, data=k,
	ylab=paste(name2, "kmer Enrichment"),
	xlab=paste(name1, "kmer Enrichment"),
	auto.key=list(corner=c(1,1), cex=1.5,font=2),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	lwd=2,cex=1.5,font=2,pch=".", 
	) +
layer(panel.key(c(sprintf('N=%d\tR=%.2f', length(x),cor(x,y, method="spearman") )),col='darkgray',points=F,cex=1.5,font=2)) +
#layer(panel.abline(lm(y ~ x), lwd=2, col = 'black')) + 
layer(panel.abline(a=0, b=1, col="darkgray",lwd=2)) +
layer(panel.text(..., labels=k$lab1, cex=1, font=2, adj=1))

#plot(p1)
#dev.off();

k$lab <- "";
tx <- tail(head(sort(k$freqs.1, decreasing=T),topn),1)
ty <- tail(head(sort(k$freqs.2, decreasing=T),topn),1)
lis <- (k$freqs.1 >= tx | k$freqs.2 >= ty)
k$lab[lis] <- gsub("kmer.","",k$kmer)[lis];

p2 <- xyplot(freqs.2 ~ freqs.1, data=k,
	ylab=paste(name2, "kmer Count"),
	xlab=paste(name1, "kmer Count"),
	auto.key=list(corner=c(1,1), cex=1.5,font=2),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	lwd=2,cex=1.5,font=2,pch=".", 
	) +
layer(panel.key(c(sprintf('N=%d\tR=%.2f', length(x),cor(x,y, method="spearman") )),col='darkgray',points=F,cex=1.5,font=2)) +
#layer(panel.abline(lm(y ~ x), lwd=2, col = 'black')) + 
layer(panel.abline(a=0, b=1, col="darkgray",lwd=2)) +
layer(panel.text(..., labels=k$lab, cex=1, font=2, adj=1))

p3 <- densityplot( ~freqs.1, data=k,
	xlab=paste(name1, "kmer Count"),
	par.settings = poster.theme,
	lwd=2,cex=1.5,font=2,pch=".")
p4 <- densityplot( ~freqs.2, data=k,
	xlab=paste(name2,"kmer Count"),
	par.settings = poster.theme,
	lwd=2,cex=1.5,font=2,pch=".")
plot.new(); pdf(paste(cargs[3],'.pdf',sep=''),h=8,w=8);
#plot(p1, split = c(1, 1, 2, 2)) 
#plot(p2, split = c(2, 1, 2, 2), newpage = FALSE) 
plot(p1)
plot(p2)
plot(p3)
plot(p4)
#title(main=cargs[1])

#update(c(p1,p2,p3,p4), main=cargs[3]);
dev.off();

quit()

plot(xyplot(1:10~1:10, auto.key=list(title="gdgdg",columns=1,space="top")), split = c(1, 1, 1, 2)) 
plot(xyplot(10:20~10:20), split = c(1, 2, 1, 2), newpage = FALSE) 

