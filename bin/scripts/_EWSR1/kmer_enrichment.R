getstats <- function(a, colpat, lencut) {
	#incomplete beta function hyperparameters
	aa <- 1
	bb <- 1

	head(a,2)
	#filter based on length
	a$seq_len <- a$end - a$start + 1 
	b <- a[a$seq_len > lencut,];

	#clipped
	tmp <- b[!is.na(b$maxpbeta),];

	d <- NULL;
	d$name <- names(b)[grep(colpat,names(b))];
	d$len <- lencut
	N <- dim(tmp)[1]
	d$N <- N
	d$kmerall <- sum(tmp$kmer.all);
	d$kmer <- colSums(tmp[,d$name])
	D <- as.data.frame(d)

	#not clipped introns of genes having at least one clipped intron
	tmp2 <- NULL
	tmp2$gene_id <- unique(tmp$gene_id)

	tmp <- subset(merge(b,tmp2), is.na(maxpbeta));
	head(tmp,2)
	D$NnotClipped <- dim(tmp)[1]
	D[paste("notclipped","background_freq",sep='.')] <- colSums(tmp[,d$name]) / sum(tmp$kmer.all)
	D[paste("notclipped","logpbeta",sep='.')] <- pbeta( D[,paste("notclipped","background_freq",sep='.'),1], D$kmer + aa, D$kmerall - D$kmer + bb, log=T);

#	D
#	head(D,2)
	#top expressed of different libraries
	for (lib in names(b)[grep(".erfc",names(b))]) {
		i <- sort(b[,c(lib)], decreasing=T, index.return = T)
		tmp <- b[i$ix,][1:N,]
		D[paste(lib,"background_freq",sep='.')] <- colSums(tmp[,d$name]) / sum(tmp$kmer.all)
		D[paste(lib,"logpbeta",sep='.')] <- pbeta( D[,paste(lib,"background_freq",sep='.'),1], D$kmer + aa, D$kmerall - D$kmer + bb, log=T);
	}
	return(D)
}


getDataFromGtfPlus <- function (filename, onames) {

	onames
	aa <- read.table(gzfile(filename), header=F);

	im <- 1:8
	io <- 9:dim(aa)[2]

	itp <- io[io %% 3 == 2]    # temporary array for the indices of punctuation marks
	la <- aa[1,itp] == ";"     # logical array for the attribute

	ns <- aa[1,itp[la]-2]      # names of attribute
	iv <- io[io %% 3 == 1][1:length(la[la])] # values of attributes

	ip <- NULL #in case the file is not gtf plus but gtf
	if ( dim(aa)[2] > itp[la][length(itp[la])] ) {
		ip <- (itp[la][length(itp[la])]+1):(dim(aa)[2]) 
	}

	a <- aa[,c(im,iv,ip)];
	#	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)));
	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)), onames)
	head(a,2)
	return(a)
}

cargs <- commandArgs(trailingOnly = T)

#cargs <- c("all.nucbed.1mers", "hg19_GMAP_GENEintrons.refseq.gtf.nucbed.werfc.wclip", "all.nucbed.1mers.onames");

plusnames <- NULL
if (length(cargs) > 2) {
	plusnames <- names(read.table(cargs[3], header=T, sep="\t"));
} 

#TODO hard coded here
plusnames2 <- c("pct_at","pct_gc","num_A","num_C","num_G","num_T","num_N","num_oth","seq_len","user_patt _count","siCTRL_1.erfc","siEWSR1_1.erfc","siCTRL_2.erfc","siEWSR1_2.erfc","undiff_osteo_P5.erfc","20daydiff_osteo_P5.erfc","maxpbeta");
if (length(cargs) > 3) {
	plusnames2 <- cargs[4:length(cargs)];
} 

x <- getDataFromGtfPlus(cargs[1], plusnames);
head(x)
y <- getDataFromGtfPlus(cargs[2], plusnames2);
a <- merge(x,y,by=c("gene_id","seqname", "start", "end", "strand"))
print(dim(a))
head(a,2)

#--------------
D <- NULL;
for (len in c(0,100,500,1000,2000,4000)) {
	print(len)
	D <- rbind(D,getstats(a, "kmer", len));
}
dim(D)

D$maxlogpbeta <- apply(D[grep(".logpbeta",names(D),fixed=T)], 1, max)
s <- sort(D$maxlogpbeta,index.return = T);
head(D[s$ix,],2)
write.table(format(D[s$ix,],digits=2), row.names=F, quote=F, sep="\t", file=paste(cargs[1],'.logpbeta', sep=''));
quit();

library(lattice)
plot.new(); pdf(paste(cargs[1],'.num.pdf',sep=''),h=8,w=8);
barchart( name ~ mu | paste("N=", N, " length cutoff=",sort(as.factor(len)),sep=""), data=D, groups=tag, auto.key=list(space="top"))
dev.off();

quit();

#--------------
library(plotrix)
library(hexbin)
plotCI(barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4)),a$mu,a$se,pch=21,xlim=c(0,l))
plot.new(); pdf(paste(cargs[1],'.pdf',sep=''), h=8, w=8);
p10 <- plotCI(barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4)),a$mu,a$se,2*a$se,add=TRUE,pch=21,xlim=c(0,l))
text(xp,(a$mu+5*max(a$se)), a$tag, srt = 90, adj=c(0,0.5))


xs <- 3;
ys <- 3;
plot.new(); pdf(paste(tag, "_params.pdf", sep=''), h=8, w=8);
plot(p1, split = c(1, 1, xs, ys))
plot(p2, split = c(1, 2, xs, ys), newpage = FALSE)
plot(p3, split = c(1, 3, xs, ys), newpage = FALSE)
plot(p4, split = c(2, 1, xs, ys), newpage = FALSE)
plot(p5, split = c(2, 2, xs, ys), newpage = FALSE)
plot(p6, split = c(2, 3, xs, ys), newpage = FALSE)
plot(p7, split = c(3, 1, xs, ys), newpage = FALSE)
plot(p8, split = c(3, 2, xs, ys), newpage = FALSE)
plot(p9, split = c(3, 3, xs, ys), newpage = FALSE)
dev.off();

quit()
#----
plotCI(x$median,x$mean,a$stdev,pt.bg=par("bg"),pch=21,
	main="plotCI with extra space on the x axis")
plotCI(x$median,x$mean,a$stdev,pt.bg=par("bg"),pch=21,err="x",add=T)
mtext("for adding horizontal error bars",3,0.5)
unique(a$median)
library(hexbin)

xyplot(log2(median) ~ log2(mean), data=x, panel=panel.hexbinplot)
densityplot(~mean, data= x, plot.points=F, nint=100)

xecdf <- ecdf(x$mean)
r <- range(x$mean)
curve(1-xecdf(x), from=r[1], to=r[2], col="red", xlim=r, ylim=0:1)
#scale data, axis and their combinations

library(plotrix)

#cargs <- c(1,1,1,"all.stat")
cargs <- commandArgs(trailingOnly = T)

a <- read.table(cargs[1], header=T)
head(a)

a$se <- a$sd / sqrt(a$N);

xp <- barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4))
l <- length(xp)

plotCI(barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4)),a$mu,a$se,pch=21,xlim=c(0,l))
plot.new(); pdf(paste(cargs[1],'.pdf',sep=''), h=8, w=8);
plotCI(barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4)),a$mu,a$se,2*a$se,add=TRUE,pch=21,xlim=c(0,l))
text(xp,(a$mu+5*max(a$se)), a$tag, srt = 90, adj=c(0,0.5))
dev.off()



