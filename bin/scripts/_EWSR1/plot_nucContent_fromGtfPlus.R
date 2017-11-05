getstats <- function(a, colpat, lencut) {

	#filter based on length
	b <- a[a$end - a$start + 1 > lencut,];

	dim(b)
	#clipped
	d <- NULL;
	d$name <- names(b)[grep(colpat,names(b))];
	d$len <- lencut
	d$tag <- "clipped";
	tmp <- b[!is.na(b$maxpbeta),];
	len <- matrix(rep(tmp$seq_len,each=length(d$name)),ncol=length(d$name),byrow=T);
	N <- dim(tmp)[1]
	d$N <- N
	d$mu <- mean(tmp[,d$name]/len)
	d$se <- sd(tmp[,d$name]/len)/sqrt(N)
	D <- as.data.frame(d)

	#top expressed of different libraries
	for (lib in names(b)[grep(".erfc",names(b))]) {
		i <- sort(b[,c(lib)], decreasing=T, index.return = T)
		tmp <- b[i$ix,][1:N,]
		len <- matrix(rep(tmp$seq_len,each=length(d$name)),ncol=length(d$name),byrow=T);
		d$tag <- gsub(".erfc","",lib)
		d$mu <- mean(tmp[,d$name]/len)
		d$se <- sd(tmp[,d$name]/len)/sqrt(N)
		D <- rbind(D,as.data.frame(d))
	}
	return(D)
}


getDataFromGtfPlus <- function (filename, onames) {
	aa <- read.table(gzfile(filename), header=F);

	head(aa,2)
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

#cargs <- c("hg19_GMAP_GENEintrons.refseq.gtf.nucbed.werfc.wclip","pct_at","pct_gc","num_A","num_C","num_G","num_T","num_N","num_oth","seq_len","user_patt_count","siCTRL_1.erfc","siEWSR1_1.erfc","siCTRL_2.erfc","siEWSR1_2.erfc","undiff_osteo_P5.erfc","20daydiff_osteo_P5.erfc","maxpbeta")

plusnames <- NULL
if (length(cargs) > 1) {
	plusnames <- cargs[2:length(cargs)];
}

x <- getDataFromGtfPlus(cargs[1], plusnames);
head(x,1);


#--------------
D <- NULL;
for (len in c(0,100,500,1000,2000,4000)) {
	print(len)
	D <- rbind(D,getstats(x, "num_[ACGT]", len));
}
head(D,2)
dim(D)
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



