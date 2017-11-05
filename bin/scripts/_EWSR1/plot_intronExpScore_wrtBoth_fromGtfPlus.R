getps <- function(x) {
	kk <- x$count / x$len;
	return(min (kk[kk>0])); 	
}

library(plyr)
start.bin <- function(i,e) {
	r <- ddply(e, .(gene_id),function(x) 
		data.frame(gene_start=min(x$start), gene_end=max(x$end)))
	r$gene_len <- r$gene_end - r$gene_start + 1;
	k <- merge(i, r, by=c("gene_id"))
	ki <- k$strand == "+"
	start_bin <- (k$start - k$gene_start + 1) / (k$gene_len)
	start_bin[ki] <- ((k$gene_end - k$end + 1) / (k$gene_len))[ki]
	return (start_bin)
}

erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
#help.search("error function") 
#? stats::Normal

getDataFromGtfPlus <- function (filename, onames) {
	aa <- read.table(gzfile(filename), header=F);

	head(aa)
	im <- 1:8
	io <- 9:dim(aa)[2]

	itp <- io[io %% 3==2]    # temporary array for the indices of punctuation marks
	la <- aa[1,itp] == ";"   # logical array for the attribute

	ns <- aa[1,itp[la]-2]    # names of attribute
	iv <- io[io %% 3 == 1][1:length(la[la])] # values of attributes

	ip <- NULL #in case the file is not gtf plus but gtf
	if ( dim(aa)[2] > itp[la][length(itp[la])] ) {
		ip <- (itp[la][length(itp[la])]+1):(dim(aa)[2]) 
	}

	a <- aa[,c(im,iv,ip)];
	head(a)
	#	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)));
	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)), onames)
	return(a)
}

cargs <- commandArgs(trailingOnly = T)

#cargs <- c("~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/271.introns_coverage.gz", "~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/271.exons_coverage.gz", "count")

tag <- paste(gsub("*(_coverage.gz)", "", cargs[1]), "combined", sep='');

print(tag)

plusnames <- NULL
if (length(cargs) > 2) {
	plusnames <- cargs[3:length(cargs)];
}

x <- getDataFromGtfPlus(cargs[1], plusnames);
y <- getDataFromGtfPlus(cargs[2], plusnames);

x$start_bin <- start.bin(x, y)

x$len <- x$end - x$start + 1;
y$len <- y$end - y$start + 1;

ps <- 1;
ps <- getps (x);
x$count <- x$count + (x$len * ps);
y$count <- y$count + (y$len * ps);

N <- sum(x$count) + sum(y$count);
print(N)

head(x,3)   
head(y,3)   

xx <- ddply(x, .(gene_id), function(x) { data.frame(count=sum(x$count), len=sum(x$end - x$start + 1) ) })
yy <- ddply(y, .(gene_id), function(x) { data.frame(count=sum(x$count), len=sum(x$end - x$start + 1) ) })
head(xx)	
head(yy,3)	

x$mu   <- x$count / x$len / N;
x$var  <- x$mu * (1 - x$mu) / N;

xx$mu  <- xx$count / xx$len / N;
xx$var <- xx$mu * (1 - xx$mu) / N;

yy$mu  <- yy$count / yy$len / N;
yy$var <- yy$mu * (1 - yy$mu) / N;

head(xx,3)
head(yy,3)
#------------
a <- merge(x, merge(xx,yy,by=c("gene_id"), suffixes=c(".introns", ".exons")), by=c("gene_id"))

a$z.exons   <- (a$mu - a$mu.exons) / sqrt(a$var + a$var.exons); 
a$z.introns <- (a$mu - a$mu.introns) / sqrt(a$var + a$var.introns); 

a$erfc.exons   <- erfc( - a$z.exons/sqrt(2)) / 2;
a$erfc.introns <- erfc( - a$z.introns/sqrt(2)) / 2;

a$erfc <- sqrt(a$erfc.exons * a$erfc.introns);

tail(a,4)
dim(a)

fields <- c("gene_id", "seqname","start","end", "strand","len", "count","count.exons","count.introns","z.exons","z.introns","erfc.exons","erfc.introns","erfc")
write.table(a[sort(a$erfc,decreasing=T,index.return=T)$ix,fields], file=paste(tag, ".txt", sep=''), quote=F, row.names=F, sep="\t")

#----
library(lattice)
library(hexbin)
#library(plotrix)

b  <- rbind(data.frame(list(type="z.exons", z=a$z.exons)), data.frame(list(type="z.introns", z=a$z.introns))) 
b2 <- rbind(data.frame(list(type="erfc", erfc=a$erfc)), data.frame(list(type="erfc.exons", erfc=a$erfc.exons)), data.frame(list(type="erfc.introns", erfc=a$erfc.introns))) 

p1 <- densityplot(~ z, groups=type, data=b, plot.points=F, type="percent", auto.key=list(columns=1, corner=c(1,1)), main=tag, ylab = "percent")
p2 <- densityplot(~ erfc, groups=type, data=b2,plot.points=F, type="percent", auto.key=list(columns=1, corner=c(1,1)), main=tag, ylab = "percent")
p3 <- xyplot(erfc.introns ~ erfc.exons, data=a, panel=panel.hexbinplot, main=tag);
p4 <- xyplot(erfc.exons ~ z.exons, data=a, panel=panel.hexbinplot, main=tag);
p5 <- xyplot(erfc.introns ~ z.introns, data=a, panel=panel.hexbinplot, main=tag);
p6 <- xyplot(erfc ~ start_bin, data=a, panel=panel.hexbinplot, main=tag);
p7 <- xyplot(erfc.exons ~ log2(len), data=a, panel=panel.hexbinplot, main=tag);
p8 <- xyplot(erfc.introns ~ log2(len), data=a, panel=panel.hexbinplot, main=tag);
p9 <- xyplot(erfc ~ log2(len), data=a, panel=panel.hexbinplot, main=tag);

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

