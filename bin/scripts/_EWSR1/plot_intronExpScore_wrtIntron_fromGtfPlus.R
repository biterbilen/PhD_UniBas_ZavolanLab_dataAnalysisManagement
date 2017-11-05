source(file.path("/import/bc2/home/zavolan/bilebi00/_EWSR1/scripts/","getDataFromGtfPlus.R"))
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))


getps <- function(x) {
	kk <- x$count / x$len;
	return(min (kk[kk>0])); 	
}

library(plyr)
get.bin <- function(i,e, strand) {
	r <- ddply(e, .(gene_id),function(x) 
		data.frame(gene_start=min(x$start), gene_end=max(x$end)))
	r$gene_len <- r$gene_end - r$gene_start + 1;
	k <- merge(i, r, by=c("gene_id"))
	#strand = "+" start_bin
	#strand = "-" end_bin
	ki <- k$strand == strand
	bin <- (k$start - k$gene_start + 1) / (k$gene_len)
	bin[ki] <- ((k$gene_end - k$end + 1) / (k$gene_len))[ki]
	return (bin)
}

erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
#help.search("error function") 
#? stats::Normal

cargs <- commandArgs(trailingOnly = T)

#cargs <- c("~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/271.introns_coverage.gz", "~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/271.exons_coverage.gz", "count")

tag <- paste(gsub("*(.coverage.gz)", "", cargs[1]), sep='');
titel <- cargs[3];

print(tag)

plusnames <- NULL
if (length(cargs) > 3) {
	plusnames <- cargs[4:length(cargs)];
}

x <- getDataFromGtfPlus(cargs[1], plusnames);
y <- getDataFromGtfPlus(cargs[2], NULL);

x$start_bin <- get.bin(x, y, "+")
x$end_bin   <- get.bin(x, y, "-")

x$len <- x$end - x$start + 1;
head(x)
ps <- 1;
ps <- getps (x);
x$count <- x$count + (x$len * ps);

N <- sum(x$count);
print(N)

head(x,3)   

xx <- ddply(x, .(gene_id), function(x) { data.frame(count=sum(x$count), len=sum(x$end - x$start + 1) ) })
head(xx)	

x$mu   <- x$count / x$len / N;
x$var  <- x$mu * (1 - x$mu) / N;

xx$mu  <- xx$count / xx$len / N;
xx$var <- xx$mu * (1 - xx$mu) / N;

head(xx,3)
#------------
a <- merge(x, xx, by=c("gene_id"), suffixes=c("", ".introns"))

a$z.introns <- (a$mu - a$mu.introns) / sqrt(a$var + a$var.introns); 

a$erfc.introns <- erfc( - a$z.introns/sqrt(2)) / 2;

tail(a,4)
dim(a)

#fields <- c("gene_id", "seqname","start","end", "strand","len", "count","count.introns","z.introns","erfc.introns")
fields <- c("seqname","start","end", "gene_id", "erfc.introns", "strand","len", "count","count.introns","z.introns")
write.table(a[sort(a$erfc.introns,decreasing=T,index.return=T)$ix,fields], file=paste(tag, ".txt", sep=''), quote=F, row.names=F, sep="\t")

#----
library(lattice)
library(hexbin)
#library(plotrix)

b  <- rbind(data.frame(list(type="z.introns", z=a$z.introns))) 
b2 <- rbind(data.frame(list(type="erfc.introns", erfc=a$erfc.introns))) 


plot.new(); pdf(paste(tag, "_params.pdf", sep=''), h=8, w=8);
#densityplot(~ z, groups=type, data=b, plot.points=F, type="percent", auto.key=list(columns=1, corner=c(1,1)), main=tag, ylab = "percent")
#xyplot(erfc.introns ~ start_bin, data=a, panel=panel.hexbinplot, main=tag);
#xyplot(erfc.introns ~ end_bin, data=a, panel=panel.hexbinplot, main=tag);
#densityplot(~ erfc, groups=type, data=b2,plot.points=F, type="percent", auto.key=list(columns=1, corner=c(1,1)), main=tag, ylab = "percent")
#xyplot(erfc.introns ~ log2(len), data=a, panel=panel.hexbinplot, main=tag);
#xyplot(erfc.introns ~ z.introns, data=a, panel=panel.hexbinplot, main=tag);
densityplot(~ z, groups=type, data=b, plot.points=F, type="percent", main=titel, ylab = "percent", par.settings = article.theme,font=article.font,cex=article.cex,lwd=article.lwd)
xyplot(erfc.introns ~ start_bin, data=a, panel=panel.hexbinplot, main=titel, par.settings = article.theme,font=article.font,cex=article.cex,lwd=article.lwd);
xyplot(erfc.introns ~ end_bin, data=a, panel=panel.hexbinplot, main=titel, par.settings = article.theme, font=article.font,cex=article.cex,lwd=article.lwd);
densityplot(~ erfc, groups=type, data=b2,plot.points=F, type="percent", main=titel, ylab = "percent",par.settings = article.theme, font=article.font,cex=article.cex,lwd=article.lwd)
xyplot(erfc.introns ~ log2(len), data=a, panel=panel.hexbinplot, main=titel,par.settings = article.theme, font=article.font,cex=article.cex,lwd=article.lwd);
xyplot(erfc.introns ~ z.introns, data=a, panel=panel.hexbinplot, main=titel,par.settings = article.theme, font=article.font,cex=article.cex,lwd=article.lwd);
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

