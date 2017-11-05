source(file.path("/import/bc2/home/zavolan/bilebi00/_EWSR1/scripts/","getDataFromGtfPlus.R"))
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly = T)
#cargs <- c("EWSR1.984.984.T.xlink_count.gtf.gz", "0.005022", "max");

plusnames <- NULL
mu_ref <- as.numeric(cargs[2]);
a <- getDataFromGtfPlus(cargs[1], plusnames);

a$len <- a$end - a$start + 1;

a$copies_count <- as.numeric(a$copies_count)
a$mut_count <- as.numeric(a$mut_count)

a$nomut_count  <- a$copies_count - a$mut_count
head(a)

#----
#a$proportion <- (a$mut_count+aa)/(a$copies_count+aa+bb); #changed on 2012.01.18 when the high mutation rates with RNAseq are observed
a$proportion <- (a$mut_count)/(a$copies_count);

head(a)
#set median mut_count and copies_count
#error rate in CLIP is usually in the order of 0.004 = 4/1000= 1/250
#set the hyperparameters accordingly when they are not set
#XXX 250 is hard coded here!!!
aa <- max(1,ceiling(250 * mu_ref)); #max is taken in case aa = 0
bb <- 250 - aa; 
#set to the Tuschl 4SU incorporation rate
#aa <- 1
#bb <- 40
#aa <- 0.5
#bb <- 0.5
#----
#set outfile tag
tag <- paste(gsub("*(_count.gtf.gz)", "", cargs[1]), "_mu_ref",as.character(format(mu_ref,nsmall=6,scientific=F)), sep='');
if (length(cargs) > 2) {
	tag <- paste(gsub("*(_count.gtf.gz)", "", cargs[1]), "_mu_ref",cargs[3], sep='');
}
#----

#calculate pbeta
a$pbeta <- pbeta(mu_ref, a$mut_count+aa, a$nomut_count+bb, log=T);

#remove NAs
b <- a[!is.na(a$pbeta),];
a <- b;

#write output to a file
options(scipen=10)
a$start <- a$start - 1; #for bed format
fields <- c("seqname", "start", "end", "gene_id", "pbeta", "strand");
write.table(a[sort(as.numeric(a$pbeta),decreasing=F,index.return=T)$ix,fields],quote=F, sep="\t", file = paste("xlinkEnrichment/", tag, ".bed", sep=""),row.names=F, col.names=F);

fields <- c("seqname", "start", "end", "gene_id", "mut_count", "strand");
write.table(a[sort(as.numeric(a$mut_count),decreasing=T,index.return=T)$ix,fields],quote=F, sep="\t", file = paste("xlinkCount/", tag, ".bed", sep=""),row.names=F, col.names=F);

fields <- c("seqname", "start", "end", "gene_id", "copies_count", "strand");
write.table(a[sort(as.numeric(a$copies_count),decreasing=T,index.return=T)$ix,fields],quote=F, sep="\t", file = paste("tagCount/", tag, ".bed", sep=""),row.names=F, col.names=F);
#----
library(latticeExtra)
library(hexbin)

tit <- paste(gsub("*(_count.gtf.gz)", "", cargs[1]), 
	"\ta:",as.character(format(aa,nsmall=0,scientific=F)),
	"\tb:",as.character(format(bb,nsmall=0,scientific=F)),
	"\tmu_ref:",as.character(format(mu_ref,nsmall=6,scientific=F)), 
	sep='');
plot.new(); pdf(paste(tag, "_params.pdf", sep=''), h=8, w=8);
hexbinplot( copies_count ~ proportion, data=a, 
	xlab = "T-to-C Proportion", ylab = "Tag Coverage",
	par.settings = poster.theme,
	scales=list(y=list(log=T)), 
	main=cargs[1],
#	style = "lattice",
	colorkey = T, aspect = 1, style = "nested.centroids", border = F,
	type = c("g"))
histogram(~ proportion, data=a, type="percent",
	xlab = "T-to-C Proportion",
	par.settings = poster.theme, col="darkgray",
 	main=cargs[1], nint=50)
levelplot( exp(pbeta) ~ proportion * copies_count, #| equal.count(copies_count,4,overlap=0), 
	data=a[as.integer(runif(500,1,dim(a)[1])),], 
	xlab = "T-to-C Proportion", ylab = "Tag Coverage",
	scales=list(y=list(log=T)),
	par.settings = poster.theme,
	main=paste(tit,"xlinkEnrichment Score",sep="\n"),
	cuts=9, 
	prepanel=prepanel.default.xyplot, 
	panel = function(...) {
		panel.grid(h=-1,v=-1)	
		panel.levelplot.points(...) 
	},
	col.regions=colorRampPalette(c("blue","white","red")), colorkey = list(space="top"))
histogram(~ exp(pbeta), data=a, type="percent",
 	xlab = "xlinkEnrichment Score",
	par.settings = poster.theme, col="darkgray",
 	main=tit, nint=50)
dev.off();
 
quit()
#----
library(plotrix)
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

