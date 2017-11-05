gcd <- function(a,b) ifelse (b==0, a, gcd(b, a %% b)) 

log_sum_exp <- function (x, y) {
	if (x > y) {
		return(x+log(1+exp(y-x)));
	} else {
		return(y+log(1+exp(x-y)));
	}
}

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
	# names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)));
	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)), onames)
	return(a)
																	
}

cargs <- commandArgs()
#cargs <- c(1,1,1,"deneme.mu-med-stdev-count-freqdesc", "dummy", "mean", "median", "stdev", "count", "freqdesc");
#cargs <- c(1,1,1,"/import/bc2/home/zavolan/bilebi00/_ARE/Analysis/xlinkEnrichment/basedOnNucleotide/271.clusters.copiesT_count.gtf.gz", "/import/bc2/home/zavolan/bilebi00/_ARE/Analysis/xlinkEnrichment/basedOnNucleotide/271.clusters.xlinkT_count.gtf.gz", "0.1", "1")

mutrate <- -1
wplot <- 0
plusnames <- NULL
if (length(cargs) > 5) {
	mutrate <- as.numeric(cargs[6]);
}

if (length(cargs) > 6) {
	wplot <- as.numeric(cargs[7]);
}

if (length(cargs) > 7) {
	plusnames <- cargs[8:length(cargs)];
}

x <- getDataFromGtfPlus(cargs[4], plusnames);
y <- getDataFromGtfPlus(cargs[5], plusnames);
dim(x)
dim(y)
a <- merge(x,y);
dim(a)
#a <- merge(x,y,by=c("seqname", "start", "end", "strand"))
a$nomut_count  <- a$copies_count - a$mut_count
a$len <- a$end - a$start + 1;
#----
#a$freq <- (a$mut_count+aa)/(a$copies_count+aa+bb); #changed on 2012.01.18 when the high mutation rates with RNAseq are observed
a$freq <- (a$mut_count)/(a$copies_count);
back <- a
a <- a[a$copies_count>0, grep("y$", names(a), perl=T, value=T,invert=T)];

#set the hyperparameters accordingly when they are not set
#error rate is usually in the order of 0.004 = 4/1000 = 1/250
if (mutrate < 0) {
	mutrate <- mean(a$freq);
}

#pbeta hyperparams default values
aa <- integer(250 * mutrate);
bb <- 250 - aa; 
print(c(mutrate, aa, bb, dim(a), mean(a$freq)))
#----
tag <- paste(gsub("*(_count.gtf.gz)", "", cargs[5]), "_mutrate",as.character(format(mutrate,nsmall=6,scientific=F)), sep='');
#----
k <- a$freq
tag1 <- gsub("*(_count.gtf.gz)", "", cargs[5]);
write.table(list(tag=tag1, mu=mean(k), sd=sd(k), N=length(k)),sep="\t",row.names=F, quote=F, file = paste(tag, ".freqstat", sep=""))
if (wplot == 0) {
	quit();
}
#----
head(a)
a$pbeta <- pbeta(mutrate, a$mut_count+aa, a$nomut_count+bb, log=F);
head(a)
fields <- c("seqname", "start", "end", "gene_id", "pbeta", "strand")
write.table(a[sort(a$pbeta,decreasing=F,index.return=T)$ix,fields],quote=F, sep="\t", file = paste(tag, ".bed", sep=""),row.names=F, col.names=F);

#----
library(lattice)
library(hexbin)
library(plotrix)

b <- rbind(data.frame(list(type="pbeta", value=a$pbeta)), data.frame(list(type="freq", value=a$freq))) 
p1 <- densityplot(~ value, groups=type, data=b,plot.points=F, type="percent", auto.key=list(columns=1, corner=c(1,1)), main=tag, ylab = "percent")
p2 <- xyplot(pbeta ~ freq, data=a, panel=panel.hexbinplot, main=tag);
p3 <- xyplot(pbeta ~ log2(len), data=a, panel=panel.hexbinplot, subset=copies_count>0, main=tag);
p4 <- xyplot(freq ~ log2(len), data=a, panel=panel.hexbinplot, main=tag);

plot.new(); pdf(paste(tag, "_params.pdf", sep=''), h=8, w=8);
plot(p1, split = c(1, 1, 2, 2))
plot(p2, split = c(1, 2, 2, 2), newpage = FALSE)
plot(p3, split = c(2, 1, 2, 2), newpage = FALSE)
plot(p4, split = c(2, 2, 2, 2), newpage = FALSE)
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

