getDataFromGtfPlus <- function (filename, onamesf) {
	aa <- read.table(gzfile(filename), header=F);
	onames <- levels(unlist(read.table(cargs[2], header=F)));

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


cargs <- commandArgs(trailingOnly = T)
#cargs <- c("hg19.clusters.transcripts.genes.2.sites.gtfplus.gz", "hg19.clusters.transcripts.genes.2.sites.onames");

a <- getDataFromGtfPlus(cargs[1], cargs[2]);

pacodes <- unique(a$PAcode)
libs <- names(a)[grep("^lib", names(a))]

D <- NULL;
pacode <- pacodes[1]
for (pacode in pacodes) {
	d <- a[a$PAcode==pacode, names(a)[grep("^lib", names(a))] ];
	d
	x <- NULL
	x$mu <- mean(d)
	x$sd <- sd(d)
	x$N  <- dim(d)[1]
	x
	x <- data.frame(x)
	x$tag <- pacode
#	x$tag <- paste(sub("lib.","",row.names(x)),pacode,sep="_")
	x$lib <- sub("lib.","",row.names(x))
	row.names(x) <- NULL
	D <- rbind(D, x)
}

D$se <- D$sd / sqrt(D$N);

a <- D
a
library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''), h=8, w=8);
#segplot(reorder(factor(tag), mu) ~ (mu-se) + (mu+se), data = a,
segplot(factor(tag) ~ (mu-se) + (mu+se) | lib, data = a,
	layout=c(1,length(libs)),
	par.settings = poster.theme,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	cex=1.5,font=2,lwd=2,
	draw.bands = F, centers = mu, xlab="TPM",
	main = levels(a$type[1]), col.symbol=T, col="black",
	segments.fun = panel.arrows, ends = "both", angle=90, length=min(a$se), unit="mm",
	panel = function(...) {
		panel.grid(v = -1, h = 0)
		panel.segplot(...)
	})
dev.off()

quit()

library(plotrix)
xp <- barplot(a$mu,col="gray") #,ylim=c(0,max(a$mu+a$se)))
l <- length(xp)

plotCI(barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4)),a$mu,a$se,pch=21,xlim=c(0,l))
plot.new(); pdf(paste(cargs[1],'.pdf',sep=''), h=8, w=8);
plotCI(barplot(a$mu,col="gray",ylim=c(0,max(a$mu+a$se)*1.4)),a$mu,a$se,2*a$se,add=TRUE,pch=21,xlim=c(0,l))
text(xp,0, a$tag, srt = 90, adj=c(1,0.5))
text(xp,a$mu+5*max(a$se), a$tag, srt = 90, adj=c(0,0.5))
dev.off()

a

