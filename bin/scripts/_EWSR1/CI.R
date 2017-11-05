cargs <- c("affinity.stats")
cargs <- commandArgs(trailingOnly = T)

a <- read.table(cargs[1], header=T)
head(a,2)


xlabel <- "mu";
ylabel <- "";
if (length(cargs)>1) {
	xlabel <- cargs[2];
} 
titel <- levels(a$type[1]);
if (length(cargs)>2) {
	titel <- cargs[3];
} 
if (length(cargs)>3) {
	ylabel<- cargs[4];
} 
linescount <- 1
if (length(cargs)>4 & as.numeric(cargs[5]) == 1) {
	library(plyr)
	Ns <- ddply(a, .(type), function(x) { sum(x$N); })
	b <- merge(a, Ns);		
	b$type <- paste(b$type,"\nN=",b$V1,sep="")
	a  <- b
	linescount <- 2
}
rel <- "same";
if (length(cargs)>5 & (cargs[6] == "same" | cargs[6] == "free")) {
	rel <- cargs[6];
}

print(titel)
a$se <- a$sd / sqrt(a$N);

source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))
library(latticeExtra)
print(unique(a$type))
cols <- 1;
rows <- length(unique(a$type));
if ( rows > 6) {
	cols <- 3;
} else if (rows > 3) {
	cols <- 2;
}
rows <- ceiling(rows/cols);
if ( rows > 3 ) {
	rows <- 3; 
}

#TODO
#hor <- TRUE
#if (hor == T) {
#	tmp <- xlabel;
#	xlabel <- ylabel;
#	ylabel <- tmp;
#}

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''), h=8, w=8);
#segplot(reorder(factor(tag), mu) ~ (mu-se) + (mu+se) | type, data = a,
segplot(factor(tag) ~ (mu-se) + (mu+se) | type, data = a,
#	horizontal = hor,
	scales=list(x=list(relation=rel)),
	main = titel,
	layout=c(cols,rows),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2,lines=linescount),
	par.settings = poster.theme,
	draw.bands = F, centers = mu, xlab=xlabel, ylab=ylabel,
	col.symbol=T,col=T,
	lwd=2, cex=1.5, font=2, pch="x",
	segments.fun = panel.arrows, ends = "both", angle=90, length=min(1,max(a$sd)), unit="mm",
	panel = function(...) {
#		panel.text(..., labels=a$N, adj=1, srt=90, pos=4)
		panel.grid(v = -1, h = -1)
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
#text(xp,a$mu+5*max(a$se), a$tag, srt = 90, adj=c(0,0.5))
dev.off()



