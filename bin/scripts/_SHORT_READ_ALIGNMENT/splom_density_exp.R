cargs <- commandArgs(trailingOnly=T);

source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))
library(latticeExtra);
library(hexbin);

a <- read.table(cargs[1], header=T, row.names="id");
rawValue <- 0;
if (length(cargs) > 2) {
	rawValue <- as.numeric(cargs[3]);	
}
lab <- "Expression";
if (length(cargs) > 3) {
	lab <- cargs[4];
}

if (rawValue == 0) {
	lab <- paste (lab, "(log2)");
}

ofiletag <- cargs[1];
getpc <- function(x,cs) { 
	colmin <- 1:cs * 0;  
	for (c in 1:cs) {
		colmin[c] <- sort(x[x[,c]>0,c])[1] #get min that is greater than 0
	}; 
	return(colmin);
};

rs <- dim(a)[1]
cs <- dim(a)[2] 
mins <- getpc(a,cs);
d <- a + matrix(rep(mins,rs),nrow=rs,byrow=T)

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

plot.new(); pdf(paste(ofiletag,'.splom.pdf', sep=''), h=10, w=10);
dd <- d;
if (rawValue == 0) {
	dd <- log2(d);
}
splom(dd, main=ofiletag, type=c("g","p"),
	varname.cex=.8, axis.text.cex=.8, 
	xlab=lab,ylab=lab,
	par.settings = poster.theme,
	upper.panel = function(x, y, ...){
		panel.hexbinplot(x, y, ...)#,
		panel.abline(a=0,b=1, ..., col = 'darkgray',lwd=2)#,
		panel.loess(x, y, ..., col = 'black',lwd=2)#,
		panel.key(c(sprintf('N=%d\nR=%.2f', length(x), cor(x,y, use = "pairwise.complete.obs", method="spearman") )),...,col='black', font=2,cex=1.5, points=F)
	},
	lower.panel = function(x,y, ...) {
#	panel.key(c(sprintf('N=%d\n%.2f', length(x), cor(x,y, use = "pairwise.complete.obs", method="spearman") )),...,col='black', font=2,cex=1.5, corner=c(0.5,0.5), points=F)
	}
	)
dev.off();

#---------------------------
cor.a <- cor(dd, use = "pairwise.complete.obs", method="spearman");
ord <- order.dendrogram(as.dendrogram(hclust(dist(cor.a))))
x <- cor.a
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)
dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)

library(latticeExtra)

plot.new(); pdf(paste(cargs[1],".heatmap.pdf",sep=""), h=8,w=8);
tit <- "Spearman Correlation";
levelplot(x[row.ord, col.ord],
	par.settings = poster.theme,
	aspect = "fill",
	xlab="",ylab="",
	main=list(label=tit),
	col.regions=colorRampPalette(c("blue","white","red")),
	scales = list(x = list(rot = 90)),
	colorkey = list(space = "top"),
	legend =
	list(right =
		list(fun = dendrogramGrob,
			args =
			list(x = dd.col, ord = col.ord,
				side = "right",
				size = 10))))
dev.off()

#---------------------------


a <- read.table(cargs[2], header=T);
if (rawValue == 0) {
	a$exp <- log2(a$exp);
}
plot.new(); pdf(paste(ofiletag,'.density.pdf', sep=''), h=8, w=8)
densityplot(exp_type ~ exp, data=a, group=exp_type, plot.points=F, 
	xlab=lab,
	par.settings = poster.theme,
	na.rm=T, auto.key=list(space="top"), main=basename(ofiletag)) +
layer(panel.refline(v=0,lwd=2,cex=1))
dev.off();

plot.new(); pdf(paste(ofiletag,'.bw.pdf', sep=''), h=8, w=8);
bwplot(exp_type ~ exp, data=a, na.rm=T, 
	xlab=lab, col="black",
	par.settings = poster.theme,
	main=basename(ofiletag)) +
layer(panel.refline(v=0,lwd=2,cex=1))
dev.off();


