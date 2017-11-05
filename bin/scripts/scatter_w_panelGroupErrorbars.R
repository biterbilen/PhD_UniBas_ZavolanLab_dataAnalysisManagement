library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("INTRON.displacement4plot","4","7","5","6","1")
cargs <- c("INTRON.displacement4plot","4","7","5","6","1")
cargs <- c("a","4","7","5","6","1")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=T,sep="\t");
indices <- as.numeric(cargs[2:5])
minvalue <- as.numeric(cargs[6])
l <- dim(a)[2]
if (indices[3] == -1) {
	indices[3] <- l + 1
	a[indices[3]] <- 0
}
if (indices[4] == -1) {
	indices[4] <- l + 2
	a[indices[4]] <- 0
}
nms <- names(a)[indices]
names(a)[indices] <- c("x","mu","sd","N");
a$se <- a$sd / sqrt(a$N)
head(a,1)

panel.ci <- function(x,y,ly,uy,subscripts,col,pch,lwd,lty,cex,fontfamily,font,...) {
	ly <- ly[subscripts]
	uy <- uy[subscripts]
	panel.grid(h = -1, v = -1)
	panel.arrows(x, ly, x, uy, length=0.0, angle=0,col="black",alpha=0.3)
	panel.xyplot(x,y, pch=pch,...)
	panel.abline(lm(y ~ x),lwd=lwd,...)
	panel.abline(h=median(y),v=median(x),lty="dashed",lwd=lwd,...)
	panel.abline(a=0,b=1,lty="dashed",lwd=lwd,col="red",...)
	#panel.loess(x,y,...)#,col=article.colors[3],lwd=article.lwd*2)
	library(grid)
	#grid.text(max(x),max(y),sprintf("N=%.0f\nR=%.2f",length(x),cor(x,y,method="spearman")),cex=cex,fontfamily=fontfamily,font=font)
	grid.text(sprintf("N=%.0f\nR=%.2f",length(x),cor(x,y,method="spearman")),x=0.85,y=0.85,just="centre")
}

prepanel.ci <- function(x, y, ly, uy, subscripts, ...) {
	x <- as.numeric(x)
	ly <- as.numeric(ly[subscripts])
	uy <- as.numeric(uy[subscripts])
	list(ylim = range(y, uy, ly, finite = TRUE))
}

lr <- 1
lc <- length(unique(a$group))

if (lc > 3) {
	lr <- 2
	lc <- ceiling(lc/lr)
}

plot.new(); pdf(paste(cargs[1],"_scatter.pdf",sep=''),h=8,w=8)
xyplot(mu ~ x|factor(panel), data=a, subset=x>minvalue,
	ly=a$mu-a$se,uy=a$mu+a$se,
	xlab=nms[1],ylab=nms[2],
	strip = strip.custom(style=1,bg="white"),
	group=group,layout=c(lr,lc),
	par.settings=article.theme,pch="o",
	auto.key=list(column=length(unique(a$group)),lines=F,points = F,rectangles=F,col=article.colors,font=article.font,fontfamily=article.fontfamily),
	xscale.components = xscale.components.log10ticks, 
	yscale.components = yscale.components.log10ticks, 
	scale=list(x=list(log=T),y=list(log=F)),main=cargs[1], 
	prepanel = prepanel.ci,
	panel = panel.superpose,
	panel.groups = panel.ci,
	)
dev.off()


