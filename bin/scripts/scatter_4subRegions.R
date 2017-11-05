library(latticeExtra)
library(reshape)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("a", "288", "289") 

a <- read.table(cargs[1],header=T)

head(a)

pos <- 1
p <- xyplot((INTRON) ~ (INTRON.1), data=a, subset=INTRON>0 & INTRON.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "INTRON", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p2 <- xyplot((X3UTR) ~ (X3UTR.1), data=a, subset=X3UTR>0 & X3UTR.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "3'UTR", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p3 <- xyplot((CDS) ~ (CDS.1), data=a, subset=CDS>0 & CDS.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "CDS", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p4 <- xyplot((X5UTR) ~ (X5UTR.1), data=a, subset=X5UTR>0 & X5UTR.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "5'UTR", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p5 <- xyplot((UGS) ~ (UGS.1), data=a, subset=UGS>0 & UGS.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "UGS", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p6 <- xyplot((DGS) ~ (DGS.1), data=a, subset=DGS>0 & DGS.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "DGS", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p7 <- xyplot((UGE) ~ (UGE.1), data=a, subset=UGE>0 & UGE.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "UGE", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p8 <- xyplot((DGE) ~ (DGE.1), data=a, subset=DGE>0 & DGE.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "DGE", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p9 <- xyplot((X3UTR+UGE) ~ (X3UTR.1+UGE.1), data=a, subset=X3UTR+UGE>0 & X3UTR.1+UGE.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "3'UTR+UGE", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p10 <- xyplot((X5UTR+DGS) ~ (X5UTR.1+DGS.1), data=a, subset=X5UTR+DGS>0 & X5UTR.1+DGS.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "5'UTR+DGS", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p11 <- xyplot((X3UTR+UGE+DGE) ~ (X3UTR.1+UGE.1+DGE.1), data=a, subset=X3UTR+UGE+DGE>0 & X3UTR.1+UGE.1+DGE.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "3'UTR+UGE+DGE", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

p12 <- xyplot((X5UTR+DGS+UGS) ~ (X5UTR.1+DGS.1+UGS.1), data=a, subset=X5UTR+DGS+UGS>0 & X5UTR.1+DGS.1+UGS.1>0,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main="Crosslinked Positions in Common Genes", 
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.abline(a=0,b=1,col="gray",lwd=article.lwd)
		panel.text(mean(x),max(y), sprintf("%s\nN=%.0f\nR=%.2f", "5'UTR+DGS+UGS", length(x),cor(x,y,method="pearson")),cex=0.8, col=article.colors[1],fontfamily=article.fontfamily,pos=pos, lwd=article.lwd)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=T)))

plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=8,w=8)
#update(c(p4,p5,p6,p,p2,p3))
update(c(p,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12))
dev.off()

