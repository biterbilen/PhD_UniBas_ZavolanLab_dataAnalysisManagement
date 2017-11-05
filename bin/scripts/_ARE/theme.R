article.colors <- c("black","slategray","darkred","blue3","gold4","darkorchid4");
article.alpha <- 0.9
article.pch <- "."
article.font <- 1
article.cex <- 1
article.lwd <- 2
article.fontfamily <- "sans"

article.theme <- list(
#	box.umbrella    = list(lty=1, col="black", pch="."),
	#box.rectangle = list(col="black",alpha=article.alpha,lwd=6),
#	box.rectangle = list(col="red",lwd=0, fill="black", pch="."),
	# strip.border  = list(col="red"),
	superpose.symbol  = list(col=article.colors),
	superpose.line    = list(col=article.colors),
	superpose.polygon = list(col=article.colors),
	axis.line     = list(col=NA), #NA
	axis.components=list(left=list(pad1=0.2),bottom=list(pad1=1)),
	# layout.heights=list(strip=1,main=0,main.key.padding=1,key.axis.padding=1,axis.xlab.padding=1,xlab.key.padding=1,bottom.padding=0,key.bottom=0,top.padding=0,main.key.padding=0,xlab=0,axis.panel=0,panel=1),
	# layout.widths=list(strip=1,main=0,main.key.padding=0,key.axis.padding=0,axis.xlab.padding=0,xlab.key.padding=0,bottom.padding=0,key.bottom=0,top.padding=0,main.key.padding=0,xlab=0,axis.panel=0,panel=1),
	# regions       = list(col="red"),
	# clip = list(panel=5),
	grid.pars     = list(fontfamily=article.fontfamily, lwd=article.lwd), #serif mono sans
	varname       = list(font=article.font, cex=article.cex),
	axis.text     = list(font=article.font, cex=article.cex, alpha=0.7, lineheight=0.1),
	strip.text    = list(font=article.font, cex=article.cex),
	par.main.text = list(font=article.font, cex=article.cex),
	par.ylab.text = list(font=article.font, cex=article.cex),
	par.xlab.text = list(font=article.font, cex=article.cex)
	)


#ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
#ltheme$strip.background$col <- "transparent" ## change strip bg
#lattice.options(default.theme = ltheme, lwd=4)      ## set as default
fnt <- 2;
poster.theme <- list(
	superpose.line = list(lwd=2),
	varname       = list(font=fnt,cex=1.5),
	axis.text     = list(font=fnt,cex=1.5),
	strip.text    = list(font=fnt,cex=1.5),
	par.main.text = list(font=fnt,cex=1.5),
	par.ylab.text = list(font=fnt,cex=1.5),
	par.xlab.text = list(font=fnt,cex=1.5)
	) 

#my.theme <- theEconomist.theme() 

