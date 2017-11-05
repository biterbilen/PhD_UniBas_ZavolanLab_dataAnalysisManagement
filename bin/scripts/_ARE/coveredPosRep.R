library(lattice)

cargs <- commandArgs()
#cargs <- c(1,1,1,"../Analysis/Reproducibility/DT_IC.stats")

a <- read.table(cargs[4],header=T);
head(a)

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)      ## set as default

#covered positional reproducibility
tit <- paste("Features Covered by MTC");
plot.new(); pdf(paste(cargs[4], '.coveredposrep','.pdf',sep=''), width=8, height=8)
barchart( 100 * coveredPositions / totalPositions ~ protein| feature, data=a, origin = 0,
	auto.key = list(title = "Feature", space="bottom",cex=1.5,font=2),
	par.strip.text = list(cex=1.5,font=2), type=c("g"),
	par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1.5,
	main=tit,
	ylab="Proportion (%)", scales=list(x=list(rot=90)))
#barchart(prop.table(b, margin = 1), xlab = "Proportion", auto.key = list(adj = 1)) 
dev.off();

