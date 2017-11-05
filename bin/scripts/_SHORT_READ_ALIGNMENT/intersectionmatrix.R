cargs <- commandArgs();

#cargs <- c(1,1,1,"trusted.matrix", "scaled", 500)

a <- read.table(cargs[4], header=T);

if (cargs[5] == "scaled") {
	x  <- t(as.matrix(scale(a)))
} else if (cargs[5] == "unscaled") {
	x  <- t(as.matrix((a)))
}
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)

library(latticeExtra)
ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

plot.new(); pdf(paste(cargs[4], ".", cargs[5], '.pdf', sep=''),h=8,w=8);
levelplot(x[row.ord, col.ord],
	par.settings=list(axis.text=list(font=2,cex=0.7),
		par.strip.text=list(font=2,cex=1.0),
		par.ylab.text=list(font=2,cex=1.0),
		par.xlab.text=list(font=2,cex=1.0),
		par.main.text=list(font=2,cex=1.0)),
	xlab="", ylab="",
	aspect = "fill", main=paste(cargs[4], cargs[5], 'top', cargs[6]),
	scales = list(x = list(rot = 90)),
	colorkey = list(space = "left"),
	legend =
	list(right =
		list(fun = dendrogramGrob,
			args =
			list(x = dd.col, ord = col.ord,
				side = "right",
				size = 10)),
		top =
		list(fun = dendrogramGrob,
			args =
			list(x = dd.row, 
				side = "top",
				type = "triangle"))))
dev.off();

