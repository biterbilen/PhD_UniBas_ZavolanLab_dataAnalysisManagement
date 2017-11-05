cargs <- commandArgs(trailingOnly = T);

#cargs <- c("2.average_clip_posterior_IG_perGeneSymbol.tab","2");

a <- read.table(cargs[1],header=F);
names(a) <- c("x","gid", "x", "tid", "x", "gname", "x","PACode", "Xlink.Score", "PAshift.Score");
head(a)

b <- a[a$PACode == paste(cargs[2],1,sep="."), c("gid","PAshift.Score","Xlink.Score")]
b[c(paste("Xlink.Score.PAsite",cargs[2], 1, sep="."))] <- b$Xlink.Score;
b$Xlink.Score <- NULL;
head(b)
for (i in 2:as.numeric(cargs[2])) {
	a2 <- a[a$PACode == paste(cargs[2],i,sep="."), c("gid","Xlink.Score")]
	a2[c(paste("Xlink.Score.PAsite", cargs[2], i, sep="."))] <- a2$Xlink.Score;
	a2$Xlink.Score <- NULL;
	tmp <- merge(b, a2, by="gid")
	b <- tmp
	head(tmp)
	head(a2)
}
head(b)

library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

d <- b[b$PAshift.Score>as.numeric(cargs[3]),];
b <- d;

n <- 64;
gcolors <- gray.colors(n, start = min(b$PAshift.Score), end = 1, gamma = 2.2)
breaks <- do.breaks(range(b$PAshift.Score), n)
b$color <- level.colors(b$PAshift.Score, at = breaks, col.regions = gcolors)

plot.new(); pdf(paste(cargs[1], '.PAs_vs_IG.pdf', sep=''), h=10, w=10);
splom(~b[,grep("Xlink.Score.PAsite",names(b))],
	groups=b$color,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	cex=1.5,font=2,
	varname.cex=1,varname.font=2,
	xlab = "",
#	ylab = "Xlink Score of PolyA Site 2.2",
 	type=c("g","p"),
	main="Score for Change of PolyA Site Usage",
	legend = list(top = list(fun = draw.colorkey, args = list(key = list(space="top", col = gcolors, at = breaks), draw = FALSE))),
	panel = function(x, y, groups, ..., subscripts) { 
		fill <- groups[subscripts] 
		panel.grid(h = -1, v = -1) 
		panel.xyplot(x, y, pch = 21, fill = fill, col="black", ...) 
	},
	lower.panel = function(x, y, ...){
#		panel.hexbinplot(x, y, ...)#,
		panel.grid(h = -1, v = -1) 
		panel.abline(a=0,b=1, ..., col = 'darkgray',lwd=2)#,
		panel.loess(x, y, ..., col = 'black',lwd=2)#,
		panel.key(c(sprintf('N = %d\nR = %.2f', length(x),cor(x,y) )),col='black',cex=1.5, font=2,points=F)
	}, 
	)

dev.off()

