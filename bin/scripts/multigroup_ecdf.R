library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("a", "all", "Likelihood of Differential Expression", "Expressed in HeLa and HEK293", "3", "0.01","less", "F");
cargs <- c("motifcount", "neg", "Motif Count", "Expressed in HeLa and HEK293", "3", "0.01","less", "F");
cargs <- c("1.$region.displacement","INTRON","x","main","6","0.000001","greater")
cargs <- c("288.closest_distances","299","x","main","3","0.01","greater")
cargs <- commandArgs(trailingOnly=T)
#cargs <- c("288.closest_distances", "238_239")
#cargs <- c("288")

id2 <- "pnorm"
xlabel <- "Closest Distance"
maintext <- "Distribution of Closest Distances to Crosslink Sites"
valIndex <- 3
pvalCut <- 0.05
alternative <- "two.sided"
clog <- T
oneSamptest <- "punif"; #pnorm etc. #test
ftag <- ""
plotref <- 1

if (length(cargs)>1) id2      <- cargs[2]
if (length(cargs)>2) xlabel   <- cargs[3]
if (length(cargs)>3) maintext <- cargs[4]
if (length(cargs)>4) valIndex <- as.numeric(cargs[5])
if (length(cargs)>5) pvalCut  <- as.numeric(cargs[6])
if (length(cargs)>6) alternative <- cargs[7]
if (length(cargs)>7) clog     <- as.logical(cargs[8])
if (length(cargs)>8) oneSamptest <- cargs[9]
if (length(cargs)>9) ftag     <- cargs[10]
if (length(cargs)>10) plotref <- cargs[11]

a <- read.table(cargs[1],header=T, stringsAsFactors=FALSE,sep="\t")
names(a)[valIndex] <- "value"; #FIXME
head(a,2)
dim(a)

#print (c(alternative,id2,xlabel,maintext,valIndex,pvalCut,clog))
groups <- unique(a$group)
if (id2 == "NULL") {
	x2 <- oneSamptest;
	id2 <- x2;
} else {
	x2 <- a[a$group==id2,"value"]
}

i <- "DGE"
b <- NULL
barinput <- NULL
barinput$name <- groups
barinput$count <- 0
barinput$plot <- T
barinput <- data.frame(barinput)

for (i in groups) {
	print (i)
	n <- which(a$group==i);
	l <- length(n)
	id <- paste(a$group[n], "\t(n=",l,")",sep="");
	x1 <- a[a$group==i,"value"]
	ks <- ks.test(x1,x2,alternative=alternative)
	#ks <- ks.test(x1,x2,alternative=alternative,exact=F)
	pv <- sprintf("%.7g",ks$p.value)
	st <- sprintf("%.2f",ks$statistic)
	write.table(format(digits=2, data.frame(list(method="KS.test", alternative=alternative,labelfore=i,labelback=id2,lengthfore=length(x1),lengthback=length(x2),pvalue=pv,statistic=st))),sep="\t",row.names=F, quote=F);
	if (ks$p.value < pvalCut | i==id2) {
		id <- paste(id, " (alternative=",alternative, ") (p-value=", pv, ") (statistic=",st,") (compared2=", id2, ")",sep="");
		a$group[n] <- id
		#b <- rbind(b,a[n,])
		if (length(n)>5000) { #for rendering more efficiently
			b <- rbind(b,a[sample(n,5000),])
		} else {
			b <- rbind(b,a[n,])
		}
		barinput[barinput$name == i,]$count = l
	} else {
		barinput[barinput$name == i,]$plot <- F
	}
}

groups <- unique(b$group)
clrs <- article.colors;
if (length(grep("median",groups))>0) { #PATCH for median graphs
	clrs <- c(article.colors[3],article.colors[3],article.colors[1],article.colors[2])
}
print(groups)
#quit()
if (clog == T) {
	b[b$value==0,"value" ] <- 1
}
if (length(groups) > 0) {
	p <- ecdfplot(~ value|factor(panel), groups=group, data=b, id2=id2,
		main=maintext,
		xlab = xlabel,ylab="Cumulative density",
		#par.settings = article.theme, lwd=article.lwd, cex=article.cex, origin=0,col=clrs,
		par.settings = article.theme, col=clrs,#origin=0,
		strip = strip.custom(style=1,bg="white"),
		par.strip.text = list(cex=article.cex,font=article.font),
		xscale.components = xscale.components.log10ticks,
		auto.key=list(space="top",columns=1, cex=article.cex,font=article.font, fontfamily=article.fontfamily,lwd=article.lwd, col=clrs,lines=F),
		panel = function(x,y,...,groups=groups,subscripts=subscripts,data2=data2,id2=id2) {
			panel.grid(h = -1, v = -1)
			#TODO not
			if (plotref == 1 && id2 == "pnorm") {
				panel.ecdfplot(rnorm(100,0.5,0.5/3),lwd=article.lwd,lty="dashed",...)
			}	else if (plotref == 1 && id2 == "punif") {
				panel.abline(a=0,b=1,lwd=article.lwd,lty="dashed",...)
			}
			panel.ecdfplot(x,groups=groups,subscripts=subscripts, ...)
		},
		scales=list(x=list(log=clog)))

	p2 <- barchart(count ~ reorder(factor(name), -count), data=barinput, drop.unused.levels=T,
		scales=list(x=list(rot=0)), stack=T,
		ylab="Count",xlab="Group",subset=plot==T, group=name,
		par.settings = article.theme, origin=0,
		strip = strip.custom(style=1,bg="white"),par.strip.text = list(cex=article.cex,font=article.font),
		auto.key = list(space="top",points=F, rectangles=F, columns=length(barinput$plot==T), col=clrs,columns=1, cex=article.cex,font=article.font, fontfamily=article.fontfamily),
		prepanel = function(x, y,...) {
			list(xlim=levels(reorder(x, -y)))
		},
		panel = function(x,y,...) {
			panel.grid(h = -1, v = -1)
			panel.barchart(x, y, ...)
		}
		)

	plot.new(); pdf(paste(cargs[1],ftag,'.pdf',sep=''),h=8,w=8)
	print(p)
	print(p2)
	dev.off()
} else{
	print("no groups are significantly different")
}


