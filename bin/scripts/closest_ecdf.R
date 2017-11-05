library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("288.closest_distances", "238_239")
#cargs <- c("288")

id2 <- NULL
ylabel <-	"Empirical cumulative distribution"
xlabel <- "Closest distance"
maintext <- "Distribution of Closest Distances to Crosslink Sites"
clog <- T
type <- "DUMMY"

if (length(cargs)>1) id2      <- cargs[2]
if (length(cargs)>2) xlabel   <- cargs[3]
if (length(cargs)>3) maintext <- cargs[4]
if (length(cargs)>4) clog     <- as.logical(cargs[5])
if (length(cargs)>5) type     <- cargs[6]

a <- read.table(cargs[1],header=T, stringsAsFactors=FALSE)
names(a)[3] <- "distance" #for backward compatibility

maxs <- aggregate(distance ~ name + tag, data=a, max)
names(maxs) <- gsub("distance", "maxdistance", names(maxs))
b <- merge(a,maxs,by=c("tag","name"))

#correction for -1 and 0 distances
s <- b$distance == -1
b[s,"distance" ] <- b[s,"maxdistance"]
b[ b$distance == 0 ,"distance" ] <- 1

head(b,2)

tags <- unique(b$tag)

plts <- NULL

if (is.null(id2) || id2 == "null" || id2 == "NULL") {
	print ("if")
	nms <- unique(b$tag)
	for (i in nms) {
		n <- which(b$tag==i);
		l <- length(n)
		id <- paste(b$tag[n], "\t(n=",l,")",sep="");
		if ((!is.null(id2)) && i == id2)  { id2 <- id[1]; }
		b$tag[n] <- id
	}
	if (length(grep("type",names(b))) == 0) { b$type <- type }
		print(head(b))
p <- ecdfplot(~ distance|type, groups=tag, data=b,
	main=maintext,
	xlab = xlabel,ylab = ylabel,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex, origin=0,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	xscale.components = xscale.components.log10ticks,
#	auto.key=list(space="top",columns=length(unique(b$name)),cex=article.cex,font=article.font, fontfamily=article.fontfamily,lwd=article.lwd, col=article.colors,lines=F),
	auto.key=list(space="top",columns=1, cex=article.cex,font=article.font, fontfamily=article.fontfamily,lwd=article.lwd, col=article.colors,lines=F),
	panel = function(x,y,...,groups=groups,subscripts=subscripts,data2=data2) {
		panel.grid(h = -1, v = -1)
		panel.ecdfplot(x,groups=groups,subscripts=subscripts, ...)
	},
	scales=list(x=list(log=clog)))

} else {
	print ("else")
	nms <- unique(b$name)
	for (i in nms) {
		n <- which(b$name==i);
		l <- length(n)
		id <- paste(b$name[n], "\t(n=",l,")",sep="");
		if ((!is.null(id2)) && i == id2)  { id2 <- id[1]; }
		b$name[n] <- id
	}
p <- ecdfplot(~ distance | factor(tag), groups=name, data=b, data2=id2,
	layout=c(1,length(unique(b$tag))), 
	main=maintext,
	xlab = xlabel, ylab=ylabel,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex, origin=0,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	xscale.components = xscale.components.log10ticks,
#	auto.key=list(space="top",columns=length(unique(b$name)),cex=article.cex,font=article.font, fontfamily=article.fontfamily,lwd=article.lwd, col=article.colors,lines=F),
	auto.key=list(space="top",columns=1,cex=article.cex,font=article.font, fontfamily=article.fontfamily,lwd=article.lwd, col=article.colors,lines=F),
	panel = function(x,y,...,groups=groups,subscripts=subscripts,data2=data2) {
		panel.grid(h = -1, v = -1)
#		panel.abline(h=0.5, ..., col=article.colors[3],lty="dashed", ...)
		panel.ecdfplot(x,groups=groups,subscripts=subscripts, ...)
		if (!is.null(data2)) {
			g1 <- groups!=data2
			g2 <- groups==data2
			x1 <- x[subscripts[g1]]
			x2 <- x[subscripts[g2]]
			ks <- ks.test(x1,x2,alternatve="l",exact=F)
			print(length(x1))
			print(length(x2))
			panel.text(1,1,sprintf("p-value=%f\nstatistic=%f\n", ks$p.value, ks$statistic),pos=1)
		}
	},
	scales=list(x=list(log=clog)))
}
plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
print(p)
dev.off()


