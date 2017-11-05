library(latticeExtra)

source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))
cargs <- commandArgs(trailingOnly=T);
#cargs <- c("bothTails/saturation")
a <- read.table(cargs[1],header=T,sep="\t")
head(a)

#p1 <- xyplot( (subGeneCount/fullGeneCount) ~ fraction |proteinName, data=a, 
p1 <- xyplot( subGeneCount ~ fraction |proteinName, data=a, 
	layout=c(1,length(unique(a$proteinName))),
	scales=list(y=list(relation="free")),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme, col=article.colors, pch=19,
	groups=sampleid, lwd=article.lwd, cex=article.cex,
	auto.key=list(space="top",points=F,columns=3,col=article.colors),
	panel =	panel.superpose, 
	xlab="Sampling Fraction", 
	ylab="Recovered Genes", main="Targeted Gene Saturation",
	#ylab="Fraction of Genes Recovered", main="Targeted Gene Saturation",
	panel.groups = function(x, y, col, ...){
		panel.xyplot(x,y, ...)
		panel.grid(h = -1, v = -1)
#		panel.ablineq(lm(y~log(x)),r.sq=T,digits=5,varNames=alist(y=y, x="log(x)"),rot=T,lwd=article.lwd,cex=1,fontfamily=article.fontfamily,at=0.5,pos=3,col.text=col,col=col)
		yy <- lm(y ~ log(x))
		a <- as.numeric(yy$coefficients[1])
		b <- as.numeric(yy$coefficients[2])
		Rs <- summary(yy)$r.squared
		print(c(Rs,a,b))
		panel.curve(a+b*log(x), type="l", lwd=article.lwd,col=col)
		panel.text(max(x),max(y),sprintf("y=%.4f+%.4f*log(x)), R^2=%.2f",a,b,Rs),col=col,col.text=col,fontfamily=article.fontfamily,cex=0.8,pos=2,adj=c(2,2)) 
	})



p2 <- xyplot( subCount ~ fraction |proteinName, data=a, 
	layout=c(1,length(unique(a$proteinName))),
	scales=list(y=list(relation="free")),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme, col=article.colors, pch=19,
	groups=sampleid, lwd=article.lwd, cex=article.cex,
	auto.key=list(space="top",points=F,columns=3,col=article.colors),
	panel =	panel.superpose, 
	xlab="Sampling Fraction", 
	ylab="Recovered Positions", main="Targeted Position Saturation",
	#ylab="Fraction of Positions Recovered", main="Targeted Position Saturation",
	panel.groups = function(x, y, col, ...){
		panel.xyplot(x,y, ...)
		panel.grid(h = -1, v = -1)
#		panel.ablineq(lm(y~log(x)),r.sq=T,digits=5,varNames=alist(y=y, x="log(x)"),rot=T,lwd=article.lwd,cex=1,fontfamily=article.fontfamily,at=0.5,pos=3,col.text=col,col=col)
		yy <- lm(y ~ log(x))
		a <- as.numeric(yy$coefficients[1])
		b <- as.numeric(yy$coefficients[2])
		Rs <- summary(yy)$r.squared
		print(c(Rs,a,b))
		panel.curve(a+b*log(x), type="l", lwd=article.lwd,col=col)
		panel.text(max(x),max(y),sprintf("y=%.4f+%.4f*log(x)), R^2=%.2f",a,b,Rs),col=col,col.text=col,fontfamily=article.fontfamily,cex=0.8,pos=2,adj=c(2,2)) 
	})


plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=8,w=8)
plot(p1)
plot(p2)
quit()


