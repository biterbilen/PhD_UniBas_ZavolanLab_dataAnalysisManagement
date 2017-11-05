#source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))
library(latticeExtra)

#x=ln(X)
#y=ln(Y)
#ln(X+Y)=ln(exp(x)+exp(y))=x+ln(1+exp(y-x))
#http://jblevins.org/notes/log-sum-exp
#x should be bigger than y
log_form_add <- function (x, y, z=NULL) {
	k <- NULL
	if (x > y) {
		k <- x+log(1+exp(y-x));
	}
	else {
		k <- y+log(1+exp(x-y));
	}
	if (is.null(z)) {
		return(k);
	} else if (length(z) == 1) {
		return(log_form_add(k, z, z=NULL));
	} else {
		#TODO this part does not work
		return(log_form_add(k, z[1],z=z[2:length(z)]));
	}
}

# Given log(x) and log(y) compute log(x - y). 
# http://snipplr.com/view/12707/
#CHECK
log_form_subtract <- function (x,y) { 
	if (x <= y) {
		return (0.0/0.0);
	}

	diff <- x - y;
	if (diff > 708.0) {
		return (x);
	}
	#We use the following trick:
	#log(x-y) = log(a(x/a - y/a)) = log(a) + log(x/a - y/a)
	#         = log(a) + log(exp(log(x)-log(a)) - exp(log(y)-log(a)))
	#We pick log(a) = (log(x) + log(y)) / 2. So
	#log(x-y) = (log(x) + log(y)) / 2 +
	#         = log(exp((log(x) - log(y)) / 2) - exp((log(y) - log(x)) / 2))
	diff <- diff / 2; 
	return (x + y) / 2.0 + log(exp(diff) - exp(-diff));
}

plotbw <- function(x,ofile,tag,lab1,lab2,x2=NULL,bin=10,yi=3,xi=1,plotfreq=T) {
	p <- bwplot(x[,yi] ~ factor(floor(x[,xi]/bin)),
		par.settings = article.theme, data2=x2,
		type=c("g"),scale=list(y=list("log"=F), x=list("rot"=90,cex=article.cex)), 
		xlab=paste(lab2,"/",bin,sep=""), 
		pch = "|",
		ylab=lab1,main=tag,cex=article.cex,fontfamily=article.fontfamily,
		panel = function(x, y, ...,data2){
			panel.superpose
			panel.bwplot(x, y, ...)
			panel.grid(h = -1, v = -1) 
			if (!is.null(data2)) {
				fid <- c("paramofGroupMedians","paramofGroupMinimums","paramofGroupMaximums")
				flab <- c("{Group Medians}", "{Group Minimums}", "{Group Maximums}")
				flab <- gsub("^",data2$param,flab)
				for (i in 1:length(fid)) 
					if (length(grep(fid[i], names(data2))>0)) 
						panel.ablineq(h=data2[fid[i]],lwd=article.lwd,label=flab[i],adj=c(0,-0.1),cex=1,fontfamily=article.fontfamily,at=0.7)
				fid <- c("FittedMean","FittedLowerBound","FittedUpperBound")
				flab <- c("Fitted Mean", "Fitted p=Lower Bound", "Fitted p=Upper Bound")
				flab <- gsub("Lower",data2$ql,flab)
				flab <- gsub("Upper",data2$qu,flab)
				for (i in 1:length(fid)) 
					if (length(grep(fid[i], names(data2))>0)) 
						panel.ablineq(h=data2[fid[i]],lwd=article.lwd,label=flab[i],adj=c(0,-0.1),cex=1,fontfamily=article.fontfamily,at=0.1,col.text="red",col="red")
			}
			#copy number frequencies for each boxplot
			if (plotfreq == T) {
				freq <- table(x)
				panel.points(x, freq[x]/max(freq), pch="o", cex=article.cex, col="gold3")
				label <- paste("[Count of {",lab2,"/",bin,"}]",sep="");
				panel.text(min(as.numeric(x)), max(y), adj=c(0,-1), paste(label, "/argmax_i", label, "_i", sep=""), cex=1, col="gold3",fontfamily=article.fontfamily)
			}
		})
	plot.new(); pdf(paste(ofile,tag,".pdf",sep=""),h=6,w=8);
	print(p)
	dev.off()
}

#The low coverage positions are hypermutated
#The beta function parameters (i.e. effective number of mutations not mutations) are regularized s.t.
#T=aa+bb(min of the least (quantile 0.25-075) deviating copy number)
#aa=muRef*T 
#bb=T-aa

#Gets the median of Group Minimums and median of Group Maximums and estimates a beta distribution with these
#Calculates the Probability that each position is coming from the estimated beta distribution [0.05-0.95]
#tail is probability calculated in
get_pbeta <- function (x, ofile, tag, lab1="", lab2="", p.tail="both",bin=10,paramfile=NULL) {
#	x <- x[s4SNP,]
	#find emprical distribution paramaters grouped by copy number
	x$V3 <- x$V1 / x$V2
	fiven <- tapply(x$V3, x$x, fivenum, na.rm=T)
	fivenm <- matrix(unlist(fiven),byrow=T,ncol=5,dimnames=list(names(fiven)))
	fivenm_params <- apply(fivenm,2,median, na.rm=T) 
	linetype <- "Median"
	if (fivenm_params[5] == 1 | fivenm_params[1] == 0) {
		fivenm_params <- apply(fivenm,2,mean, na.rm=T) 
		linetype <- "Mean"
	}

	#select trusted set for estimation of the beta distribution parameters
	library(MASS)
	minCopyNumber <- 20 #hard coded ???????????
#	sp <- x$V3 < fivenm_params[5] & x$V3 > fivenm_params[1] & x$V2 > minCopyNumber
	s <- which(x$V3 < fivenm_params[5] & x$V3 > fivenm_params[1] & x$V2 > minCopyNumber)

	shape1 <- 1
	shape2 <- shape1*(1-fivenm_params[3])/fivenm_params[3]
	
	if (length(s) > 200) {
		#beta dist params are estimated by beta.select function from LearnBayes package; 2013.03.27; Biter
		library(LearnBayes)
		probs <- c(0.5,0.9)
		sq <- as.numeric(levels(factor(quantile(x[s,"V3"], probs=probs))))
		q1 <- list(p=probs[1],x=sq[1])
		q2 <- list(p=probs[2],x=sq[2])
		distribution <- NULL
		distribution$estimate <- beta.select(q1,q2)

#		x2_3 <- as.data.frame(t(cbind(shape1=distribution$estimate[1], shape2=distribution$estimate[2])))
#		names(x2_3) <- "value"
#		x2_3$type <- "learnBayes"
#		x2_3$label <- tag 
#		write.table(x2_3[c("label", "type","value")], quote=F, sep="\t", file = paramfile, row.names=F, col.names=F,append=T)

#		#depracated way of estimating beta desnsity params	
#		res <- try( distribution <- fitdistr(x[s,"V3"],"beta",start=list(shape1=shape1,shape2=shape2)),T)
#		while(length(grep("Error",res)>0)) {
#			print("fitdistr try 2")
#			ss <- sample(s, length(s)-1)
#			shape2=1
#			while(length(grep("Error",res)>0) & shape2<500) {
#				print("fitdistr try 22")
#				res <- try( distribution <- fitdistr(x[s,"V3"],"beta",start=list(shape1=shape1,shape2=shape2)),T)
#				shape2 <- shape2 + 1
#			}
#		}
	} else  {
		distribution <- data.frame(estimate=c(shape1,shape2))
	}

	shape1 <- distribution$estimate[1]
	shape2 <- distribution$estimate[2]
	mu <- shape1/(shape1+shape2)
	ql <- 0.001
	qu <- 0.999
	muRef <- qbeta(c(ql,qu),shape1,shape2)
#lab1="haha"
#p.tail="lower"
	x2 <- list(label=lab1, p.tail=p.tail,
					ql=ql,qu=qu,
					param=linetype,paramofGroupMedians=fivenm_params[3],paramofGroupMinimums=fivenm_params[1],paramofGroupMaximums=fivenm_params[5],
					minCopyNumber=minCopyNumber, TotalN=nrow(x),
					FittedN=length(s[s]),FittedShape1=shape1,FittedShape2=shape2,FittedMean=mu,FittedLowerBound=muRef[1],FittedUpperBound=muRef[2] )

	xi <- grep("x", names(x));
	yi <- grep("V3", names(x));
	plotbw(x, ofile,paste(".",tag,".ratio",sep=""),lab1,lab2,x2=x2,xi=xi,yi=yi,bin=bin);

	x$pbetal <- pbeta(muRef[1],x$V1+shape1,x$V2-x$V1+shape2,log.p=T,lower.tail=T)
	x$pbetau <- pbeta(muRef[2],x$V1+shape1,x$V2-x$V1+shape2,log.p=T,lower.tail=F)

	#NA check
	x[which(is.na(x$pbetal)),"pbetal"] <- median(x$pbetal,na.rm=T)
	x[which(is.na(x$pbetau)),"pbetau"] <- median(x$pbetau,na.rm=T)

	if (p.tail == "both") {
		x$pbetacombined <- mapply(log_form_add, x$pbetau, x$pbetal)
	} else if(p.tail == "lower") {
		x$pbetacombined <- x$pbetal
	} else if(p.tail == "upper") {
		x$pbetacombined <- x$pbetau
	}

	xi <- grep("x", names(x));
	yi <- grep("pbetacombined", names(x));
	plotbw(x, ofile,paste(".",tag,".score",sep=""),paste("Background Probability {",lab1,"}",sep=""), lab2, plotfreq=F, xi=xi, yi=yi,bin=bin);

	mdn <- median(x$pbetacombined)
	x2_2 <- as.data.frame(t(cbind(data.frame(x2),Mdn=mdn, NlessThanMdn=length(which(x$pbetacombined < mdn)))))
	names(x2_2) <- "value"
	x2_2$type <- rownames(x2_2)
	x2_2$label <- tag 
	write.table(x2_2[c("label", "type","value")], quote=F, sep="\t", file = paramfile, row.names=F, col.names=F,append=T)

	return(x$pbetacombined)
}

cargs <- commandArgs(trailingOnly = T)
#cargs <- c("1.00.bed12.gz", "1.00", "/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/", "tmp.1.00.params.txt", "tmp.1.00.bedplus")
#cargs <- c("0.35.bed12.gz", "tmp","/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","tmp");
#cargs <- c("309.bed12.gz","309","/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","all.309.params.txt","all.309.bedplus");

paramfile <- paste(cargs[2],".params.txt", sep="")
outfile <- paste(cargs[2],".bedplus", sep="")

if (length(cargs)>2) source(file.path(cargs[3],"theme.R"))
if (length(cargs)>3) paramfile <- cargs[4]
if (length(cargs)>4) outfile <- cargs[5]

#read data
a <- read.table(gzfile(cargs[1]), sep="\t");
names(a) <- c("seqname","start","end","name","mut","strand","slop10_mut","copies","all_copies","background_mut","background_copies", "SNPdb");

#1.
x <- data.frame(matrix(ncol = 0 , nrow = nrow(a)))
x$x <- a$copies
x$V1 <- a$slop10_mut - a$mut
x$V2 <- a$slop10_mut
tag <- "NeighbourhoodNotCrosslinked.Probability"
a[,tag] <- get_pbeta(x, cargs[2],tag,"Position Neigbourhood Mutations / All Mutations", "Selected Tags",p.tail="lower",paramfile=paramfile)

#2.
x <- data.frame(matrix(ncol = 0 , nrow = nrow(a)))
x$x <- a$copies
x$V1 <- a$copies
x$V2 <- a$all_copies
tag <- "Contamination.Probability"
a[,tag] <- get_pbeta(x, cargs[2],tag,"Selected Tags / All Tags", "Selected Tags",p.tail="lower",paramfile=paramfile)

#3.
x <- data.frame(matrix(ncol = 0 , nrow = nrow(a)))
x$x <- a$copies
x$V1 <- a$background_mut
x$V2 <- a$background_copies
tag <- "SNP.Probability"
s4SNP <- which(x$V1 > 0) #at least once mutated
a[s4SNP,tag] <- get_pbeta(x[s4SNP,], cargs[2],tag,"Position Mutations / Selected Tags [RNAseq]", "Selected Tags",p.tail="upper",paramfile=paramfile)
a[-s4SNP,tag] <- min(a[s4SNP,tag]) #not calculated ones are set to the mimimum of the calculated ones (i.e. NOT SNP)

#4.
x <- data.frame(matrix(ncol = 0 , nrow = nrow(a)))
x$x <- a$copies
x$V1 <- a$mut
x$V2 <- a$copies
tag <- "NotCrosslinkedJustLower.Probability"
a[,tag] <- get_pbeta(x, cargs[2],tag,"Position Mutations / Selected Tags", "Selected Tags",p.tail="lower",paramfile=paramfile)

#5.
x <- data.frame(matrix(ncol = 0 , nrow = nrow(a)))
x$x <- a$copies
x$V1 <- a$mut
x$V2 <- a$copies
tag <- "NotCrosslinked.Probability"
a[,tag] <- get_pbeta(x, cargs[2],tag,"Position Mutations / Selected Tags", "Selected Tags",paramfile=paramfile)

SNPdbfilter <- as.numeric(a$SNPdb == T)
NeighbourhoodNotCrosslinkedFilter <- as.numeric(a$NeighbourhoodNotCrosslinked.Probability >= median(a$NeighbourhoodNotCrosslinked.Probability))
ContaminationFilter <- as.numeric(a$Contamination.Probability >= median(a$Contamination.Probability))
SNPfilter <- as.numeric(a$SNP.Probability > median(a[s4SNP,"SNP.Probability"])) #here the equality is suppressed
a$filterScore <-  apply(matrix(as.numeric(c(SNPdbfilter, NeighbourhoodNotCrosslinkedFilter,ContaminationFilter,SNPfilter )),byrow=F,nrow=nrow(a)),1,sum) 
#a$filterScore <-  mapply(paste,SNPdbfilter, NeighbourhoodNotCrosslinkedFilter,ContaminationFilter,SNPfilter, sep="");
a$Crosslinked <- a$NotCrosslinked.Probability < median(a$NotCrosslinked.Probability)

#sort and rank
a <- a[with(a, order(-NotCrosslinked.Probability, decreasing=T)),]
a$rank <- 1:nrow(a)

fields <- c(names(a)[1:3],"rank", "NotCrosslinked.Probability",names(a)[6], "Crosslinked", "filterScore", names(a)[5],names(a)[7:(length(names(a))-4)],"name")
write.table(a[fields], file=outfile, quote=F, sep="\t", row.names=F, col.names=T)
system(paste("gzip -f", outfile, sep = " "))
#---------------------------------------------------------------------------------------------

plot.new(); pdf(paste(outfile,".pdf",sep=""),h=8,w=8)
p1 <- histogram(~NotCrosslinked.Probability, data=a,nint=30)
p2 <- histogram(~NeighbourhoodNotCrosslinked.Probability, data=a,nint=30)
p3 <- histogram(~Contamination.Probability, data=a,nint=30)
p4 <- histogram(~SNP.Probability, data=a,nint=30)
plot(p1, split = c(1, 1, 2, 2)) 
plot(p2, split = c(2, 1, 2, 2), newpage = FALSE) 
plot(p3, split = c(1, 2, 2, 2), newpage = FALSE) 
plot(p4, split = c(2, 2, 2, 2), newpage = FALSE) 
dev.off()

#---------------------------------------------------------------------------------------------
quit()
#---------------------------------------------------------------------------------------------

#plot mutation/copies vs copies
tp1 <- bwplot((mut/copies) ~ factor(copies), data=a, 
	par.settings = article.theme,
	subset=s, drop=T, type=c("g"),scale=list(y=list("log"=F), x=list("rot"=60,cex=0.4)), 
	panel = function(x, y, fill="blue", ...){
		panel.bwplot(x, y,...)
		panel.grid(h = -1, v = -1)
	})
tp2 <- bwplot(pbeta ~ factor(copies), data=a, 
	par.settings = article.theme,
	subset=s,  drop=T, type=c("g"),scale=list(y=list("log"=F),x=list("rot"=60,cex=0.4)), 
	panel = function(x, y, ...){
		panel.points(mean(x), y, fill="blue", ...)
		panel.bwplot(x, y,...)
		panel.grid(h = -1, v = -1)
	})

ofile <- paste(cargs[2], ".bwplot4mutation_rate_v_copies.pdf", sep="")
plot.new(); pdf(ofile,h=8,w=8);
plot(tp1, split = c(1, 1, 1, 2)) 
plot(tp2, split = c(1, 2, 1, 2), newpage = FALSE) 
dev.off()

