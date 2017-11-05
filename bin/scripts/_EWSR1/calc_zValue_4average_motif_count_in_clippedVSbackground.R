library(latticeExtra)

cargs <- commandArgs();
#cargs <- c(1,1,1,"../Analysis/RecombHotspot/all.motiffreq");

a <- read.table(cargs[4], header=F);
ci <- which(grepl('CLIPPED', a$V1));
nci <- which(grepl('NOTCLIPPED', a$V1));

a$V5[1:dim(a)[1]] <- "GENOME";
a$V5[ci] <- "CLIPPED";
a$V5[nci] <- "NOTCLIPPED";

names(a) <- c("id", "motif", "motifFreq", "seqLen", "type");

motifs <- unique(a$motif)

countci <- length(which(a$type == "CLIPPED")) / length(motifs);
countnci <- length(which(a$type == "NOTCLIPPED")) / length(motifs);
countgi <- length(which(a$type == "GENOME")) / length(motifs);

types <- unique(a$type)

a$M <- 0;
a$N <- 0;
m <- "CCCCACCCC";
for (motif in motifs) {
	for (t in types) {
		i <- which(a$type==t & a$motif==motif);
		a$M[i] <- sum(as.numeric(a$motifFreq[i]));
		a$N[i] <- sum(as.numeric(a$seqLen[i]) - nchar(as.character(motif)));
	}
	i <- which(a$type=="CLIPPED" & a$motif==motif);
	x <- a$M[i[1]]
	N <- a$N[i[1]]
	t <- "NOTCLIPPED";
	for (t in c("NOTCLIPPED", "GENOME")) {
		i0 <- which(a$type==t & a$motif==motif);
		x0 <- a$M[i0[1]]
		y0 <- a$N[i0[1]]
		p <- x0 / y0;
		#P(x|N,p)
		m <- N * p; 
		v <- N * p * (1 - p);
		z <- (x - m) / sqrt(v);   #/ sqrt(y);
		tab <- list(foregroundVbackground=paste("CLIPPEDv",t,sep=""),motif=motif,"z-value"=z,x=x,N=N,p=p,m=m,v=v);
#		tab.fmt <- format.data.frame(tab,width=3,trim=T,scientific=F);
		write.table(tab, sep="\t",col.names=T, row.names=F, quote=F);
	}
}

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4,cex=2)      ## set as default

plot.new(); pdf(paste(cargs[4],'.pdf',sep=''), h=8, w=8);
densityplot(~ log(motifFreq)-log(seqLen-nchar(as.character(motif)))|motif, data=a, groups=type, auto.key = list(space = "top", cex=1), 
	par.settings=list(axis.text=list(font=2,cex=1.0),
		par.strip.text=list(font=2,cex=1.0),
		par.ylab.text=list(font=2,cex=1.0),
		par.xlab.text=list(font=2,cex=1.0),
		par.main.text=list(font=2,cex=1.0)),
	lwd=2,
	cex=0.5,
	#layout=c(length(motifs)/length(types),length(types)),
	xlab="Observed motif probability (log)",plot.points=F) +
layer(panel.key(c(sprintf('N(CLIPPED Superclusters) = %d\nN(NOTCLIPPED Genes) = %d\nN(GENOME overall) = %d\n', countci, countnci, countgi)),lwd=2,cex=0.75))
dev.off();


