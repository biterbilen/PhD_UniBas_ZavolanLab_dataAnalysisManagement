library(lattice)

cargs <- commandArgs()
#cargs <- c(1,1,1,"Analysis/Reproducibility/TIA1.gtf.gz");

a <- read.table(cargs[4],header=F, sep="", stringsAsFactors=F);
lib1 <- a[1,21];
lib2 <- a[1,24];
header <- c("chr","src","fea","beg","end","scr","str","frm","gene_id","gene_id_v","d","transcript_id","transcript_id_v","d","lib_count","lib_count_v","d","regularized_average_feature_value","regularized_average_feature_value_v","d",lib1,"lib1_v","d",lib2,"lib2_v","d");
names(a) <- header;

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)      ## set as default

feas <- unique(a$fea)
scrs <- sort(unique(a$scr), decreasing=T);

newheader <- c("reproducedPosProp", "fea", "scr")
aplot <- data.frame(matrix(nrow=length(feas)*length(scrs),ncol=length(newheader)))
names(aplot) <- newheader
k <- 1;
for (fea in feas) {
	for (scr in scrs) {
		suba <- a[a$scr == scr & a$fea == fea,c("lib_count_v")];
		aplot[k,] <- c(length(which(suba==2))/length(suba), fea, scr);
		k <- k+1;
	}
}

#write.table(aplot,quote=F,file=paste(cargs[4],'.cum', sep=''),sep="\t",row.names=F);

#positional reproducibility
tit = paste("Reproduced Features", gsub('.gtf.gz', '', basename(cargs[4])), sep='-');
plot.new(); pdf(paste(cargs[4], '.posrep','.pdf',sep=''), width=8, height=8)
barchart( 100*as.numeric(reproducedPosProp) ~ fea, group = as.numeric(scr), data=aplot, origin = 0,
	auto.key = list(title = "Minimum Feature Frequency", space="bottom",columns=3, cex=1.5,font=2),
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

