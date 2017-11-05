library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))
#cargs <- c("aGenePerRank/308.309.genePerRank", "title")
cargs <- c("a", "11907");
cargs <- commandArgs(trailingOnly=T)

a <- read.table(gzfile(cargs[1]), header=T, sep="\t")
totalGenes <- as.numeric(cargs[2]);

a <- a[a$rank%%500==0,]
dim(a) 
a$nf1_gene_count <- a$f1_gene_count - a$common_gene_count
a$nf2_gene_count <- a$f2_gene_count - a$common_gene_count
a$nf1f2_gene_count <- totalGenes - (a$f1_gene_count + a$f2_gene_count - a$common_gene_count)

fun <- function(x) { fisher.test(matrix(x,nr=2),alternative="g")$estimate }
a$estimate <- apply(a[c("common_gene_count","nf1_gene_count","nf2_gene_count","nf1f2_gene_count")], 1,fun)

fun <- function(x) { -log(fisher.test(matrix(x,nr=2),alternative="g")$p.value) }
a$enrichment <- apply(a[c("common_gene_count","nf1_gene_count","nf2_gene_count","nf1f2_gene_count")], 1,fun)

#fun <- function(x) { phyper(x[1]-1,x[2],x[3],x[4],log.p=F,lower.tail=F) }
#a$enrichment2 <- apply(a[c("common_gene_count","f1_gene_count","nf1_gene_count","f2_gene_count")], 1,fun)
head(a)

write.table(a,file=paste(cargs[1],".txt",sep=""),sep="\t",quote=F);

s <- is.finite(a$enrichment);
m <- max(a[s,"enrichment"])
a$enrichment[!s] <- m
max(a$enrichment)
max(a$enrichment)
#s <- is.finite(a$enrichment2);
#m <- max(a[s,"enrichment2"])
#a$enrichment2[!s] <- m
#max(a$enrichment2)
#min(a$enrichment2)


plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
#xyplot(enrichment ~ cumsum(estimate), data=a,
#	xlab="Odds Ratio",
#	ylab="Gene Overlap Enrichment (-log(p.value))",
#	scale=list(y=list("log"=F), x=list(log=F,cex=article.cex)), 
#xyplot(common_gene_count / (f1_gene_count+f2_gene_count-common_gene_count) ~ rank, data=a,
#	xscale.components = xscale.components.log10ticks,
#	scale=list(y=list("log"=F), x=list(log=F,cex=article.cex)), 
xyplot(enrichment ~ rank, data=a,
	xscale.components = xscale.components.log10ticks,
	scale=list(y=list("log"=F), x=list(log=T,cex=article.cex)), 
	panel = function(x,y,...,col=col) {
		panel.grid(h = -1, v = -1)
		panel.xyplot(x, y, ...)
	},
	par.settings = article.theme, lwd=article.lwd, cex=article.cex, pch='o',
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font)
	) # +
#	as.layer(xyplot(enrichment ~ rank, data=a, col="red"), y.same=F)
dev.off()

quit()

p1 <- xyplot((common_gene_count/total_gene_count) ~ rank | type, data=a, #[a$total_gene_count%%10==0,], 
	layout=c(1,length(unique(a$type))), main=paste("Recovered Genes per Rank", cargs[2],sep=' - '),
	ylab="Common Genes per Recovered Genes", xlab="Rank of Crosslinked Position",
	auto.key = list(space="top",points=F, col=article.colors,columns=length(unique(a$tag)), cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	group=tag, origin=0,
	ylim=c(0,1),
	scale=list(y=list("log"=F), x=list(log=T,cex=article.cex)), 
	xscale.components = xscale.components.log10ticks,
	panel = function(x,y,...,col=col) {
		panel.grid(h = -1, v = -1)
		panel.xyplot(x, y, ...)
		panel.abline(h=0.5, col="gray",...)
#		panel.text(1,0.8,"wo_filters",col=col)
#		panel.xyplot(x, y, ...,col=col)
	},
	par.settings = article.theme, lwd=article.lwd, cex=article.cex, pch='o',
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font))


plot.new(); png(paste(cargs[1],'.png',sep=''), h=720,w=720)
print(p1)
dev.off()

if (length(cargs) > 2 & cargs[3] == "pdf") {
	plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
	print(p1)
	dev.off()
}
