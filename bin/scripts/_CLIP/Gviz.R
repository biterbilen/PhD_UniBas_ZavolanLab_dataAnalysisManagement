### R code from vignette source 'vignettes/Gviz/inst/doc/Gviz.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: init
###################################################
options(width=65)
library(xtable)
source(system.file("scripts/documentation.R", package="Gviz"))
xtabDetails <- details
addParTable <- function(class, skip=c("showTitle", "size", "background.title"), add=NULL)
{
    Parameters <- data.frame("Display Parameter"=names(xtabDetails[[class]]),
                             "Description"=xtabDetails[[class]], check.names=FALSE)
    align <- "lrp{5in}"
    if(!is.null(add)){
        Parameters <- cbind(Parameters, add)
        align <- c("lp{4in}", "lp{4in}", "lp{4in}", "lp{4in}")
    }
    Parameters <- Parameters[order(Parameters[,1]),]
    Parameters <- apply(Parameters, 2, function(x) gsub("_", "\\_", x, fixed=TRUE))
    rownames(Parameters) <-  gsub("_", "\\_", rownames(Parameters), fixed=TRUE)
    sel <- Parameters[,1] %in% skip
    Parameters[,2] <- gsub("\\code{\\linkS4class{", "\\Rclass{{", Parameters[,2], fixed=TRUE)
    print(xtable(Parameters[!sel,], align=align), sanitize.text.function=function(x) x, include.rownames=FALSE,
          floating=FALSE, tabular.environment="longtable")
    return(invisible())
}

hasUcscConnection <- !is(try(rtracklayer::browserSession(), silent=TRUE), "try-error")
hasBiomartConnection <- !is(try(biomaRt::listMarts(), silent=TRUE), "try-error")



###################################################
### code chunk number 2: loadPackage
###################################################
library(Gviz)


###################################################
### code chunk number 3: AnnotationTrack
###################################################
library(GenomicRanges)
data(cpgIslands)
class(cpgIslands)
chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name="CpG")


###################################################
### code chunk number 4: plotAnnotationTrack
###################################################
plotTracks(atrack)


###################################################
### code chunk number 5: GenomeAxisTrack
###################################################
gtrack <- GenomeAxisTrack()


###################################################
### code chunk number 6: plotGenomeAxisTrack
###################################################
plotTracks(list(gtrack, atrack))


###################################################
### code chunk number 7: showIdeogramTrack (eval = FALSE)
###################################################
## itrack <- IdeogramTrack(genome=gen, chromosome=chr)


###################################################
### code chunk number 8: doIdeogramTrack
###################################################
if(hasUcscConnection){
    itrack <- IdeogramTrack(genome=gen, chromosome=chr)
}else{
    data(itrack)
}


###################################################
### code chunk number 9: plotIdeogramTrack
###################################################
plotTracks(list(itrack, gtrack, atrack))


###################################################
### code chunk number 10: GeneRegionTrack
###################################################
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))


###################################################
### code chunk number 11: zooming
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack), from=25e6, to=28e6)


###################################################
### code chunk number 12: zooming2
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome=chr)
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), from=26450430, to=26450490, cex=0.8)


###################################################
### code chunk number 13: DataTrack
###################################################
set.seed(255)
lim <- c(26463500, 26495000)
coords <- sort(c(lim[1], sample(seq(from=lim[1], to=lim[2]), 99), lim[2]))
dat <- runif(100, min=-10, max=10)
dtrack <- DataTrack(data=dat, start=coords[-length(coords)], end=coords[-1], chromosome=chr, 
                    genome=gen, name="Uniform")
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from=lim[1], to=lim[2])


###################################################
### code chunk number 14: DataTrackHist
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from=lim[1], to=lim[2], type="histogram")


###################################################
### code chunk number 15: displayPars1
###################################################
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="Gene Model", 
                           showId=TRUE, background.title="brown")
head(displayPars(grtrack))
displayPars(grtrack) <- list(background.panel="#FFFEDB")
head(displayPars(grtrack))
plotTracks(list(itrack, gtrack, atrack, grtrack), from=lim[1], to=lim[2])


###################################################
### code chunk number 16: displayPars2
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack), from=lim[1]-1000, to=lim[2], 
           background.panel="#FFFEDB", background.title="darkblue")


###################################################
### code chunk number 17: displayPars3
###################################################
dp <- availableDisplayPars(grtrack)
tail(dp)


###################################################
### code chunk number 18: GenomeAxisTrackClass1
###################################################
axisTrack <- GenomeAxisTrack()
plotTracks(axisTrack, from=1e6, to=9e6)


###################################################
### code chunk number 19: GenomeAxisTrackClass2
###################################################
axisTrack <- GenomeAxisTrack(range=IRanges(start=c(2e6, 4e6), end=c(3e6, 7e6)))
plotTracks(axisTrack, from=1e6, to=9e6)


###################################################
### code chunk number 20: GenomeAxisTrackClass3
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, add53=TRUE, add35=TRUE)


###################################################
### code chunk number 21: GenomeAxisTrackClass4
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, add53=TRUE, add35=TRUE, littleTicks=TRUE)


###################################################
### code chunk number 22: GenomeAxisTrackClass5
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, exponent=4)


###################################################
### code chunk number 23: GenomeAxisTrackClass6
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, labelPos="below")


###################################################
### code chunk number 24: GenomeAxisTrackClass7
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, scale=0.5)


###################################################
### code chunk number 25: GenomeAxisTrackClass8
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, scale=0.5, labelPos="below")


###################################################
### code chunk number 26: GenomeAxisTrackClassTable
###################################################
addParTable("GenomeAxisTrack")


###################################################
### code chunk number 27: IdeogramTrackClass1Show (eval = FALSE)
###################################################
## ideoTrack <- IdeogramTrack(genome="hg19", chromosome="chrX")
## plotTracks(ideoTrack, from=85e6, to=129e6)


###################################################
### code chunk number 28: IdeogramTrackClass1Do
###################################################
if(hasUcscConnection){
    ideoTrack <- IdeogramTrack(genome="hg19", chromosome="chrX")
}else{
    data(itrack)
}
plotTracks(ideoTrack, from=85e6, to=129e6)


###################################################
### code chunk number 29: IdeogramTrackClass2
###################################################
plotTracks(ideoTrack, from=85e6, to=129e6, showId=FALSE)


###################################################
### code chunk number 30: IdeogramTrackClassTable
###################################################
addParTable("IdeogramTrack")


###################################################
### code chunk number 31: DataClass1
###################################################
data(twoGroups)
dTrack <- DataTrack(twoGroups, name="uniform")
plotTracks(dTrack)


###################################################
### code chunk number 32: <types
###################################################
types <- data.frame(Value=c("p", "l", "b", "a", "s", "S", "g", "r", "h", "smooth", "histogram", "mountain", "boxplot", "gradient", "heatmap"),
                    Type=c("dot plot", "lines plot", "dot and lines plot", "lines plot of average (i.e., mean) values", "stair steps (horizontal first)",
                           "stair steps (vertical first)", "add grid lines", "add linear regression line", "histogram lines", "add loess curve", 
                           "histogram (bar width equal to range with)", "'mountain-type' plot relative to a baseline", "box and whisker plot",
                           "false color image of the summarized values", "false color image of the individual values"))
print(xtable(types, align="lrp{5in}"), sanitize.text.function=function(x) x, include.rownames=FALSE,
          floating=FALSE, tabular.environment="longtable")


###################################################
### code chunk number 33: typePlots
###################################################
pushViewport(viewport(layout=grid.layout(nrow=8, ncol=2)))
i <- 1
for(t in types$Value)
{
    pushViewport(viewport(layout.pos.col=((i-1)%%2)+1, layout.pos.row=((i-1)%/%2)+1))
    names(dTrack) <- t
    plotTracks(dTrack, type=t, add=TRUE, cex.title=0.8, margin=0.5)
    i <- i+1
    popViewport(1)
}
popViewport(1)
names(dTrack) <- "uniform"


###################################################
### code chunk number 34: mutitype
###################################################
plotTracks(dTrack, type=c("boxplot", "a", "g"))


###################################################
### code chunk number 35: grouping
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("a", "p"))


###################################################
### code chunk number 36: typeGroupedPlots
###################################################
pushViewport(viewport(layout=grid.layout(nrow=7, ncol=1)))
i <- 1
for(t in c("a", "s", "smooth", "histogram", "boxplot", "heatmap"))
{
    pushViewport(viewport(layout.pos.col=((i-1)%%1)+1, layout.pos.row=((i-1)%/%1)+1))
    names(dTrack) <- t
    plotTracks(dTrack, type=t, add=TRUE, cex.title=0.8, groups=rep(1:2, each=3), margin=0.5)
    i <- i+1
    popViewport(1)
}
pushViewport(viewport(layout.pos.col=((7-1)%%1)+1, layout.pos.row=((7-1)%/%1)+1))
names(dTrack) <- "hor. hist."
plotTracks(dTrack, type="histogram", stackedBars=FALSE, add=TRUE, cex.title=0.8, groups=rep(1:2, each=3), margin=0.5)
popViewport(2)
names(dTrack) <- "uniform"


###################################################
### code chunk number 37: groupingLegend
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("a", "p"), legend=TRUE)


###################################################
### code chunk number 38: biggerdata
###################################################
dat <- sin(seq(pi, 10*pi, len=500))
dTrack.big <- DataTrack(start=seq(1,100000, len=500), width=15, chromosome="chrX",
                        genome="hg19", name="sinus", 
                        data=sin(seq(pi, 5*pi, len=500))*runif(500, 0.5, 1.5))
plotTracks(dTrack.big, type="hist")


###################################################
### code chunk number 39: aggregation
###################################################
plotTracks(dTrack.big, type="hist", window=50)


###################################################
### code chunk number 40: aggregation2
###################################################
plotTracks(dTrack.big, type="hist", window=-1, windowSize=2500)


###################################################
### code chunk number 41: transformation
###################################################
plotTracks(dTrack.big, type="l", transformation=function(x){x[x<0] <- 0; x})


###################################################
### code chunk number 42: groupingAv1
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("b"), aggregateGroups=TRUE)


###################################################
### code chunk number 43: groupingAv2
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("b"), aggregateGroups=TRUE, aggregation="max")


###################################################
### code chunk number 44: DataTrackClassTable
###################################################
addParTable("DataTrack")


###################################################
### code chunk number 45: anntrack1
###################################################
aTrack <- AnnotationTrack(start=c(10, 40, 120), width=15, chromosome="chrX",
                          strand=c("+", "*", "-"), 
                          id=c("Huey", "Dewey", "Louie"), genome="hg19", name="foo")
plotTracks(aTrack)


###################################################
### code chunk number 46: anntrack2
###################################################
plotTracks(aTrack, shape="box", showFeatureId=TRUE)


###################################################
### code chunk number 47: anntrack3
###################################################
plotTracks(aTrack, shape="ellipse", showFeatureId=TRUE, fontcolor="darkblue")


###################################################
### code chunk number 48: anntrack4
###################################################
aTrack.groups <- AnnotationTrack(start=c(50, 180, 260, 460, 860, 1240), width=c(15,20,40,100,200, 20), 
                                 chromosome="chrX",
                                 strand=rep(c("+", "*", "-"), c(1,3,2)), 
                                 group=rep(c("Huey", "Dewey", "Louie"), c(1,3,2)), 
                                 genome="hg19", name="foo")
plotTracks(aTrack.groups, showId=TRUE)


###################################################
### code chunk number 49: stacking1
###################################################
aTrack.stacked <- AnnotationTrack(start=c(50, 180, 260, 800, 600, 1240), width=c(15,20,40,100,500, 20), 
                                 chromosome="chrX",
                                 strand="*",
                                 group=rep(c("Huey", "Dewey", "Louie"), c(1,3,2)), 
                                 genome="hg19", name="foo")
plotTracks(aTrack.stacked, showId=TRUE)


###################################################
### code chunk number 50: stacking2
###################################################
plotTracks(aTrack.stacked, stacking="dense")


###################################################
### code chunk number 51: features
###################################################
feature(aTrack.stacked)
feature(aTrack.stacked)[1:4] <- c("foo", "bar", "bar", "bar")


###################################################
### code chunk number 52: featuresPlot
###################################################
plotTracks(aTrack.stacked, showId=TRUE, foo="darkred", bar="darkgreen")


###################################################
### code chunk number 53: overplotting
###################################################
data("denseAnnTrack")
plotTracks(denseAnnTrack, showOverplotting=TRUE)


###################################################
### code chunk number 54: collapse1
###################################################
data(collapseTrack)
plotTracks(ctrack, extend.left=1800)


###################################################
### code chunk number 55: collapse2
###################################################
plotTracks(ctrack, extend.left=1800, min.width=1)


###################################################
### code chunk number 56: collapse3
###################################################
plotTracks(ctrack, extend.left=1800, min.width=1, collapse=TRUE)


###################################################
### code chunk number 57: collapse4
###################################################
plotTracks(ctrack, extend.left=1800, min.width=3, min.distance=5, collapse=TRUE)


###################################################
### code chunk number 58: collapse5
###################################################
plotTracks(ctrack, extend.left=1800, min.width=3, min.distance=5, collapse=TRUE,
           mergeGroups=TRUE)


###################################################
### code chunk number 59: AnnotationTrackClassTable
###################################################
addParTable("AnnotationTrack")


###################################################
### code chunk number 60: generegtrack
###################################################
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="foo")
head(gene(grtrack))
head(transcript(grtrack))
head(exon(grtrack))
head(symbol(grtrack))
plotTracks(grtrack, extend.left=20000, showId=TRUE)


###################################################
### code chunk number 61: generegtrack2
###################################################
plotTracks(grtrack, extend.left=20000, showId=TRUE, geneSymbols=FALSE)


###################################################
### code chunk number 62: generegtrack3
###################################################
plotTracks(grtrack, collapseTranscripts=TRUE, shape="arrow", showId=TRUE, extend.left=20000)


###################################################
### code chunk number 63: tdb2grt1
###################################################
library(GenomicFeatures)
samplefile <- system.file("extdata", "UCSC_knownGene_sample.sqlite", package="GenomicFeatures")
txdb <- loadDb(samplefile)
GeneRegionTrack(txdb)


###################################################
### code chunk number 64: tdb2grt2
###################################################
txTr <- GeneRegionTrack(txdb, chromosome="chr6", start=300000, end=350000)


###################################################
### code chunk number 65: generegtrack4
###################################################
feature(txTr)
plotTracks(txTr, showId=TRUE, extend.left=1000)


###################################################
### code chunk number 66: GeneRegionTrackClassTable
###################################################
addParTable("GeneRegionTrack")


###################################################
### code chunk number 67: BiomartGeneRegionTrackShow (eval = FALSE)
###################################################
## biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
##                                   name="ENSEMBL")
## plotTracks(biomTrack)


###################################################
### code chunk number 68: BiomartGeneRegionTrackDo
###################################################
if(hasBiomartConnection){
    biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
                                        name="ENSEMBL")
}else{
    data("biomTrack")
}
plotTracks(biomTrack)


###################################################
### code chunk number 69: BiomartGeneRegionTrackCol
###################################################
plotTracks(biomTrack, col.line=NULL)


###################################################
### code chunk number 70: BiomartGeneRegionTrackClassTable
###################################################
addInfo <- t(data.frame(displayPars(biomTrack, names(details[["BiomartGeneRegionTrack"]]))))
colnames(addInfo) <- "Color"
addParTable("BiomartGeneRegionTrack", add=addInfo)


###################################################
### code chunk number 71: DetailsAnnotationTrack1
###################################################
library(GenomicRanges)
probes <- GRanges(seqnames="chr7", ranges=IRanges(start=c(2000000, 2070000, 2100000, 2160000), end=c(2050000, 2130000, 2150000, 2170000)),
                  strand=c("-", "+", "-", "-"))


###################################################
### code chunk number 72: DetailsAnnotationTrack2
###################################################
methylation <- matrix(c(rgamma(400, 1)), ncol=100, dimnames=list(paste("probe", 1:4, sep=""), NULL))
methylation[,51:100] <- methylation[,51:100] + 0:3
sgroups <- rep(c("grp1","grp2"), each=50)


###################################################
### code chunk number 73: DetailsAnnotationTrack3
###################################################
library(lattice) 
details <- function(identifier, ...) {
    d <- data.frame(signal=methylation[identifier,], group=sgroups)
    print(densityplot(~signal, group=group, data=d, main=list(label=identifier, cex=0.7),
                      scales=list(draw=FALSE, x=list(draw=TRUE)), ylab="", xlab="",
                      ), newpage=FALSE, prefix="plot")
}


###################################################
### code chunk number 74: DetailsAnnotationTrack4
###################################################
deTrack <- AnnotationTrack(range=probes, genome="hg19", chromosome=7, id=rownames(methylation),
                           name="probe details", stacking="squish",
                           fun=details)
plotTracks(deTrack)


###################################################
### code chunk number 75: DetailsAnnotationTrack5
###################################################
selFun <- function(identifier, start, end, track, GdObject, ...){
    gcount <- table(group(GdObject))
    ## This computes the width of 2 pixels in genomic coordinates
    pxRange <- Gviz:::.pxResolution(min.width=20, coord="x")
    return((end-start)<pxRange && gcount[identifier]==1)
}


###################################################
### code chunk number 76: DetailsAnnotationTrack6
###################################################
detFun <- function(identifier, GdObject.original, ...){
    plotTracks(list(GenomeAxisTrack(scale=0.3, size=0.2, cex=0.7), GdObject.original[group(GdObject.original)==identifier]), 
               add=TRUE, showTitle=FALSE)
}


###################################################
### code chunk number 77: DetailsAnnotationTrack7
###################################################
data(geneDetails)
deTrack2 <- AnnotationTrack(geneDetails, fun=detFun, selectFun=selFun,
                            groupDetails=TRUE, details.size=0.5, detailsConnector.cex=0.5, detailsConnector.lty="dotted",
                            shape=c("smallArrow", "arrow"), showId=TRUE)
plotTracks(deTrack2, extend.left=90000)


###################################################
### code chunk number 78: DetailsAnnotationTrack5
###################################################
plotTracks(deTrack, details.size=0.75, detailsConnector.pch=NA, detailsConnector.col="darkred", 
           detailsBorder.fill="#FFE3BF", detailsBorder.col="darkred", shape="box", detailsConnector.lty="dotted")


###################################################
### code chunk number 79: DetailsAnnotationTrackClassTableSec
###################################################
addParTable("DetailsAnnotationTrack")


###################################################
### code chunk number 80: SequenceTrack1
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
sTrack <- SequenceTrack(Hsapiens)
sTrack


###################################################
### code chunk number 81: SequenceTrack2
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20050)


###################################################
### code chunk number 82: SequenceTrack3
###################################################
fcol <- c(A="darkgray", C="darkgray", T="darkgray", G="darkgray")
plotTracks(sTrack, chromosome=1, from=20000, to=20050, fontcolor=fcol)


###################################################
### code chunk number 83: SequenceTrack4
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20050, add53=TRUE)


###################################################
### code chunk number 84: SequenceTrack5
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20050, add53=TRUE, complement=TRUE)


###################################################
### code chunk number 85: SequenceTrack6
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20100)


###################################################
### code chunk number 86: SequenceTrack7
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=201000)


###################################################
### code chunk number 87: SequenceTrack8
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20100, cex=0.5)


###################################################
### code chunk number 88: ucscTrack1 (eval = FALSE)
###################################################
## from <- 65921878
## to <- 65980988
## knownGenes <- UcscTrack(genome="mm9", chromosome="chrX", track="knownGene", from=from, to=to,
##                         trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
##                         symbol="name", transcript="name", strand="strand", fill="#8282d2", name="UCSC Genes")


###################################################
### code chunk number 89: ucscTrack2 (eval = FALSE)
###################################################
## refGenes <- UcscTrack(genome="mm9", chromosome="chrX", track="xenoRefGene", from=from, to=to,
##                       trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
##                       symbol="name2", transcript="name", strand="strand", fill="#8282d2",
##                       stacking="dense", name="Other RefSeq")
## 
## ensGenes <- UcscTrack(genome="mm9", chromosome="chrX", track="ensGene", from=from, to=to,
##                       trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
##                       symbol="name2", transcript="name", strand="strand", fill="#960000",
##                       name="Ensembl Genes")


###################################################
### code chunk number 90: ucscTrack3 (eval = FALSE)
###################################################
## cpgIslands <- UcscTrack(genome="mm9", chromosome="chrX", track="cpgIslandExt", from=from, to=to,
##                         trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",
##                         shape="box", fill="#006400", name="CpG Islands")
## 
## snpLocations <-  UcscTrack(genome="mm9", chromosome="chrX", track="snp128", from=from, to=to,
##                            trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",
##                            feature="func", strand="strand", shape="box", stacking="dense", fill="black",
##                            name="SNPs")


###################################################
### code chunk number 91: ucscTrack4 (eval = FALSE)
###################################################
## conservation <- UcscTrack(genome="mm9", chromosome="chrX", track="Conservation", table="phyloP30wayPlacental",
##                           from=from, to=to, trackType="DataTrack", start="start", end="end", data="score",
##                           type="hist", window="auto", col.histogram="darkblue", fill.histogram="darkblue", 
##                           ylim=c(-3.7, 4), name="Conservation")
## 
## gcContent <- UcscTrack(genome="mm9", chromosome="chrX", track="GC Percent", table="gc5Base",
##                        from=from, to=to, trackType="DataTrack", start="start", end="end", data="score",
##                        type="hist", window=-1, windowSize=1500, fill.histogram="black", col.histogram="black",
##                        ylim=c(30, 70), name="GC Percent")


###################################################
### code chunk number 92: ucscTrack5 (eval = FALSE)
###################################################
## axTrack <- GenomeAxisTrack()
## idxTrack <- IdeogramTrack(genome="mm9", chromosome="chrX")


###################################################
### code chunk number 93: ucscTrackLoad
###################################################
data(ucscItems)


###################################################
### code chunk number 94: ucscTrack6
###################################################
plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, ensGenes, cpgIslands,
                gcContent, conservation, snpLocations), from=from, to=to, showTitle=FALSE)


###################################################
### code chunk number 95: multPlot1
###################################################
chroms <- c("chr1", "chr2", "chr3", "chr4")
maTrack <- AnnotationTrack(range=GRanges(seqnames=chroms, ranges=IRanges(start=1, width=c(100,400,200,1000)),
                                         strand=c("+", "+", "-", "+")), genome="mm9", chromosome="chr1", name="foo")

mdTrack <- DataTrack(range=GRanges(seqnames=rep(chroms, c(10, 40, 20, 100)),
                                   ranges=IRanges(start=c(seq(1,100,len=10), seq(1,400,len=40), seq(1, 200, len=20), 
                                                          seq(1,1000, len=100)), width=9), values=runif(170)), 
                     data="values", chromosome="chr1", genome="mm9", name="bar")


###################################################
### code chunk number 96: multPlot2
###################################################
mgTrack <- GenomeAxisTrack(scale=50, labelPos="below", exponent=3)
chromosome(itrack) <- "chr1"


###################################################
### code chunk number 97: multPlot3
###################################################
ncols <- 2
nrows <- length(chroms)%/%ncols
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrows, ncols)))
for(i in seq_along(chroms)){
    pushViewport(viewport(layout.pos.col=((i-1)%%ncols)+1, layout.pos.row=(((i)-1)%/%ncols)+1))
    plotTracks(list(itrack, maTrack, mdTrack, mgTrack), chromosome=chroms[i], add=TRUE)
    popViewport(1)
}


###################################################
### code chunk number 98: multPlot4
###################################################
library(lattice)
chroms <- data.frame(chromosome=chroms)
xyplot(1~chromosome|chromosome, data=chroms, panel=function(x){plotTracks(list(itrack , maTrack, mdTrack, mgTrack), chromosome=x, add=TRUE, showId=FALSE)}, 
       scales=list(draw=FALSE), xlab=NULL, ylab=NULL)


###################################################
### code chunk number 99: session-info
###################################################
  sessionInfo()


