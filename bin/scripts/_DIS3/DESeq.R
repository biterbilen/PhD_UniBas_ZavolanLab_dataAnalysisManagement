###################################################
### chunk number 1: countsTable
###################################################
#line 62 "vignettes/DESeq/inst/doc/DESeq.Rnw"
library("DESeq")
countsTable <- read.delim(system.file("extra/TagSeqExample.tab", package="DESeq"), 
                          header=TRUE, stringsAsFactors=TRUE, row.names="gene")
conds <- factor(c("T", "T", "T", "Tb", "N", "N"))


###################################################
### chunk number 2: countsTable
###################################################
#line 69 "vignettes/DESeq/inst/doc/DESeq.Rnw"
head(countsTable)
conds


###################################################
### chunk number 3: demo1
###################################################
#line 75 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
res <- nbinomTest( cds, "T", "N")


###################################################
### chunk number 4: res
###################################################
#line 84 "vignettes/DESeq/inst/doc/DESeq.Rnw"
head(res)


###################################################
### chunk number 5: 
###################################################
#line 107 "vignettes/DESeq/inst/doc/DESeq.Rnw"
library( DESeq )
exampleFile = system.file( "extra/TagSeqExample.tab", package="DESeq" )
exampleFile


###################################################
### chunk number 6: 
###################################################
#line 114 "vignettes/DESeq/inst/doc/DESeq.Rnw"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE )
head( countsTable )


###################################################
### chunk number 7: 
###################################################
#line 122 "vignettes/DESeq/inst/doc/DESeq.Rnw"
rownames( countsTable ) <- countsTable$gene
countsTable <- countsTable[ , -1 ]


###################################################
### chunk number 8: 
###################################################
#line 132 "vignettes/DESeq/inst/doc/DESeq.Rnw"
conds <- c( "T", "T", "T", "Tb", "N", "N" )


###################################################
### chunk number 9: 
###################################################
#line 143 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds <- newCountDataSet( countsTable, conds )


###################################################
### chunk number 10: 
###################################################
#line 150 "vignettes/DESeq/inst/doc/DESeq.Rnw"
head( counts(cds) )


###################################################
### chunk number 11: 
###################################################
#line 156 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds <- cds[ ,-1 ]


###################################################
### chunk number 12: 
###################################################
#line 168 "vignettes/DESeq/inst/doc/DESeq.Rnw"
libsizes <- c( T1a=6843583, T1b=7604834, T2=13625570, T3=12291910,
   N1=12872125, N2=10502656 )
sizeFactors(cds) <- libsizes[-1]


###################################################
### chunk number 13: 
###################################################
#line 176 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds <- estimateSizeFactors( cds )
sizeFactors( cds )


###################################################
### chunk number 14: 
###################################################
#line 190 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds <- estimateVarianceFunctions( cds )


###################################################
### chunk number 15: 
###################################################
#line 212 "vignettes/DESeq/inst/doc/DESeq.Rnw"
countValue <- 123
baseLevel <- countValue / sizeFactors(cds)["T1b"]
rawVarFuncForGB <- rawVarFunc( cds, "T" )
rawVariance <- rawVarFuncForGB( baseLevel )
fullVariance <- countValue + rawVariance * sizeFactors(cds)["T1b"]^2
sqrt( fullVariance )


###################################################
### chunk number 16: figSCV
###################################################
#line 239 "vignettes/DESeq/inst/doc/DESeq.Rnw"
scvPlot( cds, ylim=c(0,2) )


###################################################
### chunk number 17: 
###################################################
#line 265 "vignettes/DESeq/inst/doc/DESeq.Rnw"
rawVarFuncTable( cds )


###################################################
### chunk number 18: 
###################################################
#line 287 "vignettes/DESeq/inst/doc/DESeq.Rnw"
diagForT <- varianceFitDiagnostics( cds, "T" )
head( diagForT )


###################################################
### chunk number 19: figFit
###################################################
#line 301 "vignettes/DESeq/inst/doc/DESeq.Rnw"
smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )
lines( log10(fittedBaseVar) ~ log10(baseMean), 
   diagForT[ order(diagForT$baseMean), ], col="red" )


###################################################
### chunk number 20: figECDF
###################################################
#line 325 "vignettes/DESeq/inst/doc/DESeq.Rnw"
par( mfrow=c(1,2 ) )
residualsEcdfPlot( cds, "T" )
residualsEcdfPlot( cds, "N" )


###################################################
### chunk number 21: 
###################################################
#line 354 "vignettes/DESeq/inst/doc/DESeq.Rnw"
res <- nbinomTest( cds, "N", "T" )
head(res)


###################################################
### chunk number 22: figDE
###################################################
#line 377 "vignettes/DESeq/inst/doc/DESeq.Rnw"
plotDE <- function( res )
   plot( 
      res$baseMean, 
      res$log2FoldChange, 
      log="x", pch=20, cex=.1, 
      col = ifelse( res$padj < .1, "red", "black" ) )

plotDE( res )


###################################################
### chunk number 23: 
###################################################
#line 391 "vignettes/DESeq/inst/doc/DESeq.Rnw"
resSig <- res[ res$padj < .1, ]


###################################################
### chunk number 24: 
###################################################
#line 395 "vignettes/DESeq/inst/doc/DESeq.Rnw"
head( resSig[ order(resSig$pval), ] )


###################################################
### chunk number 25: 
###################################################
#line 401 "vignettes/DESeq/inst/doc/DESeq.Rnw"
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )


###################################################
### chunk number 26: 
###################################################
#line 405 "vignettes/DESeq/inst/doc/DESeq.Rnw"
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )


###################################################
### chunk number 27: figResVarDens
###################################################
#line 427 "vignettes/DESeq/inst/doc/DESeq.Rnw"
plot( density( res$resVarA, na.rm=TRUE, from=0, to=20 ), col="red" )
lines( density( res$resVarB, na.rm=TRUE, from=0, to=20 ), col="blue" )
xg <- seq( 0, 20, length.out=1000 ); lines( xg, dchisq( xg, df=1 ), col="grey" )


###################################################
### chunk number 28: 
###################################################
#line 444 "vignettes/DESeq/inst/doc/DESeq.Rnw"
table( res$resVarA > 15 | res$resVarB > 15)


###################################################
### chunk number 29: 
###################################################
#line 449 "vignettes/DESeq/inst/doc/DESeq.Rnw"
( 1 - pchisq( 15, df=1) ) * nrow(counts(cds))


###################################################
### chunk number 30: 
###################################################
#line 517 "vignettes/DESeq/inst/doc/DESeq.Rnw"
resTbvsN <- nbinomTest( cds, "N", "Tb" )


###################################################
### chunk number 31: figDE_Tb
###################################################
#line 521 "vignettes/DESeq/inst/doc/DESeq.Rnw"
plot( 
   resTbvsN$baseMean, 
   resTbvsN$log2FoldChange, 
   log="x", pch=20, cex=.1, 
   col = ifelse( resTbvsN$padj < .1, "red", "black" ) )


###################################################
### chunk number 32: 
###################################################
#line 559 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds2 <- cds[ ,c( "T1b", "N1" ) ]


###################################################
### chunk number 33: 
###################################################
#line 566 "vignettes/DESeq/inst/doc/DESeq.Rnw"
cds2 <- estimateVarianceFunctions( cds2, method="blind" )


###################################################
### chunk number 34: 
###################################################
#line 571 "vignettes/DESeq/inst/doc/DESeq.Rnw"
res2 <- nbinomTest( cds2, "N", "T" )


###################################################
### chunk number 35: figDE2
###################################################
#line 583 "vignettes/DESeq/inst/doc/DESeq.Rnw"
plot( 
   res2$baseMean,
   res2$log2FoldChange, 
   log="x", pch=20, cex=.1, 
   col = ifelse( res2$padj < .1, "red", "black" ) )


###################################################
### chunk number 36: 
###################################################
#line 591 "vignettes/DESeq/inst/doc/DESeq.Rnw"
addmargins( table( res_sig = res$padj < .1, res2_sig = res2$padj < .1 ) )


###################################################
### chunk number 37: 
###################################################
#line 604 "vignettes/DESeq/inst/doc/DESeq.Rnw"
colsN <- conditions(cds) == "N"
colsT <- conditions(cds) == "T"
baseMeansNT <- getBaseMeansAndVariances( 
   counts(cds)[ , colsN|colsT ], 
   sizeFactors(cds)[ colsN|colsT ] )$baseMean
pvals2b <- nbinomTestForMatrices(
   counts(cds)[ ,colsN ],
   counts(cds)[ ,colsT ],
   sizeFactors(cds)[ colsN ],
   sizeFactors(cds)[ colsT ],
   rawVarFunc( cds2, "_blind", TRUE )( baseMeansNT ),
   rawVarFunc( cds2, "_blind", TRUE )( baseMeansNT ) )


###################################################
### chunk number 38: padj
###################################################
#line 620 "vignettes/DESeq/inst/doc/DESeq.Rnw"
padj2b <- p.adjust( pvals2b, method="BH" )
notNAinRes2 <- !is.na( res2$padj )
addmargins( table( 
   res_sig = res$padj[notNAinRes2] < .1, 
   res2b_sig = padj2b[notNAinRes2] < .1 ) )


###################################################
### chunk number 39: vsd
###################################################
#line 672 "vignettes/DESeq/inst/doc/DESeq.Rnw"
vsd <- getVarianceStabilizedData( cds )


###################################################
### chunk number 40: modlr
###################################################
#line 682 "vignettes/DESeq/inst/doc/DESeq.Rnw"
mod_lfc = (rowMeans( vsd[, conditions(cds)=="T"] ) - 
           rowMeans( vsd[, conditions(cds)=="N"] ))


###################################################
### chunk number 41: dah
###################################################
#line 693 "vignettes/DESeq/inst/doc/DESeq.Rnw"
lfc = res$log2FoldChange
finite = is.finite(lfc)
table(as.character(lfc[!finite]), useNA="always")


###################################################
### chunk number 42: repl
###################################################
#line 700 "vignettes/DESeq/inst/doc/DESeq.Rnw"
LargeNumber = 10
lfc = ifelse(finite, lfc, sign(lfc)*LargeNumber)


###################################################
### chunk number 43: figscp1 eval=FALSE
###################################################
## #line 705 "vignettes/DESeq/inst/doc/DESeq.Rnw"
## plot( lfc, mod_lfc, pch=16, 
##       col = ifelse(finite, "#80808040", "red"))
## abline(a=0, b=1, col="#40404040")


###################################################
### chunk number 44: figscp2
###################################################
#line 710 "vignettes/DESeq/inst/doc/DESeq.Rnw"
png(file="DESeq-figmodlr.png", width=800, height=800, pointsize=24)
plot( lfc, mod_lfc, pch=16, 
      col = ifelse(finite, "#80808040", "red"))
abline(a=0, b=1, col="#40404040")
dev.off()


###################################################
### chunk number 45: dists
###################################################
#line 742 "vignettes/DESeq/inst/doc/DESeq.Rnw"
dists <- dist( t( vsd ) )


###################################################
### chunk number 46: figHeatmap1
###################################################
#line 747 "vignettes/DESeq/inst/doc/DESeq.Rnw"
heatmap(as.matrix( dists ), 
        symm=TRUE,
        scale="none",
        col = colorRampPalette(c("darkblue","white"))(100))


###################################################
### chunk number 47: figHeatmap2a
###################################################
#line 762 "vignettes/DESeq/inst/doc/DESeq.Rnw"
select = order(res$pval)[1:100]
colors = colorRampPalette(c("white","darkblue"))(100)
heatmap( vsd[select,], 
         col = colors, scale = "none")


###################################################
### chunk number 48: figHeatmap2b
###################################################
#line 771 "vignettes/DESeq/inst/doc/DESeq.Rnw"
heatmap( counts(cds)[select,], 
         col = colors, scale = "none")


###################################################
### chunk number 49: sessi
###################################################
#line 845 "vignettes/DESeq/inst/doc/DESeq.Rnw"
sessionInfo()


