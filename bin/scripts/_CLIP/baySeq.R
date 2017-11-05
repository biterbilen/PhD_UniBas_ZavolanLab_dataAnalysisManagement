### R code from vignette source 'vignettes/baySeq/inst/doc/baySeq.Rnw'

###################################################
### code chunk number 1: baySeq.Rnw:31-33
###################################################
set.seed(102)
options(width = 90)


###################################################
### code chunk number 2: baySeq.Rnw:37-38
###################################################
library(baySeq)


###################################################
### code chunk number 3: baySeq.Rnw:42-44 (eval = FALSE)
###################################################
library(snow)
cl <- makeCluster(4, "SOCK")


###################################################
### code chunk number 4: baySeq.Rnw:48-49
###################################################
##  cl <- NULL



###################################################
### code chunk number 5: baySeq.Rnw:57-59
###################################################
data(simData)
simData[1:10,]


###################################################
### code chunk number 6: baySeq.Rnw:63-64
###################################################
  replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")


###################################################
### code chunk number 7: baySeq.Rnw:73-75
###################################################
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1),
               DE = c(1,1,1,1,1,2,2,2,2,2))


###################################################
### code chunk number 8: baySeq.Rnw:86-87
###################################################
CD <- new("countData", data = simData, replicates = replicates, groups = groups)


###################################################
### code chunk number 9: baySeq.Rnw:92-93
###################################################
  CD@libsizes <- getLibsizes(CD)


###################################################
### code chunk number 10: plotMA
###################################################
plotMA.CD(CD, samplesA = "simA", samplesB = "simB", col = c(rep("red", 100), rep("black", 900)))


###################################################
### code chunk number 11: figPlotMA
###################################################
plotMA.CD(CD, samplesA = "simA", samplesB = "simB", col = c(rep("red", 100), rep("black", 900)))


###################################################
### code chunk number 12: baySeq.Rnw:115-116
###################################################
CD@annotation <- data.frame(name = paste("count", 1:1000, sep = "_"))


###################################################
### code chunk number 13: baySeq.Rnw:122-123
###################################################
CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)


###################################################
### code chunk number 14: baySeq.Rnw:129-133
###################################################
CD <- getLikelihoods.NB(CD, pET = 'BIC', cl = cl)
CD@estProps
CD@posteriors[1:10,]
CD@posteriors[101:110,]


###################################################
### code chunk number 15: baySeq.Rnw:138-139
###################################################
CD@estProps[2]


###################################################
### code chunk number 16: baySeq.Rnw:146-147
###################################################
topCounts(CD, group = "DE")  


###################################################
### code chunk number 17: plotPosteriors
###################################################
plotPosteriors(CD, group = "DE", col = c(rep("red", 100), rep("black", 900)))


###################################################
### code chunk number 18: figPlotPosteriors
###################################################
plotPosteriors(CD, group = "DE", col = c(rep("red", 100), rep("black", 900)))


###################################################
### code chunk number 19: baySeq.Rnw:168-170
###################################################
if(!is.null(cl))
  stopCluster(cl)


###################################################
### code chunk number 20: baySeq.Rnw:182-184
###################################################
data(mobData)
data(mobAnnotation)


###################################################
### code chunk number 21: baySeq.Rnw:195-197
###################################################
seglens <- mobAnnotation$end - mobAnnotation$start + 1
cD <- new("countData", data = mobData, seglens = seglens, annotation = mobAnnotation)


###################################################
### code chunk number 22: baySeq.Rnw:201-202
###################################################
cD@libsizes <- getLibsizes(cD, estimationType = "quantile")


###################################################
### code chunk number 23: baySeq.Rnw:210-211
###################################################
cDPair <- cD[,1:4]


###################################################
### code chunk number 24: baySeq.Rnw:215-216
###################################################
replicates(cDPair) <- as.factor(c("D3/D3", "D3/D3", "WT/D3", "WT/D3"))


###################################################
### code chunk number 25: baySeq.Rnw:222-223
###################################################
NDE <- c(1,1,1,1)


###################################################
### code chunk number 26: baySeq.Rnw:228-229
###################################################
mobile <- c("non-mobile","non-mobile","mobile","mobile")


###################################################
### code chunk number 27: baySeq.Rnw:234-235
###################################################
groups(cDPair) <- list(NDE = NDE, mobile = mobile)


###################################################
### code chunk number 28: baySeq.Rnw:239-240
###################################################
cDPair <- getPriors.NB(cDPair, samplesize = 1e4, cl = NULL)


###################################################
### code chunk number 29: plotPriors
###################################################
plotPriors(cDPair, group = "NDE")


###################################################
### code chunk number 30: figPlotPriors
###################################################
plotPriors(cDPair, group = "NDE")


###################################################
### code chunk number 31: baySeq.Rnw:264-265
###################################################
cDPair <- getLikelihoods.NB(cDPair, nullData = TRUE, cl = NULL)


###################################################
### code chunk number 32: baySeq.Rnw:269-270
###################################################
cDPair


###################################################
### code chunk number 33: baySeq.Rnw:275-276
###################################################
cDPair@estProps


###################################################
### code chunk number 34: baySeq.Rnw:283-284
###################################################
topCounts(cDPair, group = 2)


###################################################
### code chunk number 35: baySeq.Rnw:289-290
###################################################
topCounts(cDPair, group = 2, normaliseData = TRUE)


###################################################
### code chunk number 36: baySeq.Rnw:295-296
###################################################
topCounts(cDPair, group = NULL, number = 500)


###################################################
### code chunk number 37: plotPairPosteriors
###################################################
plotPosteriors(cDPair, group = 2, samplesA = 1:2, samplesB = 3:4)


###################################################
### code chunk number 38: figPlotPairPosteriors
###################################################
plotPosteriors(cDPair, group = 2, samplesA = 1:2, samplesB = 3:4)


###################################################
### code chunk number 39: plotMAPost
###################################################
plotMA.CD(cDPair, samplesA = c(1,2), samplesB = c(3,4),
          col = rgb(red = exp(cDPair@posteriors[,2]), green = 0, blue = 0))


###################################################
### code chunk number 40: figPlotMAPost
###################################################
plotMA.CD(cDPair, samplesA = c(1,2), samplesB = c(3,4),
          col = rgb(red = exp(cDPair@posteriors[,2]), green = 0, blue = 0))


###################################################
### code chunk number 41: baySeq.Rnw:335-336
###################################################
cD@replicates <- as.factor(c("D3/D3", "D3/D3", "WT/D3", "WT/D3", "WT/WT", "WT/WT"))


###################################################
### code chunk number 42: baySeq.Rnw:342-343
###################################################
NDE <- factor(c(1,1,1,1,1,1))


###################################################
### code chunk number 43: baySeq.Rnw:348-349
###################################################
d3dep <- c("wtRoot","wtRoot","wtRoot","wtRoot","dicerRoot","dicerRoot")


###################################################
### code chunk number 44: baySeq.Rnw:354-355
###################################################
mobile <- c("dicerShoot","dicerShoot","wtShoot","wtShoot","wtShoot","wtShoot")


###################################################
### code chunk number 45: baySeq.Rnw:360-361
###################################################
groups(cD) <- list(NDE = NDE, d3dep = d3dep, mobile = mobile)  


###################################################
### code chunk number 46: baySeq.Rnw:367-369
###################################################
cD <- getPriors.NB(cD, cl = NULL)
cD <- getLikelihoods.NB(cD, nullData = TRUE, cl = NULL)


###################################################
### code chunk number 47: baySeq.Rnw:373-374
###################################################
topCounts(cD, group = "mobile", normaliseData = TRUE)  


###################################################
### code chunk number 48: baySeq.Rnw:378-379
###################################################
topCounts(cD, group = "d3dep", normaliseData = TRUE)  


