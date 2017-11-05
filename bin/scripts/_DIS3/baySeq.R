###################################################
### chunk number 1: 
###################################################
#line 32 "vignettes/baySeq/inst/doc/baySeq.Rnw"
set.seed(102)


###################################################
### chunk number 2: 
###################################################
#line 37 "vignettes/baySeq/inst/doc/baySeq.Rnw"
library(baySeq)


###################################################
### chunk number 3:  eval=FALSE
###################################################
## #line 42 "vignettes/baySeq/inst/doc/baySeq.Rnw"
##   library(snow)
## cl <- makeCluster(4, "SOCK")


###################################################
### chunk number 4: 
###################################################
#line 48 "vignettes/baySeq/inst/doc/baySeq.Rnw"
  cl <- NULL


###################################################
### chunk number 5: 
###################################################
#line 57 "vignettes/baySeq/inst/doc/baySeq.Rnw"
data(simCount)
data(libsizes)
simCount[1:10,]
libsizes


###################################################
### chunk number 6: 
###################################################
#line 65 "vignettes/baySeq/inst/doc/baySeq.Rnw"
  replicates <- c(1,1,1,1,1,2,2,2,2,2)


###################################################
### chunk number 7: 
###################################################
#line 75 "vignettes/baySeq/inst/doc/baySeq.Rnw"
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1),
               DE = c(1,1,1,1,1,2,2,2,2,2))


###################################################
### chunk number 8: 
###################################################
#line 88 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CD <- new("countData", data = simCount, replicates = replicates, libsizes = as.integer(libsizes), groups = groups)


###################################################
### chunk number 9: plotMA
###################################################
#line 94 "vignettes/baySeq/inst/doc/baySeq.Rnw"
plotMA.CD(CD, samplesA = 1:5, samplesB = 6:10, col = c(rep("red", 100), rep("black", 900)))


###################################################
### chunk number 10: figPlotMA
###################################################
#line 100 "vignettes/baySeq/inst/doc/baySeq.Rnw"
#line 94 "vignettes/baySeq/inst/doc/baySeq.Rnw#from line#100#"
plotMA.CD(CD, samplesA = 1:5, samplesB = 6:10, col = c(rep("red", 100), rep("black", 900)))
#line 101 "vignettes/baySeq/inst/doc/baySeq.Rnw"


###################################################
### chunk number 11: 
###################################################
#line 111 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CD@annotation <- data.frame(name = paste("count", 1:1000, sep = "_"))


###################################################
### chunk number 12: 
###################################################
#line 131 "vignettes/baySeq/inst/doc/baySeq.Rnw"
  CDP.Poi <- getPriors.Pois(CD, samplesize = 20, takemean = TRUE, cl = cl)


###################################################
### chunk number 13: 
###################################################
#line 137 "vignettes/baySeq/inst/doc/baySeq.Rnw"
  CDP.Poi@priors


###################################################
### chunk number 14: 
###################################################
#line 145 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CDPost.Poi <- getLikelihoods.Pois(CDP.Poi, pET = "BIC", cl = cl)
CDPost.Poi@estProps
CDPost.Poi@posteriors[1:10,]
CDPost.Poi@posteriors[101:110,]


###################################################
### chunk number 15: 
###################################################
#line 155 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CDPost.Poi@estProps[2]


###################################################
### chunk number 16: 
###################################################
#line 163 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CDP.NBML <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)


###################################################
### chunk number 17: 
###################################################
#line 170 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = 'BIC', cl = cl)
CDPost.NBML@estProps
CDPost.NBML@posteriors[1:10,]
CDPost.NBML@posteriors[101:110,]


###################################################
### chunk number 18: 
###################################################
#line 179 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CDPost.NBML@estProps[2]


###################################################
### chunk number 19: 
###################################################
#line 187 "vignettes/baySeq/inst/doc/baySeq.Rnw"
topCounts(CDPost.NBML, group = 2)  


###################################################
### chunk number 20: 
###################################################
#line 193 "vignettes/baySeq/inst/doc/baySeq.Rnw"
NBML.TPs <- getTPs(CDPost.NBML, group = 2, TPs = 1:100)
Poi.TPs <- getTPs(CDPost.Poi, group = 2, TPs = 1:100)


###################################################
### chunk number 21: plotPosteriors
###################################################
#line 199 "vignettes/baySeq/inst/doc/baySeq.Rnw"
plotPosteriors(CDPost.NBML, group = 2, samplesA = 1:5, samplesB = 6:10, col = c(rep("red", 100), rep("black", 900)))


###################################################
### chunk number 22: figPlotPosteriors
###################################################
#line 206 "vignettes/baySeq/inst/doc/baySeq.Rnw"
#line 199 "vignettes/baySeq/inst/doc/baySeq.Rnw#from line#206#"
plotPosteriors(CDPost.NBML, group = 2, samplesA = 1:5, samplesB = 6:10, col = c(rep("red", 100), rep("black", 900)))
#line 207 "vignettes/baySeq/inst/doc/baySeq.Rnw"


###################################################
### chunk number 23: FPsPlot
###################################################
#line 215 "vignettes/baySeq/inst/doc/baySeq.Rnw"
plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(0, log(1000)), xlab = "Number of counts selected", ylab = "(Log) False Positives")
lines(x = 1:1000, y = log(1:1000 - NBML.TPs[1:1000]), type = "l", col = "red")
lines(x = 1:1000, y = log(1:1000 - Poi.TPs[1:1000]), type = "l", col = "blue")

legend(x = "topleft", lty = c(1,1), col = c("red", "blue"), legend = c("Negative Binomial", "Poisson-Gamma"))


###################################################
### chunk number 24: figFPs
###################################################
#line 226 "vignettes/baySeq/inst/doc/baySeq.Rnw"
#line 215 "vignettes/baySeq/inst/doc/baySeq.Rnw#from line#226#"
plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(0, log(1000)), xlab = "Number of counts selected", ylab = "(Log) False Positives")
lines(x = 1:1000, y = log(1:1000 - NBML.TPs[1:1000]), type = "l", col = "red")
lines(x = 1:1000, y = log(1:1000 - Poi.TPs[1:1000]), type = "l", col = "blue")

legend(x = "topleft", lty = c(1,1), col = c("red", "blue"), legend = c("Negative Binomial", "Poisson-Gamma"))
#line 227 "vignettes/baySeq/inst/doc/baySeq.Rnw"


###################################################
### chunk number 25: 
###################################################
#line 237 "vignettes/baySeq/inst/doc/baySeq.Rnw"
if(!is.null(cl))
  stopCluster(cl)


###################################################
### chunk number 26: 
###################################################
#line 259 "vignettes/baySeq/inst/doc/baySeq.Rnw"
set.seed(101)


###################################################
### chunk number 27: 
###################################################
#line 265 "vignettes/baySeq/inst/doc/baySeq.Rnw"
data(factData)
data(factlibsizes)


###################################################
### chunk number 28: 
###################################################
#line 272 "vignettes/baySeq/inst/doc/baySeq.Rnw"
  replicates <- c(1,1,2,2,3,3,4,4)
factgroups <- list(NDE = c(1,1,1,1,1,1,1,1),
                   DE.A.B = c(1,1,1,1,2,2,2,2),
                   DE.C.D = c(1,1,2,2,1,1,2,2))


###################################################
### chunk number 29: 
###################################################
#line 285 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CDfact <- new("countData", data = factCount, replicates = replicates, libsizes = factlibsizes, groups = factgroups)
CDfact@annotation <- data.frame(name = paste("count", 1:1000, sep = "_"))
CDfactP.NBML <- getPriors.NB(CDfact, samplesize = 1000, estimation = "QL", cl = cl)
CDfactPost.NBML <- getLikelihoods.NB(CDfactP.NBML, pET = "BIC", cl = cl)
CDfactPost.NBML@estProps


###################################################
### chunk number 30: 
###################################################
#line 294 "vignettes/baySeq/inst/doc/baySeq.Rnw"
topCounts(CDfactPost.NBML, group = 2)  


###################################################
### chunk number 31: 
###################################################
#line 299 "vignettes/baySeq/inst/doc/baySeq.Rnw"
topCounts(CDfactPost.NBML, group = 3)  


###################################################
### chunk number 32: 
###################################################
#line 311 "vignettes/baySeq/inst/doc/baySeq.Rnw"
data(simSeg)
data(libsizes)
simSeg[1:10,]
libsizes


###################################################
### chunk number 33: 
###################################################
#line 320 "vignettes/baySeq/inst/doc/baySeq.Rnw"
replicates <- c(1,1,1,1,1,2,2,2,2,2)
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1),
               DE = c(1,1,1,1,1,2,2,2,2,2))


###################################################
### chunk number 34: 
###################################################
#line 328 "vignettes/baySeq/inst/doc/baySeq.Rnw"
SD <- new("countData", data = simSeg[,-1], replicates = replicates, libsizes = libsizes, groups = groups, seglens = simSeg[,1])


###################################################
### chunk number 35: 
###################################################
#line 334 "vignettes/baySeq/inst/doc/baySeq.Rnw"
SD@annotation <- data.frame(name = paste("gene", 1:1000, sep = "_"))


###################################################
### chunk number 36: 
###################################################
#line 340 "vignettes/baySeq/inst/doc/baySeq.Rnw"
SDP.NBML <- getPriors.NB(SD, samplesize = 1000, estimation = "QL", cl = cl)
SDP.Pois <- getPriors.Pois(SD, samplesize = 20, cl = cl)


###################################################
### chunk number 37: 
###################################################
#line 346 "vignettes/baySeq/inst/doc/baySeq.Rnw"
SDPost.Pois <- getLikelihoods.Pois(SDP.Pois, pET = "BIC", cl = cl)  
SDPost.NBML <- getLikelihoods.NB(SDP.NBML, pET = "BIC", cl = cl)


###################################################
### chunk number 38: 
###################################################
#line 353 "vignettes/baySeq/inst/doc/baySeq.Rnw"
CSD <- new("countData", data = simSeg[,-1], replicates = replicates, libsizes = libsizes, groups = groups)
CSD@annotation <- data.frame(name = paste("gene", 1:1000, sep = "_"))

CSDP.NBML <- getPriors.NB(CSD, samplesize = 1000, estimation = "QL", cl = cl)
CSDPost.NBML <- getLikelihoods.NB(CSDP.NBML, pET = "BIC", cl = cl)
CSDP.Pois <- getPriors.Pois(CSD, samplesize = 20, cl = cl)
CSDPost.Pois <- getLikelihoods.Pois(CSDP.Pois, pET = "BIC", cl = cl)


###################################################
### chunk number 39: FPseglen
###################################################
#line 367 "vignettes/baySeq/inst/doc/baySeq.Rnw"
SD.NBML.FPs <- 1:nrow(simSeg) - getTPs(SDPost.NBML, group = 2, TPs = 1:100)
CSD.NBML.FPs <- 1:nrow(simSeg) - getTPs(CSDPost.NBML, group = 2, TPs = 1:100)
SD.Pois.FPs <- 1:nrow(simSeg) - getTPs(SDPost.Pois, group = 2, TPs = 1:100)
CSD.Pois.FPs <- 1:nrow(simSeg) - getTPs(CSDPost.Pois, group = 2, TPs = 1:100)

plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(0, log(100)), xlab = "Number of counts selected", ylab = "(Log) False Positives")
lines(x = 1:1000, y = log(CSD.NBML.FPs[1:1000]), type = "l", col = "red")
lines(x = 1:1000, y = log(SD.NBML.FPs[1:1000]), type = "l", col = "red", lty = 2)
lines(x = 1:1000, y = log(CSD.Pois.FPs[1:1000]), type = "l", col = "blue")
lines(x = 1:1000, y = log(SD.Pois.FPs[1:1000]), type = "l", col = "blue", lty = 2)
legend(x = "topleft", lty = c(1,2,1,2), col = c("red", "red", "blue", "blue"), legend = c("Negative Binomial (ignoring segment lengths)", "Negative Binomial (including segment lengths)", "Poisson-Gamma (ignoring segment lengths)", "Poisson-Gamma (including segment lengths)"))
                                                                               


###################################################
### chunk number 40: 
###################################################
#line 390 "vignettes/baySeq/inst/doc/baySeq.Rnw"
NSDPost.NBML <- getLikelihoods.NB(SDP.NBML, pET = "BIC", nullData = TRUE, bootStraps = 1, cl = cl)
NCSDPost.NBML <- getLikelihoods.NB(CSDP.NBML, pET = "BIC", nullData = TRUE, bootStraps = 1, cl = cl)


###################################################
### chunk number 41: 
###################################################
#line 397 "vignettes/baySeq/inst/doc/baySeq.Rnw"
topCounts(NSDPost.NBML, group = NULL)
topCounts(NCSDPost.NBML, group = NULL)


###################################################
### chunk number 42: FPNull
###################################################
#line 404 "vignettes/baySeq/inst/doc/baySeq.Rnw"
NSD.FPs <- 1:nrow(simSeg) - getTPs(NSDPost.NBML, group = NULL, TPs = 101:200, decreasing = TRUE)
NCSD.FPs <- 1:nrow(simSeg) - getTPs(NCSDPost.NBML, group = NULL, TPs = 101:200, decreasing = TRUE)
plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(0, log(1000)), xlab = "Number of counts selected", ylab = "(Log) False Positives")
lines(x = 1:1000, y = log(NCSD.FPs[1:1000]), type = "l", col = "red")
lines(x = 1:1000, y = log(NSD.FPs[1:1000]), type = "l", lty = 2, col = "red")
legend(x = "topleft", lty = c(1,2), col = c("red", "red"), legend = c("Negative Binomial (ignoring segment lengths)", "Negative Binomial (including segment lengths)"))


###################################################
### chunk number 43: pltPriors
###################################################
#line 419 "vignettes/baySeq/inst/doc/baySeq.Rnw"
plotPriors(SDP.NBML, group = 1)


###################################################
### chunk number 44: figPriors
###################################################
#line 426 "vignettes/baySeq/inst/doc/baySeq.Rnw"
#line 419 "vignettes/baySeq/inst/doc/baySeq.Rnw#from line#426#"
plotPriors(SDP.NBML, group = 1)
#line 427 "vignettes/baySeq/inst/doc/baySeq.Rnw"


