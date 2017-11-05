library(baySeq)

#TODO debug me
cargs <- c("input_4_bayseq","representative.id.gid");
cargs <- c("a","null","T",1,1)
cargs <- commandArgs(trailingOnly=T)

a              <- read.delim(cargs[1], header=T, sep="\t")
gid            <- NULL
nullData       <- T
bootStraps     <- 100
ncl            <- 1
nullLikelihood <- 0.8 #min likelihood to call not expressed; the higher the more permissive to be expressed
if (length(cargs) > 1 && !(cargs[2]=="null" || cargs[2]=="NULL")) gid <- read.table(cargs[2],header=F,col.names=c("ids","gid"))
if (length(cargs) > 2) nullData       <- as.logical(cargs[3])
if (length(cargs) > 3) bootStraps     <- as.numeric(cargs[4])
if (length(cargs) > 4) ncl            <- as.numeric(cargs[5])
if (length(cargs) > 5) nullLikelihood <- as.numeric(cargs[6])

print(c(nullData,bootStraps,ncl,nullLikelihood))
print(head(gid,1))

#filter out the 0 expression regions
a <- a[rowSums(a[grep("exp",names(a))])>0,];

#hard coded might be problematic
replicates <- c(1,1,2,2) #this corresponds to the conditions (i.e. ctrl1 ctrl2 trtmt1 trtmt2)
#groups
groups <- list(NDE = c(1,1,1,1), DE = replicates)

#ids vector
ids  <- a[,1]

#length matrix 
lens <- as.matrix(a[,grep("len",names(a))])

#read count matrix
all  <- a[,grep("exp",names(a))]

#dimensions of data
ngenes   <- dim(all)[1]
nsamples <- dim(all)[2]

#total expression
sums     <- colSums(all)

#remove tag from names
names(all) <- gsub(".exp","",names(all))

CD <- new("countData", data = as.matrix(all), replicates = replicates, libsizes = as.integer(sums), groups = groups, seglens = lens, annotation=as.data.frame(ids))
#CD@annotation <- as.data.frame(ids)

cl <- NULL
if ("snow" %in% installed.packages()[, 1]) { 
	library(snow)
	cl <- makeCluster(ncl, "SOCK") 
}

#priors
CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)

#plot prior for group NDE and DE(xn)
plot.new();pdf(paste(cargs[1],"_density.pdf",sep=''),h=8,w=8)
plotPriors(CDP.NBML, group = "NDE")
plotPriors(CDP.NBML, group = "DE")
dev.off()

#posteriors
CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = 'BIC', cl = cl, nullData = nullData, bootStraps = bootStraps)
print(CDPost.NBML@estProps)

if (!is.null(cl)) stopCluster(cl)

#all likelihoods
bayseq_de   <- topCounts(CDPost.NBML, group="DE", number=ngenes)
#get null model likelihoods and ids of expressed genes
if (nullData == T) {
	bayseq_null <- topCounts(CDPost.NBML, group=NULL, number=ngenes)
	expressed <- as.data.frame(bayseq_null[bayseq_null$Likelihood<nullLikelihood,"ids"])
	names(expressed) <- "ids"
	bayseq_de_exp <- merge(bayseq_de, expressed)
} else {
	bayseq_de_exp <- bayseq_de
}
head(bayseq_de_exp,2)
#----------------------------------------------------------------------------------
#output prep
#----------------------------------------------------------------------------------
#matrix of read counts of expressed genes
sums <- colSums(bayseq_de_exp[2:(nsamples+1)])
sums <- matrix(rep(sums,dim(bayseq_de_exp)[1]),byrow=T,ncol=nsamples)

#matrix of lengths of expressed genes with ids attached
lens <- a[1:2]; names(lens) <- c("ids","len")

#attached length
bayseq_de_exp <- merge(bayseq_de_exp, lens)
head(bayseq_de_exp)

#expressed gene length matrix
lens <- matrix(rep(bayseq_de_exp$len,nsamples),byrow=F,nrow=dim(bayseq_de_exp)[1])
head(lens)

#matrix of pernuc read counts of expressed genes
pernucExp <- bayseq_de_exp[2:(nsamples+1)] / sums / lens
names(pernucExp) <- sub("$",".pernucExp",names(pernucExp))

#attach pernucExp
bayseq_de_exp <- cbind(bayseq_de_exp,pernucExp)

#attach fold change from means of groups
#fold change between group means
e1 <- rowMeans(pernucExp[grep(levels(CD@groups$DE)[1],replicates)])
e2 <- rowMeans(pernucExp[grep(levels(CD@groups$DE)[2],replicates)])
bayseq_de_exp$log2FC <- log2(e1) - log2(e2) 

#attach gene names
if (!is.null(gid)) bayseq_de_exp <- merge(bayseq_de_exp,gid)

head(bayseq_de_exp,2)

#write
write.table(bayseq_de_exp[order(bayseq_de_exp$FDR),], file=paste(cargs[1],"_DE.txt",sep=""),sep="\t",quote=F,row.names=F)

