###################################################
### chunk number 1: load <- library
###################################################
#line 39 "vignettes/goseq/inst/doc/goseq.Rnw"
library(goseq)


###################################################
### chunk number 2: set <- width
###################################################
#seq.Rnw"
## gene.vector=as.integer(assayed.genes%in%de.genes)
## names(gene.vector)=assayed.genes
## head(gene.vector)


###################################################
### chunk number 4: supported <- genomes eval=FALSE
###########################
## #line 86 "vignettes/goseq/inst/doc/goseq.Rnw"
## supportedGeneIDs()


###################################################
### chunk number 6: getLengthDataFromUCSC eval=FALSE
###################################################
## #li######
#line 140 "vignettes/goseq/inst/doc/goseq.Rnw"
library(edgeR)
table.summary=read.table(system.file("extdata","Li_sum.txt",package='goseq'),sep='\t',header=TRUE,stringsAsFactors=FALSE)
counts=table.summary[,-1]
rownames(counts)=table.summary[,1]
grp=factor(rep(c("Control","Treated"),times=c(4,3)))
summarized=DGEList(counts,lib.size=colSums(counts),group=grp)


#########)
topTags(tested)


###################################################
### chunk number 9: edger <- 3
###################################################
#line 160 "vignettes/goseq/inst/doc/goseq.Rnw"
genes=as.integer(p.adjust(tested$table$p.value[tested$table$logFC!=0],method="BH")<.05)
names(genes)=row.names(tested$table[tested$table$logFC!=0,])
table(genes)


################# chunk number 11: head <- geneids
###################################################
#line 176 "vignettes/goseq/inst/doc/goseq.Rnw"
head(supportedGeneIDs(),n=12)


###################################################
### chunk number 12: pwf
###################################################
#line 188 "vignettes/g
GO.wall=goseq(pwf,"hg18","ensGene")
head(GO.wall)


###################################################
### chunk number 14: GO.samp
###################################################
#line 215 "vignettes/goseq/inst/doc/goseq.Rnw"
GO.samp=goseq(pwf,"hg18","ensGene",method="Sampling",repcnt=1000)


###################################################
### chunk number 15: head <- samp
###################################################
#line 218 "vignettes/goseq/inst/doc/goseq.Rnw"
head(GO.samp)


###################################################
### chunk number 16: plot <- wal <- v <- samp
##################################################################################
### chunk number 17: GO.nobias
###################################################
#line 235 "vignettes/goseq/inst/doc/goseq.Rnw"
GO.nobias=goseq(pwf,"hg18","ensGene",method="Hypergeometric")
head(GO.nobias)


#################################
#line 262 "vignettes/goseq/inst/doc/goseq.Rnw"
enriched.GO=GO.wall$category[p.adjust(GO.wall$over <- represented <- pvalue,method="BH")<.05]
head(enriched.GO)


###################################################
### chunk number 21: GO <- explained
###################################################
#line 269 "vignettes/goseq/inst/doc/goseq.Rnw"
library(GO.db)
for(go in enriched.GO[1:10]){
		print(GOTERM[[go]])
			cgo(names(genes),"hg18","ensGene")
			length(go)
			length(genes)
			head(go)


			###################################################
			### sis eval=FALSE
			###################################################
			## #line 307 "vignettes/goseq/inst/doc/goseq.Rnw"
			## pwf=nullp(genes,"hg18","ensGene")
			## go=goseq(pwf,"hg18","ensGene")


			###################################################
			### chunk number 26: verbose <- analysis eval=FALSE
			###################################################
			## #line 312 "vignettes/goseq/inst/doc/goseq.Rnw"
			## gene <- lengths=getlength(names(genes),"hg18","ensGene")
			## pwf=nullp(genes,bias.data=gene <- lengths)
			## go <- map=getgo(names(genes),"hg18","ensGene")
			## go=goseq(pwf,"hg18","ensENSEMBL 2 Entrez
			## en2eg=as.list(org.Hs.egENSEMBL2EG)
			## #Get the mapping from Entrez 2 KEGG
			## eg2kegg=as.list(org.Hs.egPATHmbine the two maps
			## kegg=lapply(en2eg,grepKEGG,eg2kegg)
			## hene')
			KEGG=goseq(pwf,'hg18','ensGene',test.cats="KEGG")
			head(KEGG)


			###################################################
			### chunk number 30: KEGG <- from <- db
			###################################################
			#line 361 "vignettes/goseq/inst/doc/goseq.Rnw"
			kegg=as.list(org.Hs.egPATH)
			head(kegg)


			###################################################
			### chunk number 31: countbias
			###################################################
			#line 376 "v##################################
			#line 395 "vignettes/goseq/inst/doc/goseq.Rnw"
			sessionInfo()



