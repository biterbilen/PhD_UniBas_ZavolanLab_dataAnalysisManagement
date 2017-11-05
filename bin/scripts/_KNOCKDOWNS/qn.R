qn <- function(file, ofile,wheader) {
#D.HuR=read.delim("~bilebi00/_KNOCKDOWNS/data/MRNAExpression/si_HuR_mRNASEQ/transcript_expression", 
#	skip=1);
#colnames(D.HuR)=c("id", "counts", "length", "expression");

#controlFile="~bilebi00/_KNOCKDOWNS/data/MRNAExpression/si_GFP_mRNASEQ/transcript_expression";
#D.GFP=read.delim(controlFile, skip=1);
#colnames(D.GFP)=c("id", "counts", "length", "expression");

#D=merge(D.HuR, D.GFP, by.x=1, by.y=1, all=T);
#D[is.na(D)]=0;

#D.qn=normalizeQuantiles(D[, c(4, 7)]); # columns 4 and 7 are the ones containing the expression levels

#file <- "293wt-RNAseq_HEK293_Nanofectin_RNAseq_HEK293_siAUF1_RNAseq_HEK293_siGFP_RNAseq_HEK293_siTIA1_RNAseq_HEK293_sihnRNPC_RNAseq_mRNASeq_No_4SU_No_XL_rep_A_Clip13_mRNASeq_No_4SU_No_XL_rep_B_Clip13_si_GFP_mRNASEQ_si_HuR_mRNASEQ";

	if (length(grep(".raw", file)) == 0 ) {
		file = paste(file,'.raw', sep='');
	}
	library(limma);
	D=read.table(file, header=T, stringsAsFactors=T, check.names=F);

	D.qn=normalizeQuantiles(D[,2:dim(D)[2]]);
	header <- F;
	if (wheader == 1) {
		header <- colnames(D);
	}

	write.table(cbind(D[,1],D.qn), sep="\t", row.names=F, col.names=header, quote=F, file=ofile);
	#write.table(cbind(D[,1],D.qn), sep="\t", row.names=F, col.names=colnames(D), quote=F, file=ofile);
}

cargs <- commandArgs()
file <- cargs[4];
ofile <- paste(file, '.qn', sep='');
wheader <- 1;
if (length(cargs) == 5) {
	wheader <- cargs[5];
}
qn(file,ofile,wheader)

