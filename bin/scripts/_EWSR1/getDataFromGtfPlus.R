getDataFromGtfPlus <- function (filename, onames) {
	aa <- read.table(gzfile(filename), header=F,sep="\t");

	groups <- matrix(aa$V9)
	aa$V9 <- NULL

	names(aa) <- c("seqname","source","feature","start","end","score","strand","frame", onames)

	#separate the groups
	bb <- matrix(unlist(lapply(groups, function(x) { strsplit(x,"; *")[[1]] } )), nrow=length(groups), byrow=T)

	#group names
	gnames <- unlist(lapply(t(bb[1,]), function(x) { sub(" +(.+)","",x,perl=T) }))

	#add groups to aa
	count <- 1
	for (g in gnames) {
		pat <- paste(g, " +", sep="");
		aa[,g] <- unlist(lapply(matrix(bb[, count]), function(x) { sub(pat,"",x,perl=T) }))
		count <- count + 1
	}
	return (aa)
}

getDataFromGtfPlus.back <- function (filename, onames) {
	aa <- read.table(gzfile(filename), header=F);

	head(aa)
	im <- 1:8
	io <- 9:dim(aa)[2]

	itp <- io[io %% 3==2]    # temporary array for the indices of punctuation marks
	la <- aa[1,itp] == ";"   # logical array for the attribute

	ns <- aa[1,itp[la]-2]    # names of attribute
	iv <- io[io %% 3 == 1][1:length(la[la])] # values of attributes

	ip <- NULL #in case the file is not gtf plus but gtf
	if ( dim(aa)[2] > itp[la][length(itp[la])] ) {
		ip <- (itp[la][length(itp[la])]+1):(dim(aa)[2])
	}
	a <- aa[,c(im,iv,ip)];
	head(a)
	# names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)));
	names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)), onames)
	return(a)

}

#cargs <- commandArgs(trailingOnly = T)

#mutrate <- -1
#plusnames <- NULL
#if (length(cargs) > 1) {
#	mutrate <- as.numeric(cargs[2]);
#}

#a <- getDataFromGtfPlus(cargs[1], plusnames);
