rename <- function(nms) {
	nms.new <- nms[1]
	for (i in 2:length(nms)) {
		count <- 0;
		for (j in 1:(i-1)) {
			if (nms[i] == nms[j]) {
				count <- count + 1;
			}
		}
		nms.new <- c(nms.new, ifelse(count>0, paste(nms[i], as.character(count),sep="."), nms[i]))
	}
	return(nms.new)
}

getDataFromGtfPlus <- function (filename, onames=NULL) {
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
	# names(a) <- c("seqname","source","feature","start","end","score","strand","frame", levels(unlist(ns)));
	names(a) <- rename(c("seqname","source","feature","start","end","score","strand","frame", as.vector(unlist(ns)), onames))
	return(a)
																	
}

#cargs <- commandArgs(trailingOnly = T)

#a <- getDataFromGtfPlus(cargs[1], plusnames);
#head(a)

