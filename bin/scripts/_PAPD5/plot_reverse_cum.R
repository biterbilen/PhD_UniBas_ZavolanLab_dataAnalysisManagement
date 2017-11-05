cargs <- commandArgs();
library(lattice)

cargs <- c(1,1,1,"raw_");

files <- dir(path=".", pattern=cargs[4], all.files=T);
files 

