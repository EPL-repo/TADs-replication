# 
#  Example code to generate a tab-separated file with sigma values (with noise subtraction,
#  first variable/column for each chromosome) from an EdUseq-HU CSV output file
#  
#  The file NEcycEHU14_bin-size_10000_quality_37_chr1-X_sigma_all_0b.csv was generated   
#  after processing EdUseq-HU data from U2OS cells expressing normal levels of cyclinE
#  and can be found in https://github.com/EPL-repo/TADs-replication/CSV_files/EdUseq/U2OS
#
#  It contains 25,000 lines (10-kb bins) and 14 variables per chromosome; chromosome
#  data blocks are separated by a blank column
#

tab <- read.csv("NEcycEHU14_bin-size_10000_quality_37_chr1-X_sigma_all_0b.csv",h=F)

tbl <- data.frame(matrix(0,ncol=23,nrow=nrow(tab)))
for(i in 1:23){
	tbl[,i] <- tab[,(i-1)*15+1]
}

colnames(tbl) <- paste0("sigma",1:23)

tbl[,24] <- seq(from=0,to=nrow(tab)-1,by=1)
tbl[,25] <- seq(from=0,to=(nrow(tab)*10000)-10000,by=10000)
tbl[,26] <- seq(from=9999,to=(nrow(tab)*10000),by=10000)

colnames(tbl)[24:26] <- c("bin","start","end")

start <- rep(tbl$start,23)
end <- rep(tbl$end,23)
bin <- rep(tbl$bin,23)
chr <- NULL
sig <- NULL
for(i in 1:23){
	chr <- c(chr,rep(paste0("chr",i),nrow(tbl)))
	sig <- c(sig,tbl[,i])
} 

new.tab <- data.frame(chr = chr,
					  start = start,
					  end = end,
					  sigma = sig,
					  bin = bin)

options(scipen = 999)
 
write.table(new.tab,
			"NEcycEHU14_bin-size_10000_quality_37_chr1-X_sigma_all_0b_test.tsv",
			quote=F,
			sep="\t",
			row.names=FALSE)

