# 
#  Example R code to find local extrema in a vector
#  
#  Can be used to process EdUseq-HU outputs, e.g. by creating a vector with previously
#  calculated sigma values and then looking for local maximal sigma values; the obtained
#  peaks can be further filtered using a predefined sigma threshold value (such as 
#  positive values only)
#
#
#  Default step = 2 (bw=2) 
#

find.peaks <- function(x,bw=2,x.coo=c(1:length(x))){

	pos.x.max <- NULL
	pos.y.max <- NULL
	pos.x.min <- NULL
	pos.y.min <- NULL 	

	for(i in 1:(length(x)-1)){ 		
		if((i+1+bw)>length(x)){
			sup.stop <- length(x)}
		else{sup.stop <- i+1+bw}
		
		if((i-bw)<1){inf.stop <- 1}
		else{inf.stop <- i-bw}

		subset.sup <- x[(i+1):sup.stop]
		subset.inf <- x[inf.stop:(i-1)]

		is.max <- sum(subset.inf > x[i]) == 0
		is.nomin <- sum(subset.sup > x[i]) == 0

		no.max <- sum(subset.inf > x[i]) == length(subset.inf)
		no.nomin <- sum(subset.sup > x[i]) == length(subset.sup)

		if(is.max & is.nomin){
			pos.x.max <- c(pos.x.max,x.coo[i])
			pos.y.max <- c(pos.y.max,x[i])
		}

		if(no.max & no.nomin){
			pos.x.min <- c(pos.x.min,x.coo[i])
			pos.y.min <- c(pos.y.min,x[i])
		}
	}
	return(list(pos.x.max,pos.y.max,pos.x.min,pos.y.min))
}
