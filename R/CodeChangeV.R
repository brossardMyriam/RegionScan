CodeChangeV <- function(data, binlist){ 
	mvar <- NULL
	rin <-cor(data[,binlist, drop=F])
	diag(rin) <-0
	nr <-nrow(rin)
	code <- rep(0,nr)
	names(code) <- binlist
	iter <-0
	change <- 1

	if(length(binlist)>1) {
	#	if(!is.null(bin)) { message(paste("Recoding of the SNPs within bin", bin,"for MLC tests")) }
		while(change==1){
			change <- 0
			cneg <- apply(rin,1,function(x){ sum(x<0)})
			maxneg <- which.max(cneg) 
			
			if(cneg[maxneg]>(nr-1)/2){ 
				rin[maxneg,] <- -1*rin[maxneg,] 
				rin[,maxneg] <- -1*rin[,maxneg]
				code[maxneg] <- 1-code[maxneg] 
				change <- 1
			}
			iter <- iter+1
		}
		if(iter> 2*nr ) {
		  message(paste("WARNING::: Please check variants in bin ", bin, ". Some SNPs appear negatively correlated with more than 50% of the SNPs",sep="")) 
		}
	}
	diag(rin) <-1
	return(list(rout=rin,codechange=code))
}