MLC<-function(beta=NULL, sigma=NULL, sigmainv, Z=NULL, invcor=NULL, 
	binvector, codechange, tol){
	V_delta_inverse<-NULL
	bin<-unique(binvector)
	nbin<-length(bin)
	J <- matrix( 0, ncol=nbin, nrow=length(binvector) )
	binsize <- unlist( lapply(unique(binvector), 
			function(x) { length(binvector[which(binvector==x)])}) )
	cumsize <- c( 0, cumsum(binsize) )
		
	for(b in 1:nbin) {
		sizec <- binsize[b]
		sizecum <- cumsize[b]
		J[ (sizecum+1) : (sizecum+sizec), b] <- rep(1, sizec) 
	}
	
	if(!is.null(Z) && !is.null(invcor)){
		if(any(codechange==1)) {
        Z[names(which(codechange==1))]<- -Z[names(which(codechange==1))]
		    invcor[names(which(codechange==1)),]<- -invcor[names(which(codechange==1)),]
		    invcor[,names(which(codechange==1))]<- -invcor[,names(which(codechange==1))]
		}		
		W <- (invcor%*%J) %*% solve(t(J)%*% invcor%*% J,tol=tol)
		deltabin <- t(W)%*%Z
		V_delta <- solve(t(J)%*%invcor%*%J)
		deltabin_names <- c("bin","NSNPs.kept","deltaZ","deltaZ.se","deltaZ.pvalue")
	
	} else {
    if(any(codechange==1)) {
    		beta[names(which(codechange==1))] <- -beta[names(which(codechange==1))]
		    sigmainv[names(which(codechange==1)),] <- -sigmainv[names(which(codechange==1)),]
		    sigmainv[,names(which(codechange==1))] <- -sigmainv[,names(which(codechange==1))]
	  }
     
		W <- (sigmainv%*%J)%*%solve( (t(J)%*%sigmainv%*%J) )
		deltabin <- t(W)%*%beta
		V_delta <- solve(t(J)%*%sigmainv%*%J)
		deltabin_names <- c("bin","NSNPs.kept","deltaB","deltaB.se","deltaB.pvalue")
	}
	
	tryCatch({V_delta_inverse=solve( V_delta, tol=tol )}, error=function(e){cat("ERROR :",	conditionMessage(e), "\n")})
	if(!is.null(V_delta_inverse) ){ 
			stat <- t(deltabin)%*%V_delta_inverse%*%deltabin
			pvalue <- pchisq(stat, df=nbin, lower.tail=F)
			
			deltabin.se <- sqrt(diag(V_delta))
			LCbin <- (deltabin/deltabin.se)^2
			deltabin.pvalue <- pchisq(LCbin,df=1,lower.tail=F)
				
			deltabin<-cbind(bin=bin,binsize=binsize, deltabin=deltabin,
			deltabin.se=deltabin.se,deltabin.pvalue=deltabin.pvalue)
			colnames(deltabin)<-deltabin_names
	} else { 
		 deltabin<-NULL
	}
		
   return(list(stat=stat, df=nbin, pvalue=pvalue, deltabin=deltabin))
}
