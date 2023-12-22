LCbin <- function( beta, cov, cor, binlist, code=rep(1,length(beta)), tol=tol)
{
	LCbinout <- matrix(0,length(binlist),5)
	LCBout <- LCZout <- list()
 	for( i in seq_along(binlist) ){
		betasub <- beta[binlist[[i]]]
		covsub <- cov[binlist[[i]],binlist[[i]]]
		corsub <- cor[binlist[[i]],binlist[[i]]]
		codesub <- code[binlist[[i]]]
		Bout <- MLCB( beta=betasub, cov=covsub, 
      code=codesub, tol=tol )
		Zout <- MLCZ(  beta=betasub, cov=covsub, 
      cor=corsub, code=codesub, tol=tol )
		LCbinout[i,] <- c( length(binlist[[i]]), Bout$stat, 
      Bout$pvalue, Zout$stat, Zout$pvalue ) 
		LCBout <- rbind(LCBout,Bout) 
	}
	colnames(LCbinout) <- c( "nvariants", 
    "LCBbin", "LCBbin.p", 
    "LCZbin", "LCZbin.p" )
	LCBout <- base::cbind(1:length(binlist),
    length(binlist), LCBout)
    return(list(LCbinout,LCBout))
}