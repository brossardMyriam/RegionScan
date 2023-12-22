aliasing <- function(pheno, data, binvector){
	calias <- NULL
  lmodel <- as.formula(paste(pheno, "~", paste(names(binvector), collapse= "+")))
	calias <- stats::alias(lmodel, data)$Complete
	if( !is.null(calias) ) {
		valiases <- rownames( calias )
		final <- setdiff( names(binvector),valiases )
		aliasout <- binvector[ which(names(binvector)%in%valiases) ]
		binvector <- binvector[ final ]
		final<-binvector
		#bin <- unique(binvector)
		#size <- unlist(lapply(bin, function(x) { length(binvector[which(binvector==x)])}))
		#names(size) <- bin
		#binsize <- sort(size, decreasing=T)
		#binsize2 <- 1:length(binsize)
		#names(binsize2) <- names(binsize)
		#binvectorf <- lapply(as.numeric(names(binsize2)), 
		#	function(bi) { names(which(binvector==bi)) })
		#binvectorf2 <- unlist(lapply(1:length(binvectorf), 
		#	function(c) { rep(c,length(binvectorf[[c]])) } ))
		#final <- binvectorf2
		#names(final) <- unlist(binvectorf)
	} else { final <- binvector ; aliasout <-NULL } 
	return(list( final=final, aliasout=aliasout ))
}