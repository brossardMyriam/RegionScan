clustering<-function(data, edgecut) {
		binvector <- NULL
		snplist <- as.character(colnames(data))
		binlist <- CLQ(data, edgecut=edgecut) 
		binvector <- unlist(lapply(1:length(binlist), function(c) { rep(c,length(binlist[[c]])) }))
		names(binvector) <- unlist(binlist)
		return(list(binvector=binvector, binlist=binlist))
} 
		
