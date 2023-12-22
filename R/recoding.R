recoding <- function(binlist, data) {
	cout <- lapply(seq_along(binlist), function(x) { CodeChangeV(data, binlist[[x]])}) 
	newcode <- unlist(lapply(cout,function(l) {  l$codechange }))
	rout<-lapply(cout,function(l) {  l$rout})
	return(list(newcode=newcode, rout=rout))
}