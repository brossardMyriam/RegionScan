assigning <-  function( binvector, codechange, 	multialSNP.kept,  data, edgecut ) {
	rcutprunout	<- binmap <- NULL
  nbin <- max(binvector)
  recbial<-data[, c(names(binvector),multialSNP.kept) ]
  if(any(codechange==1)) recbial[,names(codechange==1)]<-2-recbial[,names(codechange==1)]
  rmat<-cor(recbial)
	for(snp in multialSNP.kept)  { 
			rbin<-lapply(unique(binvector), 
        function(bin) { 	binn<-names(binvector[which(binvector==bin)]) 
        return(rmat[snp,binn,drop=FALSE])})
			rbinmin<-unlist(lapply(rbin,function(x) x[which.min(x)]))
			names(rbinmin)<-unique(binvector)
      maxrbin<-max(rbinmin)>=edgecut
   if(!isTRUE(maxrbin)) {  
        binmap<-which.max(rbinmin)
        binvector<-c(binvector, binmap) 
        names(binvector)[length(binvector)]<-snp
			} else {
			  binvector<-c(binvector, nbin+1) 
        names(binvector)[length(binvector)]<-snp	
      }
   }
  return(binvector=sort(binvector))
} 