CLQ <- function(data, edgecut=edgecut)
{
	binlist <- NULL
	snplist<-colnames(data)
	r<-cor(data)
	absr<- abs(r)
	absr[absr<=edgecut] <- 0
	
	g <- igraph::graph.adjacency(absr, mode="undirected", weighted=TRUE, diag=FALSE)
	clique <- igraph::maximal.cliques(g)
	wcl <- NULL
	
	# cliques for each snp
	for (snp in snplist){
		k <- lapply(clique, function(x) is.element(snp,names(x)))
		cl <- clique[k==TRUE]
		ll <- sapply(cl, length)
		mc <- cl[which.max(ll)]
		wcl <- append(wcl,mc)
	}
	wcl <- unique(wcl)
	wcl <-lapply(wcl, function(x) { unlist(names(x))})
	
	while(length(snplist)>0){
    Csize <- unlist(lapply(wcl, length))
	if(max(Csize)>1){
		maxC <- wcl[which(Csize==max(Csize))]
		coderst <- lapply(maxC,function(x){ChooseMaximal(data,x,edgecut)}) 
		clstsize <- lapply(coderst,function(x) length(unique(x)))
		mccs <- max(unlist(clstsize)) #maximal coded cluster size
	if (length(clstsize[which(clstsize==mccs)])==1){
        maxclq <- coderst[which(clstsize==mccs)]
      }else{
        clqsum <- lapply(coderst,function(x){sum(r[unlist(x),unlist(x)])})
        maxclq <- coderst[which.max(clqsum)]
      }
		binlist <- append(binlist,maxclq)
		snplist <- setdiff(snplist,unlist(maxclq))			
		wcl <- lapply(wcl, function(x) setdiff(unlist(x),unlist(maxclq)))
		wcl[which(sapply(wcl,length)==0)] <- NULL
		wcl <- unique(wcl)
    }else{
	  	binlist <- append(binlist,wcl)
	  	snplist <- setdiff(snplist,unlist(wcl))
    }
  }
	return(binlist)
}