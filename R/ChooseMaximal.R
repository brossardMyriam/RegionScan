ChooseMaximal <- function(data, vt, edgecut)
{
	codeW <- CodeChangeV(data,vt)$rout
	codeW[codeW<edgecut] <- 0
	subg <- igraph::graph.adjacency( codeW, mode="undirected", 
		weighted=TRUE, diag=FALSE, add.colnames=NULL )
	lgstcliq <- igraph::largest.cliques(subg)
	if(length(lgstcliq)==1){
		FC <- unlist(igraph::largest.cliques(subg))
	}else{
		sumvt <- sapply(lgstcliq,function(x) {sum(codeW[x,x])})
		cliqno <- which(sumvt==max(sumvt))
		FC <- lgstcliq[[cliqno]]
	}
  return(vt[FC])
}