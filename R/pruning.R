pruning <- function( binvector, data, rcut,	multialSNP.setA=NULL, multialSNP.setB=NULL )
{ 
	removed <- newbinvector <- NULL
   
	for( bin in unique(binvector) ) { 
  	bremoved <- NULL
  	bkept <- names( binvector[ binvector==as.character(bin) ] )
     
  	if( length( bkept ) > 1 ) {
  			rbin <- abs( cor(data[, bkept, drop=FALSE], use="pairwise.complete.obs") )
  			diag( rbin ) <- 0
  			prox_length <- lapply( bkept, function(snp) { length(rbin[snp,][which(rbin[snp,]>=rcut)]) } )
   			names( prox_length ) <- bkept
 
        bkept.setA <- bkept[which(bkept%in%multialSNP.setA)]
        bkept.setB <- bkept[which(bkept%in%multialSNP.setB)]
        
        if( length( bkept.setA ) > 0 ) {
          prox.setA <- prox_length[ bkept.setA ]
          names(prox.setA) <- bkept.setA
          maxcs <- unlist( prox.setA )[ which.max( unlist( prox.setA )) ]
  			  while(length(maxcs)>0 && maxcs>0){
  				  bremoved <- c( bremoved , names(maxcs) )
  				  bkept <- setdiff( bkept , bremoved )
            bkept.setA <- setdiff( bkept.setA , bremoved )        
  				  rbin <- rbin[ bkept, bkept, drop=FALSE ]
                   
  				  prox.setA <- lapply( bkept.setA, function(snp) { 
  					  length( rbin[snp,][ which(rbin[snp,]>=rcut) ]) })
  				  names( prox.setA ) <- bkept.setA
  				  maxcs <- unlist( prox.setA )[which.max( unlist( prox.setA )) ]
  			  }
        }
        
   			# then prun biallelic SNPs
        rbin <- abs( cor(data[,bkept, drop=FALSE], use="pairwise.complete.obs") )
  			diag( rbin ) <- 0
  			prox_length <- lapply( bkept, function(snp) { length(rbin[snp,][which(rbin[snp,]>=rcut)]) })
  			names( prox_length ) <- bkept
        maxcs <- unlist( prox_length )[ which.max( unlist( prox_length)) ]
         
       	  if( length( bkept.setB ) == 0 ) {
          while(maxcs>0){
    				bremoved <- c( bremoved , names(maxcs) )
    				bkept <- setdiff( bkept , bremoved )
    				rbin <- rbin[ bkept, bkept, drop=FALSE ]
    				prox_length <- lapply( bkept, function(snp) { 
    					length( rbin[snp,][ which(rbin[snp,]>=rcut) ]) })
    				names( prox_length ) <- bkept
    				maxcs <- unlist( prox_length )[which.max( unlist( prox_length )) ]
    			}
        } else if (length( bkept.setB ) > 0) {
             prox_length <- lapply( bkept.setB, function(snp) { length(rbin[snp,bkept.setB][which(rbin[snp,bkept.setB]>=rcut)]) })
		         names( prox_length ) <- bkept.setB
             maxcs <- unlist( prox_length )[ which.max( unlist( prox_length ) ) ]
             while(maxcs>0){
    			    	bremoved <- c( bremoved , names(maxcs) )
    				    bkept <- setdiff( bkept , bremoved )
    				    rbin <- rbin[ bkept, bkept, drop=FALSE ]
    				    prox_length <- lapply( bkept, function(snp) { 
    					      length( rbin[snp,][ which(rbin[snp,]>=rcut) ]) })
    				    names( prox_length ) <- bkept
    				    maxcs <- unlist( prox_length )[which.max( unlist( prox_length )) ]
 			        }
            rbin <- abs( cor(data[,bkept], use="pairwise.complete.obs") )
  			    diag( rbin ) <- 0
  			    prox_length <- lapply( bkept, function(snp) { length(rbin[snp,][which(rbin[snp,]>=rcut)]) })
  			    names( prox_length ) <- bkept
            maxcs <- unlist( prox_length )[ which.max( unlist( prox_length)) ]
            while( maxcs>0 ){
    				    bremoved <- c( bremoved , names(maxcs) )
    				    bkept <- setdiff( bkept , bremoved )
    				    rbin <- rbin[ bkept, bkept, drop=FALSE ]
    				    prox_length <- lapply( bkept, function(snp) { 
    					    length( rbin[snp,][ which(rbin[snp,]>=rcut) ]) })
    				    names( prox_length ) <- bkept
    				    maxcs <- unlist( prox_length )[which.max( unlist( prox_length )) ]
    			  }
  	}
	}
  bkeptvector <- rep(bin, length(bkept))
	names(bkeptvector) <- bkept
	newbinvector <- c(newbinvector, bkeptvector)
	bremovedvector <- rep(bin, length(bremoved))
	names(bremovedvector)<-bremoved
	removed<-c(removed,bremovedvector)
 }
	return( list(binvector=newbinvector , removed=removed ))
}
