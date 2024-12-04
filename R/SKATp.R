SKATp <- function( data, pheno, covlist, binvector, pheno_type, geno_type ,is_dosage, SKAT_kernel, SKAT_weights, SKAT_weights_beta )
{
	if(geno_type=="D") { is_dosage=TRUE 
    } else { is_dosage=FALSE }
	genotype <- as.matrix(data[,names(binvector),drop=FALSE])
	if(!is.null(covlist)) {
		n.model <- as.formula(paste(pheno, "~" ,	paste(covlist,collapse="+")))
	} else {
		n.model <- as.formula(paste(pheno, "~ 1"))
	}
	
	SKAT.null <- SKAT::SKAT_Null_Model( n.model,   data=as.data.frame(data), out_type=pheno_type )
	SKAT.r <- SKAT::SKAT( genotype, SKAT.null, is_dosage=is_dosage,  kernel=SKAT_kernel, weights.beta=SKAT_weights_beta, weights=SKAT_weights )
	SKATw.r <- SKAT::SKAT( genotype, SKAT.null, is_dosage=is_dosage, kernel=SKAT_kernel,  weights.beta=SKAT_weights_beta, weights=SKAT_weights, method="optimal.adj")	
	return( list( SKAT.p=SKAT.r$p.value, SKATO.p=SKATw.r$p.value ) )
}