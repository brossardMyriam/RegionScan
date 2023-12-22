Wald <- function( beta, invcov, df=NULL )
{
	if( is.null( df ) ) { df <- length( beta ) } 
	if( !is.null( invcov ) ){
		ans <- t(beta) %*% invcov %*% beta
		pvalue <- pchisq( ans, df=df, lower.tail=F )	
	}
	else {
		ans <- NA
		pvalue <- NA
		}
	return( list( stat=ans, df=df, pvalue=pvalue ) )
}