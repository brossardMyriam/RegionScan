gates_v2<- function(x, cor_G)
{
   
	if (corpcor::is.positive.definite(cor_G) == FALSE) { cor_G<-COMBAT::ld.Rsquare(cor_G) }
   pval_sort <- sort(x)
    pval_order <- order(x)
    n_snps <- length(x)
    cor_P <- cor_G[pval_order, pval_order]
    cor_P <- 0.2982 * cor_P^6 - 0.0127 * cor_P^5 + 0.0588 * cor_P^4 +
        0.0099 * cor_P^3 + 0.6281 * cor_P^2 - 9e-04 * cor_P
    p_gates <- COMBAT::ext_simes(pval_sort, cor_P)
    return(p_gates)
}