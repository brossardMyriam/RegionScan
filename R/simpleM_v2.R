simpleM_v2 <-function (x, cor_G, pca_cut_perc = 0.995)
{
    if (corpcor::is.positive.definite(cor_G) == FALSE) { cor_G<-COMBAT::ld.Rsquare( cor_G) }
    min_p_obs <- min(x)
    num_of_snps <- length(x)
    cor_r <- cor_G
    eigen_values <- eigen(cor_r, only.values = TRUE)$values
    eigen_values_sorted <- sort(eigen_values, decreasing = TRUE)
    sum_eigen_values <- sum(eigen_values_sorted)
    M_eff_G <- 1
    for (k in 1:num_of_snps) {
        temp <- sum(eigen_values_sorted[1:k])/sum_eigen_values
        if (temp >= pca_cut_perc) {
            M_eff_G <- k
            break
        }
    }
    p_simpleM <- 1 - (1 - min_p_obs)^M_eff_G
	# if(p_simpleM==0) { p_simpleM<- min_p_obs*M_eff_G }
	return(list(simpleM.df=M_eff_G, simpleM.p=p_simpleM))
}