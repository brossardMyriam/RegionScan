# Modified version from package 'COMBAT'
## added on 2022-6-18 to deal with special cases of non positive def matrix & Meff>nsnp
ext_simes_v2<-function(x, rmat) {
    eff.fun <- function(ldmat) {
        ldmat <- as.matrix(ldmat)
        eff.local <- nrow(ldmat)
        if (eff.local <= 1) return(1)
        ev <- eigen(ldmat, only.values = TRUE)$values
        if (sum(ev < 0) != 0) {
            ev <- ev[ev > 0]
            ev <- ev/sum(ev) * eff.local
        }
        ev <- ev[ev > 1]
        if(sum(ev - 1)<eff.local)  { 
			return(eff.local - sum(ev - 1)) 
		} else { 
			return(eff.local)
		}
	}
    eff.global <- eff.fun(rmat)
    n_values <- length(x)
    candid <- sapply(1:n_values,
        function(i) {
        (eff.global * x[i])/eff.fun(rmat[1:i,1:i]) })
    p_ext_simes <- min(candid)
    return(list(p=p_ext_simes, Meff=eff.global))
}