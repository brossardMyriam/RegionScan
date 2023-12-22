MinPJ <- function(tvalue, pvalue, corr)
{ #MinP from James
	n <- length(pvalue)
	if(n==1)return(list(statp=c(tvalue,pvalue),adjp=pvalue))
	gvector <- c(0.28209479, 0.14104740, 0.08578128, 0.05814822, 0.04224021,
			0.03219472, 0.02542143, 0.02062518, 0.01709725, 0.01442215,
            0.01234263, 0.01069224, 0.00935924, 0.00826625, 0.00735830,
            0.00659537, 0.00594782, 0.00539322, 0.00491441)
	gvalue <- gvector[length(gvector)]
	if(n<19)gvalue <- gvector[length(pvalue)]
	m <- n*(n-1)/2 
	rij <- abs(corr)
	rm <- ((sum(rij)-sum(diag(rij)))/2)/m	
	tempsum <- 0
	
	for(s in 2:n){ # general correlation structure ( Armitage and Parmap approx)
		for(t in 1:(s-1)){ tempsum <- tempsum+(rij[s,t]-rm)^2 }
	}		
	rmv <- rm+2*tempsum/m # rho 
 	nval <- qnorm(pvalue/2,0,1,lower.tail=F)
	pdf <- (1/exp(2*3.141592654))*exp(-nval*nval/2)
	d2 <- 1-pvalue
	d1 <- d2^n
	d3 <- 0
	d4 <- n*(n-1)*pdf*gvalue
	if( rmv>=1 ) { jp <- pvalue*n }
	if(rmv<1){
		jp <- 1-d1*(1-rmv**2)-d2*(rmv**2)-d3*rmv*(1-rmv)-d4*(2-2*sqrt(1-rmv)-rmv-rmv**2)
	} # approximation pf multi normal probabilities not guaranteed after more than 5 tests
	Pmin <- min(pvalue)
	jpmin <- min(jp)
	Tmax <- max(abs(tvalue))
	return(list(statp=c(Tmax,jpmin),adjp=jp))
}