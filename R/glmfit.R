glmfit<-function(data, pheno, covlist=NULL, binvector, family, tol, firthreg) {
	invcov_g<-invcor_g<-NULL
	lmodel <- as.formula(paste(pheno, "~", 
		paste(c(names(binvector),covlist), collapse= "+")))
	lmodelvif <- as.formula(paste(pheno, "~", paste(c(names(binvector)), collapse= "+")))
			
	glmfit <- stats::glm(lmodel, data, family=family)
	
	if(firthreg==TRUE && family=="binomial") { glmfit <- glm(lmodel, data, family=family, na.action=na.omit, method = "brglmFit") }
	
	beta_g <- stats::coef(glmfit)[names(binvector)]
	vcov_g <- stats::vcov(glmfit)[names(binvector),names(binvector),drop=FALSE]
	cor_g<-  summary(glmfit,correlation=T)$correlation[names(binvector),names(binvector),drop=FALSE]
	beta_SE <- sqrt(diag(vcov_g))
	Z_g<-beta_g/beta_SE
	tryCatch({ invcov_g <- solve(vcov_g, tol=tol)},  error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
	tryCatch({ invcor_g <- solve(cor_g, tol=tol)},  error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
	p_g <- coef(summary(glmfit))[names(binvector),4]
	
	if(!is.null(covlist)) {
		beta_covlist <- stats::coef(glmfit)[covlist]
		vcov_covlist <- diag(stats::vcov(glmfit)[covlist,covlist])
		beta_covlist_SE<-sqrt(vcov_covlist)
		p_covlist <- coef(summary(glmfit))[covlist,4]
	} else {
		beta_covlist<- beta_covlist_SE<-p_covlist<- NA
	}
	
	vif_g <- rms::vif(glm(lmodelvif, data, family=family))
	
	return(list(beta_g=beta_g, beta_SE=beta_SE, Z_g=Z_g,
		vcov_g=vcov_g, invcov_g=invcov_g,invcor_g=invcor_g, p_g=p_g, vif_g=vif_g, beta_covlist=beta_covlist,
		beta_covlist_SE=beta_covlist_SE, p_covlist=p_covlist))
}