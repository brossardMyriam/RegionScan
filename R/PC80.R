PC80 <- function( data=data, pheno=pheno, covlist=covlist, binvector=binvector, family=family, tol=tol) 
{
	gen <- data[,names(binvector),drop=F]
	pc <- prcomp(gen, scale.=T) 
	corone <- cor(gen,use="pairwise.complete.obs")
	ei <- eigen(corone)
	cumei <- NULL
	for(m in 1:length(ei$values)){
		cumei[m]<- sum(ei$values[1:m])
	}
	cumei <- cumei/cumei[m]
	ind <- 1
	if(cumei[1]<0.8) { 
		ind <- c(1:length(names(gen)))[cumei<0.8]
		ind <- c(ind,length(ind)+1)   
	}
	pcsub <- pc$x[,ind,drop=F]
	pcload <- pc$rotation[,ind,drop=F]
	PCnames80 <- paste("Region_PC",ind,sep="")
	colnames(pcsub)<-PCnames80
	
	data80 <- data.frame(cbind(data,pcsub))
	lmodel <- as.formula(paste(pheno, "~", paste(c(PCnames80,covlist), collapse= "+")))
	glmfit <- glm(lmodel, data80, family=family)
	beta_g <- coef(glmfit)[PCnames80]
	vcov_g <- vcov(glmfit)[PCnames80,PCnames80]
	beta_SE <- sqrt(diag(vcov_g))
	invcov_g <- solve(vcov_g, tol=tol)
	p_g <- coef(summary(glmfit))[PCnames80,4]
	
	PC80one <- Wald(beta_g,invcov_g)
	return(list(stat=PC80one$stat,df=PC80one$df,pvalue=PC80one$pvalue, PC80load=pcload))
}