processing<-function(data, SNPinfo, pheno, covlist, mafcut, multiallelic) { 

	removed.multiallelic<-removed.mafcut<-bialSNP<-NULL
	multialSNP.setA<-multialSNP.setB<-multialSNP.left<-check<-multialSNP.kept<-NULL
	nSNPs<-nSNPs.mafcut<-NA
	SNPinfo$alfreq <- apply(data[,SNPinfo$variant], 2, function(x) { mean(x)/2 })
	SNPinfo$maf <- ifelse(SNPinfo$alfreq>0.5, 1-SNPinfo$alfreq, SNPinfo$alfreq)
	nSNPs <- length(SNPinfo$variant)
	removed.mafcut <- subset(SNPinfo, maf<mafcut)$variant
	SNPinfo <- subset(SNPinfo, maf>=mafcut)
	nSNPs.mafcut <- length(SNPinfo$variant)
	
	if(!is.null(SNPinfo)) { 
		bialSNP<-subset(SNPinfo, multiallelic==0)$variant
	
		if(isFALSE(multiallelic)) {
			removed.multiallelic<-subset(SNPinfo, multiallelic==1)$variant
      SNPinfo<-subset(SNPinfo, variant%in%bialSNP)
		}
	
		data <- data[,c(pheno, covlist, SNPinfo$variant), drop=FALSE]
		torec <- subset(SNPinfo, alfreq>0.5 & multiallelic==0)$variant
		if(length(torec)>0) { data[,torec] <- 2-data[,torec] }
    SNPinfo$maf <- apply(data[,SNPinfo$variant], 2, function(x) { mean(x)/2 })
	  		
  	if(isTRUE(multiallelic)) { 	#setA_multi, setB_multi=NULL
  		multialSNP.left <- subset(SNPinfo, multiallelic==1)
  		multialSNP.kept <- multialSNP.left$variant
  		if(!is.null(multialSNP.left)){
  			multialSNP.setB <- multialSNP.left[c(duplicated(multialSNP.left$bp), 
          duplicated(multialSNP.left$bp, fromLast=TRUE)), "variant"]
  			multialSNP.setA <- setdiff(multialSNP.left$variant, multialSNP.setB) #1 variable left
  		}
  	  } 
  }
	return(list(data=data, SNPinfo=SNPinfo, bialSNP=bialSNP, 
		nSNPs=nSNPs, nSNPs.mafcut=nSNPs.mafcut,
		removed.multiallelic=removed.multiallelic, removed.mafcut=removed.mafcut,
		multialSNP.kept=multialSNP.kept, multialSNP.setA=multialSNP.setA, multialSNP.setB=multialSNP.setB))
}