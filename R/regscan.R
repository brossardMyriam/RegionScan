#' regscan : main function for genome-wide region-level analysis
#' 
#' @param phenocov a dataframe containing the phenotype(s) and  covariates if applicable
#' @param pheno the name of the phenotype column in phenocov to be analyzed
#' @param regionlist list of characters with region names
#' @param REGIONinfo dataframe including region definition, with 4 columns: chr”, “region”, “start.bp”, “end.bp”, as provided for example by the output of “BigLD” function from the gpart R package
#' @param geno_type (required): geno_type=“D” if genotypes are in allele dosage format , or genotype format (geno_type=“G”) ; this argument required for SKAT/SKATO tests
#' @param pheno (required): name of the phenotype column in phenocov to be analyzed; must be included in phenocov input
#' @param pheno_type (required): pheno_type=“C” if phenotype is continuous, or pheno_type=“D” if dichotomous
#' @param regionlist (optional): subset or list of the region names to analyze
#' @param covlist (optional): list of the covariate names to include in the analyses (must be part of phenocov input)
#' @param covout (default=FALSE): specify if covariates estimates should be reported ; will produce an additional output
#' @param SNPinfo: dataframe including the SNP information,  must include the following columns: "chr", "pos", "variant", "ref", "alt", “multialSNP”, "multialSNP.rec","freq.alt","multialSNP": binary variable=0, if the SNP is bi-allelic (two alleles), or 1 if more than two alleles"multialSNP.rec": indicates if the multi-allelic SNP has been recoded b recodeVCF (such that baseline allele = major allele), “freq.alt": frequency of the alternate allele (minor allele). This input can be generated using the auxiliary function recodeVCF.
#' @param vcfname:NULL 
#' @param machr2:NULL
#' @param qcinput:NULL 
#' @param info_score:NULL
#' @param multiallelic (default=FALSE): If FALSE, extract & process only the biallelic SNPs; if TRUE: include multiallelic SNPs in addition of the biallelic SNPs
#' @param multial_nmaxalleles:For extraction of SNPs with less than multial_nmaxalleles alleles (e.g. 2 to extract biallelic SNPs (default), 3 to extract bi-allelic and triallelic SNPs, etc)
#' @param mafcut:0.05
#' @param rcut:parameter to prun out the variants within a region (default set to 0.99)
#' @param firthreg: (default is FALSE)
#' @param MLCheatmap:TRUE/FALSE to produce heatmap of the correlation within region (default is FALSE)
#' @param regionlist:specify a list of regions to analyze (optional)
#' @param alltests:by default include all the default region-level tests ; to include other region-level tests (set to TRUE)
#' @param edgecut:clustering CLQ cut parameter (default is 0.5)
#' @param tol:1e-16
#' @param qcmachr2: threshold used to filter out SNPs with low Mach R2 imputation quality score
#' @param covout:FALSE 
#' @param LDpruning:TRUE
#' @param singleSNPall:FALSE
#' @param SKAT_kernel: type of kernel used for SKAT and SKATO tests. 
#' @param SKAT_weights_beta: a numeric vector of parameters for the beta weights for the weighted kernels. If you want to use your own weights, please use the “weights” parameter. It will be ignored if “weights” parameter is not null. 
#' @param verbose:FALSE
#' @param debug:FALSE
#' @param parallel:FALSE
#' @param region_tests: specify the list of region-level tests to plot (default is NULL, all tests will be plotted)
#' @param outname: name of the output file produced 
#' @export
#' @return
#' @export
#'
#' @examples
# analysis for a quantitative trait with default parameter values

regscan <- function(phenocov=NULL, pheno, REGIONinfo, geno_type, pheno_type,  
	covlist=NULL, data=NULL, SNPinfo=NULL, vcfname=NULL, machr2=NULL, qcinput=NULL, info_score=NULL,
	multiallelic=FALSE, multial_nmaxalleles=2, mafcut=0.05, rcut=0.99, firthreg=FALSE, MLCheatmap=FALSE, 
	regionlist=NULL, alltests=FALSE, edgecut=0.5, tol=1e-16, qcmachr2=NULL, 
	covout=FALSE, LDpruning=TRUE, singleSNPall=FALSE, 
	SKAT_kernel="linear.weighted", SKAT_weights=NULL, SKAT_weights_beta=c(1,25),
	verbose=FALSE, debug=FALSE, parallel=FALSE)
{ 
	if(is.null(pheno_type)) { return("Please specify if the phenotype is continuous (pheno_type='C') or dichotomous (pheno_type='D')") }
	if(is.null(geno_type)) { return("Please specify if the genotypes are in allele dosage format (geno_type='D') or genotype format (geno_type='G')") }
	if( pheno_type=="C" ) family="gaussian"
	if( pheno_type=="D" ) family="binomial"
	options(warn=-1)
	# regionout <- binout <- variantout <- rmvariantout <- res <- NULL
	message("*********************************") 
	message("Input parameters") 
	message(paste("pheno:",pheno,sep=""))
	message(toString(paste("covlist:",toString(covlist))))
	message(paste("mafcut:",mafcut,sep=""))
	message(paste("edgecut:",edgecut,sep=""))
	message(paste("rcut:",rcut,sep=""))
	message(paste("geno_type:",geno_type,sep=""))
	message(paste("pheno_type:",pheno_type,sep=""))
	message(paste("multiallelic:",multiallelic,sep=""))
	message(paste("firthreg:",firthreg,sep=""))
    message(paste("covout:",covout,sep=""))
	message(paste("LDpruning:",LDpruning,sep=""))
	message(paste("alltests:",alltests,sep=""))
	message(paste("tol:",tol,sep=""))
	message(paste("singleSNPall:",singleSNPall,sep=""))
	message(paste("MLCheatmap:",MLCheatmap,sep=""))
	message(paste("debug:",debug,sep=""))
	message("*********************************") 
	#checking inputs
	if(is.null(vcfname) && is.null(data)) {	message(cat("ERROR: vcfname or data is required")) ; break  } 
	if(is.null(SNPinfo)) { message(cat("ERROR: SNPinfo is required")) ; break  } 
	if(is.null(REGIONinfo)) { message(cat("ERROR: REGIONinfo is required")) ; break  } 
	if(is.null(phenocov)) {	message(cat("ERROR: phenocov is required")) ; break  } 
	if(is.null(pheno)) { message(cat("ERROR: pheno is required")) ; break  } 
	if(is.null(pheno_type)) { message(cat("ERROR: pheno_type is required")) ; break  } 
	if(is.null(geno_type)) { message(cat("ERROR: geno_type is required")) ; break  } 
	if(!is.null(regionlist)) { REGIONinfo <-subset(REGIONinfo, region%in%regionlist) } 
   
  inputs<-lapply(REGIONinfo$region, function(regionc) { 
      tryCatch({      
			regioncur <-  REGIONinfo[which(REGIONinfo$region==regionc),]
		    region <- as.character(regioncur["region"])
		    chr <- as.numeric(regioncur["chr"])
		    start <- as.numeric(regioncur["start.bp"])
		    end <- as.numeric(regioncur["end.bp"])
		    SNPinfo_region <- subset(SNPinfo, SNPinfo$bp>=start & SNPinfo$bp<=end)
			datapcov <- phenocov[,c(pheno,covlist),drop=FALSE]
			data_region <- na.omit(cbind(datapcov,data[,SNPinfo_region$variant,drop=F]))
			return(regionc=list(SNPinfo_region=SNPinfo_region, data_region=data_region, regioninfo=regioncur))
     }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  })
  rm(data); rm(REGIONinfo) 
  if (isTRUE(parallel)) {
	n.cores <- parallel::detectCores() - 1
	my.cluster <- parallel::makeCluster( n.cores,  type = "PSOCK")
	#print(my.cluster) #check cluster definition 
	doParallel::registerDoParallel(cl = my.cluster) #register it to be used by %dopar%
	foreach::getDoParRegistered() #check if it is registered (optional)
	foreach::getDoParWorkers() 
  } 

	rscan <- foreach::foreach( regionc=inputs, .combine="rbind", .multicombine=T, .verbose = verbose,
    .packages=c("igraph", "rms", "SKAT", "foreach", "COMBAT", "brglm2", "ggplot2", 
      "reshape2", "doSNOW", "karyoploteR" )) %dopar% {
    regionout <- binout <- snpout <- outremoved <- filterout <- outsingleSNPall <- prunout <- assignout <- covoutput<-covarout<-NULL
    regioninfo <-  regionc$regioninfo
	message(paste("Region", regioninfo$region, sep=":"))
    region <- as.character(regioninfo["region"])
	chr <- as.numeric(regioninfo["chr"])
	start <- as.numeric(regioninfo["start.bp"])
	end <- as.numeric(regioninfo["end.bp"])
	SNPinfo_region <- regionc[["SNPinfo_region"]]
	data_region <- regionc[["data_region"]]
    tryCatch({
   		if(!is.null(SNPinfo_region)) {  
			processout<- RegionScan:::processing(data=data_region, pheno=pheno, covlist=covlist, 
				SNPinfo=SNPinfo_region, mafcut=mafcut, multiallelic=multiallelic)
			nSNPs <-nSNPs.kept <- length(processout$SNPinfo$variant)
			if(isTRUE(debug)) { save(processout, pheno, covlist, 
			file=paste("chr", chr,"_", region, "_checking.Rdata", sep="")) } 
			if(isTRUE(nSNPs>1)) {
				clustout <- RegionScan:::clustering(data=processout$data[,processout$bialSNP], edgecut=edgecut) #clustering (biallelic SNPs only) 
			if(length(processout$multialSNP.kept)>0) { #multi-allelic
            codechange.bial <- RegionScan:::recoding(binlist=clustout$binlist, data=processout$data) # recoding (for assignation multial SNP)
		        assignout <- RegionScan:::assigning(binvector=clustout$binvector, codechange=codechange.bial$newcode, multialSNP.kept=processout$multialSNP.kept, data=processout$data, edgecut=edgecut)
            
            if(isTRUE(LDpruning)) {  # optional
				        prunout <- RegionScan:::pruning(binvector=assignout, data=processout$data, rcut=rcut, 
                    multialSNP.setA=processout$multialSNP.setA, 
                    multialSNP.setB=processout$multialSNP.setB ) 
					      aliasout <- RegionScan:::aliasing( pheno=pheno, data=processout$data, binvector=prunout$binvector )
				    } else {	# no pruning on LD
					    aliasout <- RegionScan:::aliasing( pheno=pheno, data=processout$data, binvector=assignout)
				    } 
            
		    codechange.MLC <- RegionScan:::codechange.bial$newcode[names(codechange.bial$newcode)%in%names(aliasout$final)]
            mual <- setdiff(names(aliasout$final),names(codechange.bial$newcode))
           	binvector.LC <- rep(1,length(codechange.MLC))
        	names(binvector.LC) <- names(codechange.MLC)
        	codechange.LC <- RegionScan:::recoding(binlist=list(names(binvector.LC)), data=processout$data)$newcode # recoding (for LC test)
            
            if(length(mual)>0) {
                codechange.mual <- rep(0, length(mual))
                names(codechange.mual) <- mual
                codechange.MLC <- c(codechange.MLC, codechange.mual)
                codechange.MLC <- codechange.MLC[ names(aliasout$final) ]
                codechange.LC <- c(codechange.LC, codechange.mual)
                codechange.LC <- codechange.LC[ names(aliasout$final) ]
                binvector.LC <- rep(1,length(codechange.MLC))
           }  
        } else { # no multi-allelic SNPs 
            if(isTRUE(LDpruning)) {  # optional
						prunout <- RegionScan:::pruning(binvector=clustout$binvector, data=processout$data, rcut=rcut) 
					    aliasout <- RegionScan:::aliasing( pheno=pheno, data=processout$data, binvector=prunout$binvector )
				    } else {	# no pruning on LD
					    aliasout <- RegionScan:::aliasing( pheno=pheno, data=processout$data, binvector=clustout$binvector )
            }
				binlist.MLC <- lapply(unique(aliasout$final), function(bin) { names(which(aliasout$final==bin)) })
		        codechange.MLC <- RegionScan:::recoding(binlist.MLC, processout$data)$newcode # recoding (for MLC test)
	       		binlist.LC <- list(unlist(names(aliasout$final)))
        		binvector.LC <- rep(1,length(unlist(aliasout$final)))
        		names(binvector.LC) <- unlist(names(aliasout$final))
        		codechange.LC <- RegionScan:::recoding(binlist.LC, processout$data)$newcode # recoding (for LC test)
       }
        
			  nSNPs.kept <- length(aliasout$final)
			  if(isFALSE(singleSNPall)) { 
					sgout <- data.frame(do.call("rbind", lapply(names(aliasout$final), function(x) {  
						binvector<- 1 ; names(binvector) <- x
						sglmout<-RegionScan:::glmfit(data=processout$data, pheno=pheno, binvector=binvector, covlist=covlist, family=family, tol=tol,firthreg=firthreg)
						c(variant=x,sglm.beta=unname(sglmout$beta_g), sglm.se=unname(sglmout$beta_SE), sglm.pvalue=sglmout$p_g) })))
                                                                                            
				} else {
			if(length(assignout)>0) { binvectorsg <- assignout 
			} else { 
				binvectorsg <- clustout$binvector }
				sgout <- data.frame(do.call("rbind", lapply(names(binvectorsg), function(x) {  
						  binvector <- 1 ; names(binvector) <- x
						  sglmout <- RegionScan:::glmfit(data=processout$data, pheno=pheno, binvector=binvector, covlist=covlist, family=family,tol=tol,firthreg=firthreg)
						  c(variant=x,sglm.beta=unname(sglmout$beta_g), sglm.se=unname(sglmout$beta_SE), sglm.pvalue=sglmout$p_g)})))
						
					outsingleSNPall <- data.frame(cbind(chr=chr, region=region, start.bp=start, end.bp=end, bin=binvectorsg, 
							processout$SNPinfo[match(names(binvectorsg),processout$SNPinfo$variant),c("bp","multiallelic","ref","alt","maf")], sgout))
				}
		
		#region-level tests
		glmout <- RegionScan:::glmfit(data=processout$data, pheno=pheno, covlist=covlist, binvector=aliasout$final, family=family, tol=tol, firthreg=firthreg)
   		MLCBout <- RegionScan:::MLC(beta=glmout$beta_g, sigmainv=glmout$invcov_g, binvector=aliasout$final, codechange=codechange.MLC, tol=tol)
		LCBout <- RegionScan:::MLC(beta=glmout$beta_g, sigmainv=glmout$invcov_g, binvector=binvector.LC,	codechange=codechange.LC, tol=tol)
		Waldout <- RegionScan:::Wald(glmout$beta_g, glmout$invcov_g)
		PC80out <-RegionScan:::PC80( data=processout$data, pheno=pheno, covlist=covlist, binvector=aliasout$final, family=family, tol=tol) 
		SKATout <-RegionScan:::SKATp( data=processout$data, pheno=pheno, covlist=covlist,pheno_type=pheno_type, geno_type=geno_type, binvector=aliasout$final,
		  SKAT_kernel=SKAT_kernel, SKAT_weights=SKAT_weights, SKAT_weights_beta=SKAT_weights_beta)
		
		regionout <-c(chr=chr, region=region, start.bp=start, end.bp=end, nSNPs=nSNPs, nSNPs.kept=nSNPs.kept, maxVIF=max(glmout$vif),
				Wald=Waldout$stat, Wald.df=Waldout$df, Wald.p=Waldout$pvalue, 
				PC80=PC80out$stat, PC80.df=PC80out$df, PC80.p=PC80out$pvalue,
				MLCB=MLCBout$stat, MLCB.df=MLCBout$df, MLCB.p=MLCBout$pvalue,
				LCB=LCBout$stat, LCB.df=LCBout$df, LCB.p=LCBout$pvalue, 
				SKAT.p=SKATout$SKAT.p, SKATO.p=SKATout$SKATO.p)
	    
    if( isTRUE(nSNPs.kept>1)) {	
  		sgsel<-subset(sgout,variant%in%names(aliasout$final))
  		corsel<-cor(processout$data[,sgsel$variant,drop=F])
  		pval<-as.numeric(as.character(sgsel$sglm.pvalue))
  		names(pval)<-sgsel$variant
  		out_simes <-  RegionScan:::ext_simes_v2(pval, corsel)
  		out_simpleM <- RegionScan:::simpleM_v2(pval, corsel) 
  		out_gates <- RegionScan:::gates_v2(pval, corsel)
  		minp_uncorected<-min(pval)
			regionout <- c(regionout,	simes.p=out_simes$p, simpleM.df=out_simpleM$simpleM.df, 
  			simpleM.p=out_simpleM$simpleM.p, GATES.p=out_gates, single_Wald.p=minp_uncorected)
     } else {
		regionout$maxVIF<-NA
     	regionout <- c(regionout,	simes.p=glmout$p_g , simpleM.df=1 , 
				simpleM.p=glmout$p_g ,GATES.p=glmout$p_g, single_Wald.p=glmout$p_g)
     }   

   	 if ( length(assignout)>0 ) { 
        bin.bfP<-lapply(unique(assignout), function(x) { names(which(assignout==x)) })
     } else {
        bin.bfP<-lapply(unique(clustout$binvector), function(x) { names(which(clustout$binvector==x)) })
     }
	 binstart.bfP.bp <- unlist( lapply(bin.bfP, function(x) { min(as.numeric(as.character(subset(processout$SNPinfo,variant%in%x)$bp)))}))
     binend.bfP.bp <- unlist( lapply(bin.bfP, function(x) {   max(as.numeric(as.character(subset(processout$SNPinfo, variant%in%x)$bp)))}))
     bin.afP<-lapply(unique(aliasout$final), function(x) { names(which(aliasout$final==x)) })
     binstart.afP.bp <- unlist( lapply(bin.afP, function(x) { min(as.numeric(as.character(subset(processout$SNPinfo, variant%in%x)$bp)))}))
     binend.afP.bp <- unlist( lapply(bin.afP, function(x) { max(as.numeric(as.character(subset(processout$SNPinfo, variant%in%x)$bp)))}))
     binsize.bfP<-unlist(lapply(bin.bfP, length))
     binsize.afP<-unlist(lapply(bin.afP, length))
     binout <- data.frame(cbind(chr=chr, region=region, start.bp=start, end.bp=end, 
     binstart.bfP.bp=binstart.bfP.bp, binend.bfP.bp=binend.bfP.bp,
     binstart.afP.bp=binstart.afP.bp, binend.afP.bp=binend.afP.bp,
       binsize.bfP=binsize.bfP, binsize.afP=binsize.afP, MLCBout$deltabin))
	 if(isTRUE(alltests)) {
			MLCZout <- RegionScan:::MLC(Z=glmout$Z_g, invcor=glmout$invcor_g,  
					binvector=aliasout$final, codechange=codechange.MLC, tol=tol)
			LCZout <- RegionScan:::MLC(Z=glmout$Z_g, invcor=glmout$invcor_g, binvector=binvector.LC, 
				codechange=codechange.LC, tol=tol)
			regionout <- c(regionout,	MLCZ=MLCZout$stat,MLCZ.df=MLCZout$df,MLCZ.p=MLCZout$pvalue,
				LCZ=LCZout$stat,LCZ.df=LCZout$df,LCZ.p=LCZout$pvalue)
			binout <- c(binout,MLCZout$deltabin)
		}
		snpout<-data.frame(cbind( chr=chr, region=region, start.bp=start, end.bp=end, 
				bin=aliasout$final, processout$SNPinfo[match(names(aliasout$final),
				processout$SNPinfo$variant),c("bp","multiallelic","ref","alt","maf")],
				MLC.codechange=codechange.MLC[names(glmout$beta_g)], 
				LC.codechange=codechange.LC[names(glmout$beta_g)],
        		subset(sgout,variant%in%names(aliasout$final)), 
				mglm.vif=glmout$vif, mglm.beta=glmout$beta_g, 
				mglm.se=glmout$beta_SE, mglm.pvalue=glmout$p_g ))
		
    	# removed list
			allremoved<-filterout<-NULL
			if(length(processout$removed.mafcut)>0) { allremoved<-data.frame(cbind(variant=processout$removed.mafcut, bin="NA", reason="mafcut")) }
			if(length(processout$removed.multiallelic)>0) { allremoved<-data.frame(rbind(allremoved, 
				cbind(variant=processout$removed.multiallelic, bin="NA", reason="multial"))) }
			if(length(names(prunout$removed))>0) { allremoved<-data.frame(rbind(allremoved, 
				cbind(variant=names(prunout$removed),bin=prunout$removed, reason="rcut"))) }
			if(length(names(aliasout$aliasout))>0) { allremoved<-data.frame(rbind(allremoved,	
				cbind( variant=names(aliasout$aliasout),bin=aliasout$aliasout, reason="alias"))) }
			if(!is.null(allremoved)) { 
				  bp <- SNPinfo_region[match(allremoved$variant,SNPinfo_region$variant),]$bp
				  names(bp) <- SNPinfo_region[match(allremoved$variant,SNPinfo_region$variant),]$variant
				  filterout <- cbind(chr=chr, region=region, start=start, end=end, bp=bp, allremoved)
			}
		
			if(!is.null(filterout) && nrow(filterout)>1) { filterout<-filterout[order(filterout$bin,filterout$bp),] }
			if(nrow(binout)>1) { binout<-binout[order(binout$bin),] }
			if(nrow(snpout)>1) { snpout<-snpout[order(snpout$bin,snpout$bp),] } 
			if(!is.null(outsingleSNPall)) { outsingleSNPall<-outsingleSNPall[order(outsingleSNPall$bin,outsingleSNPall$bp),] } 
    
      
     if (isTRUE(MLCheatmap)) {
        if(length(assignout)>0) { 
            binvector<-assignout
        } else {
            binvector<-clustout$binvector
        }
        
        nbins <- max(binvector)    
        colvect <- scales::hue_pal()(nbins)  
        outfile<-paste("MLC_heatmap_chr",chr,".region",region, sep="")
        
        tmp0<- subset(processout$SNPinfo, variant%in%names(binvector)) # before pruning by bp
        tmp0<-tmp0[order(tmp0$bp),]
        cormat1<-apply(cor(processout$data[,tmp0$variant,drop=F]),2, function(x) { round(x,2) } ) # all SNPs by bp
        tdata1<-processout$data[,tmp0$variant,drop=F]
        tdata1[,names(which(codechange.MLC==1))]<-2-tdata1[, names(which(codechange.MLC==1))] # recoded
        cormat<- apply(cor(tdata1),2, function(x) { round(x,2) } )
        tmp<-subset(processout$SNPinfo, variant%in%names(aliasout$final)) # after pruning & but not rec ordered by bp
        tmp<-tmp[order(tmp$bp),]
        tdata2<-processout$data[,names(aliasout$final)] # after pruning, ordered by bin
        tdata2[,names(which(codechange.MLC==1))]<-2-tdata2[, names(which(codechange.MLC==1))]
        cormat3<-apply(cor(tdata2),2, function(x) { round(x,2) } ) #ordered by bin & recoded after prunning
        RegionScan:::RegionHeatmap(cormat1, type_mat="Correlation,\n" , 
          ptitle="Within region correlation (before pruning & recoding),\n SNPs ordered by pos", 
            colvect=colvect, clustorder=F)
			
        ggplot2::ggsave(paste(outfile,"_before_pruning_and_recoding_ordered_by_pos.pdf",sep=""))
        RegionScan:::RegionHeatmap(cormat1[names(binvector),names(binvector)],
          binvector=binvector, type_mat="Correlation", 
            ptitle="Within region correlation (before pruning & recoding),\n SNPs ordered by LDbin", 
            colvect=colvect, clustorder=F)
        ggplot2::ggsave(paste(outfile,"_before_pruning_and_recoding_ordered_by_LDbin.pdf",sep=""))
        
        RegionScan:::RegionHeatmap(cormat1[tmp$variant,tmp$variant], type_mat="Correlation", 
          ptitle="Within region correlation (after pruning & before recoding),\n SNPs ordered by pos", 
          colvect=colvect, clustorder=F)
         ggplot2::ggsave(paste(outfile,"_after_pruning_and_before_recoding_ordered_by_pos.pdf",sep=""))
       
        RegionScan:::RegionHeatmap(cormat1[names(aliasout$final),names(aliasout$final)], type_mat="Correlation", 
          ptitle="Within region correlation (after pruning & before recoding),\n SNPs ordered by LDbin", 
           binvector=aliasout$final, colvect=colvect, clustorder=F)
         ggplot2::ggsave(paste(outfile,"_after_pruning_and_before_recoding_ordered_by_LDbin.pdf",sep=""))  
       
       RegionScan:::RegionHeatmap(cormat3[tmp$variant,tmp$variant], type_mat="Correlation" , 
          ptitle="Within region correlation (after pruning & recoding),\n SNPs ordered by pos",
          colvect=colvect, clustorder=F)
        ggplot2::ggsave(paste(outfile,"_after_pruning_and_after_recoding_ordered_by_pos.pdf",sep=""))
       
        RegionScan:::RegionHeatmap(cormat3, type_mat="Correlation" , 
          ptitle="Within region correlation (after pruning & recoding),\n SNPs ordered by LDbin",
          binvector=aliasout$final, colvect=colvect, clustorder=F)
        ggplot2::ggsave(paste(outfile,"_after_pruning_and_after_recoding_ordered_by_LDbin.pdf",sep=""))
       }
     
		} else { # 1 SNP only
			binvector<- 1 
			names(binvector) <-processout$SNPinfo[,"variant"]
			mglmtime<-glmout <- RegionScan:::glmfit(data=processout$data, pheno, binvector=binvector,
					covlist=covlist,family=family,tol=tol,firthreg=firthreg)
			MLCBout <- RegionScan:::MLC(beta=glmout$beta_g, sigmainv=glmout$invcov_g,  
				binvector=binvector, codechange=0, tol=tol)
			PC80out<- Waldout<- LCBout <- MLCBout
			SKATout<-RegionScan:::SKATp( data=processout$data, pheno=pheno, 
				covlist=covlist, pheno_type=pheno_type,
				geno_type=geno_type,binvector=binvector,SKAT_kernel=SKAT_kernel, 						
				SKAT_weights=SKAT_weights, SKAT_weights_beta=SKAT_weights_beta)
			
			sgout<-	c(variant=processout$SNPinfo[,"variant"],
				sglm.beta=unname(glmout$beta_g),
						sglm.se=unname(glmout$beta_SE), 
						sglm.pvalue=glmout$p_g)
				
			binout <- c(chr=chr, region=region, start.bp=start, end.bp=end, 
				binstart.bfP.bp=unique(processout$SNPinfo["bp"]), binend.bfP.bp=unique(processout$SNPinfo["bp"]),
				binstart.afP.bp=unique(processout$SNPinfo["bp"]), binend.afP.bp=unique(processout$SNPinfo["bp"]), 
				binsize.bfP=1, binsize.afP=1, MLCBout$deltabin)
   
			snpout<-c(chr=chr, region=region, start.bp=start, end.bp=end, 
				bin=1, processout$SNPinfo[,c("bp","multiallelic","ref","alt","maf")],
				MLC.codechange=0, LC.codechange=0,sgout, 
				mglm.vif="NA", mglm.beta=glmout$beta_g, 
				mglm.se=glmout$beta_SE, mglm.pvalue=glmout$p_g )
				
			regionout<-c(chr=chr, region=region, start.bp=start, end.bp=end, 
				nSNPs=nSNPs, nSNPs.kept=nSNPs.kept, maxVIF="NA",
				Wald=Waldout$stat,Wald.df=Waldout$df, Wald.p=Waldout$pvalue,   
				PC80=PC80out$stat, PC80.df=PC80out$df,PC80.p=PC80out$pvalue,
				MLCB=MLCBout$stat,MLCB.df=MLCBout$df,MLCB.p=MLCBout$pvalue,
				LCB=LCBout$stat,LCB.df=LCBout$df,LCB.p=LCBout$pvalue,
				SKAT.p=SKATout$SKAT.p,SKATO.p=SKATout$SKATO.p,
				simes.p=glmout$p_g , simpleM.df=1 , 
				simpleM.p=glmout$p_g ,GATES.p=glmout$p_g, single_Wald.p=glmout$p_g)
        
		if(isTRUE(alltests)) {
			MLCZout <- LCZout <- RegionScan:::MLC(Z=glmout$Z_g, invcor=glmout$invcor_g,  binvector=binvector, codechange=0, tol=tol)
			regionout <- c(regionout,MLCZ=MLCZout$stat, MLCZ.df=MLCZout$df, MLCZ.p=MLCZout$pvalue,LCZ=LCZout$stat,LCZ.df=LCZout$df,LCZ.p=LCZout$pvalue)	
			binout<- c(binout,MLCZout$deltabin)
 		}		 
	} # ENF IF (regions with 1 SNP)
    if (isTRUE(covout) && length(covlist)>0) { 
      covoutput<-data.frame(cbind(chr=chr, region=region, start.bp=start, end.bp=end, covlist=covlist,
        mglm_beta=glmout$beta_covlist, mglm_pvalue=glmout$p_covlist)) }
        
	} # END SNP INFO EXISTING
 }, error=function(e){cat("ERROR :",	conditionMessage(e), "\n")})
  
   	return(list(regionout=regionout, binout=binout, snpout=snpout, filterout=filterout, 
   outsingleSNPall=outsingleSNPall, covarout=covoutput))
  }  
  #print(rscan)

  if(length(rscan)>6) {
	    regionout<-as.data.frame(do.call("rbind", lapply(rscan[,"regionout"],function(x) { unlist(x)})))
	    binout<-as.data.frame(do.call("rbind", rscan[,"binout"]))
	    snpout<-as.data.frame(do.call("rbind", rscan[,"snpout"]))
	    filterout<-as.data.frame(do.call("rbind", rscan[,"filterout"]))
		outsingleSNPall<-as.data.frame(do.call("rbind",rscan[,"outsingleSNPall"]))
		covarout<-as.data.frame(do.call("rbind", rscan[,"covarout"]))
               
  } else {
		
	    regionout<-as.data.frame(do.call("rbind", lapply(rscan["regionout"],function(x) { unlist(x)})))
	    binout<-as.data.frame(do.call("rbind", rscan["binout"]))
	    snpout<-as.data.frame(do.call("rbind", rscan["snpout"]))
	    filterout<-as.data.frame(do.call("rbind", rscan["filterout"]))
		outsingleSNPall<-as.data.frame(do.call("rbind", rscan["outsingleSNPall"]))
		covarout<-as.data.frame(do.call("rbind", rscan["covarout"]))
  }
	regionout$chr<-as.numeric(as.character(regionout$chr))
	if (isTRUE(parallel)) { parallel::stopCluster(cl = my.cluster) }
	
	rownames(regionout)<-rownames(binout)<-rownames(snpout)<-rownames(filterout)<-rownames(outsingleSNPall)<-rownames(covarout)<-NULL
	
	return( list(regionout=regionout, binout=binout, snpout=snpout, 
    filterout=filterout, outsingleSNPall=outsingleSNPall, covout=covarout) )
} 
