#' recodeVCF : prepare input file for region-level analysis
#'
#' @vcfname : name of the vcf file (without extension)
#' @multiallelic: if TRUE includes multi-allelic variants (default is FALSE)
#' @multial_nmaxalleles: number of maximum alleles per variant kept (default is 2)
#' @qc_machr2: threshold to filter out variants with imputation quality score < qc_machr2
#' @info_score: threshold to filter out variants with imputation quality score < info_score
#' @qcfiler : data.frame including imputation quality scores
#' @chr=chr: chromosome # to extract information from vcf file
#' @start= start position to extract information from vcf file 
#' @end= end position to extract information from vcf file 
#' @export

recodeVCF <- function(vcfname, multiallelic, multial_nmaxalleles, qc_machr2, 
	info_score, qcfiler, chr=chr, start=start, end=end) {

	if (isFALSE(multiallelic)) { try(system(sprintf("vcftools --vcf %s --chr %i --from-bp %i --to-bp %i --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ./tmp_dir/%s", vcfname, chr, start, end, vcfregion), intern = TRUE))
	} else { try(system(sprintf("vcftools --vcf %s --chr %i --from-bp %i --to-bp %i --remove-indels --max-alleles %i --recode --recode-INFO-all --out ./tmp_dir/%s", 
		vcfname, chr, start, end.bp, multial_nmaxalleles, vcfregion), intern = TRUE)) 
	}
		 	  
	vcfpath <-paste("./tmp_dir/",vcfregion,".recode.vcf",sep="")
	pos <- dup.pos.count <- row.multpos<- NULL
	vcfile <- read.vcfR(vcfpath)
	infoI <- data.frame(vcfile@fix, stringsAsFactors =F)
	infoI$variant <-paste("X",infoI$CHROM,".",infoI$POS,".",
	infoI$REF,".",infoI$ALT,sep="")	
	probsI <- t(extract.gt(vcfile,element='GP',IDtoRowNames=F))
	dosageI <- t(extract.gt(vcfile,element='DS',IDtoRowNames=F))
		
  	if (!is.null(qcfiler)) { 
		qcfiler[,"POS"] <- as.character(qcfiler[,"POS"])
		qcfiler[,"CHROM"] <- as.character(qcfiler[,"CHROM"])
		qcfiler[,"REF"] <- as.character(qcfiler[,"REF"])
		qcfiler[,"ALT"] <- as.character(qcfiler[,"ALT"])
		infoI <- merge(infoI[,-8],qcfiler, by.x=c("CHROM","POS","REF","ALT"), by.y=c("CHROM","POS","REF","ALT" ) )
		infoI <-infoI[,c("CHROM","POS","ID","REF","ALT", "QUAL","FILTER","INFO","variant")]
	}	
		
	if (nrow(na.omit(dosageI))==0) { #IF no dosage file
		ones <- apply(probsI,2,function (x) { matrix(as.numeric(unlist(strsplit(x,split=","))),	ncol=3,byrow=T)[,2]})
		twos <- apply(probsI,2,function (x) { matrix(as.numeric(unlist(strsplit(x,split=","))),	ncol=3,byrow=T)[,3]})
		dosageI <- ones+2*twos
		rownames(dosageI) <- rownames(probsI)
	} 
	colnames(dosageI)<-colnames(probsI)<-as.character(infoI$variant)
		
   # identification of indels and multiallelic variants
    indels <- unique(sort(c(which(nchar(infoI$REF)>1),which(nchar(infoI$ALT)>1)))) 
	if (length(indels)>0) { infoI<-infoI[-indels,,drop=F] }
	dupos <- infoI[duplicated(infoI$POS),"POS"] 
	dupos.nb <- unlist(lapply(unique(dupos), function(x) {length(which(infoI$POS==x))})) 

 	if (isFALSE(multiallelic)) {
		if (length(dupos)>0) { 
          infoI <- infoI[-which(infoI[duplicated(infoI$POS),"POS"]%in%infoI$POS),] 
          dosageI <- dosageI[,colnames(dosageI)%in%infoI$variant, drop=F]
		}
    } else { 
		if (length(dupos)>0){
			names(dupos.nb) <- unique(dupos)
			rem <- which(dupos.nb>multial_nmaxalleles) # remove SNPs with more than multial_nmaxalleles 
			if(length(rem)>0) { infoI<- infoI[-which(infoI$POS%in%unique(names(rem))),] }
		}
		dosageI <- dosageI[,colnames(dosageI)%in%infoI$variant, drop=F]
		probsI <- probsI[,colnames(probsI)%in%infoI$variant, drop=F]
		pos <- unique(infoI[duplicated(infoI$POS),"POS"], drop=F)
				
		infoI$multialSNP <- 0
		infoI$multialSNP.rec <- 0
		infoI$old.ref<-infoI$REF
		infoI$old.alt<-infoI$ALT
				
		if(length(pos)>0) {
		  for (var in 1:length(pos)) {
			alt <- which(infoI[,"POS"]%in%pos[var])
			infoI[alt,"multialSNP"] <- 1
			num.alt <- length(alt)
			dosageI.mv <- dosageI[,alt, drop=F]
			probsI.mv <- probsI[,alt, drop=F]
			if ( is.vector(dosageI.mv)==T ) { 
				freq_alt <- mean(as.numeric(as.character(dosageI.mv)))/2 
			} else {
				freq_alt <- apply(dosageI.mv,2, function(x) { mean(as.numeric(as.character(x)))/2 })
			}		
			freq <- c(1-sum(freq_alt),freq_alt)
			names(freq) <- c("ref",paste("a",1:num.alt,sep=""))
			v <- freq[which.max(freq)]
							
			if(names(v)!="ref") { #flip
				ones <- apply(probsI.mv,2, function (x) { 
							  matrix(as.numeric(unlist(strsplit(x,split=","))), ncol=3,byrow=T)[,2]})
				twos <- apply(probsI.mv,2,
						  function (x) { 
							matrix(as.numeric(unlist(strsplit(x,split=","))),ncol=3,byrow=T)[,3]})
							dosageI[,alt[which.max(freq_alt)],drop=F] <- 2*(1-rowSums(twos))-rowSums(ones)
							newref <- infoI[alt[which.max(freq_alt)],"old.alt",drop=F]
							newalt <- infoI[alt[which.max(freq_alt)],"old.ref",drop=F]
							infoI[alt[which.max(freq_alt)],"ALT",drop=F] <- newalt
							infoI[alt,"REF",drop=F] <- newref
							infoI[alt,"multialSNP.rec",drop=F] <- 1
			}
		}
	  } 
	}

	dosageF<- cbind.data.frame(ID=rownames(dosageI), apply(dosageI,2,function(x){as.numeric(as.character(x))}))
    freq.alt <- apply(dosageF[,-1], 2, mean)/2
	R2_mach.alt <- apply(dosageF[,-1], 2, function(x) { 	
					z <- var(x)/(2*(mean(x)/2)*(1-(mean(x)/2))) 
					z[which(is.nan(z))] <- 0
					return(z)
			})
    colnames(infoI)[c(1,2,4,5)] <- c("chr","pos","ref","alt")
	infoF <- cbind(infoI, freq.alt, R2_mach.alt)
		 	
   if(isFALSE(multiallelic)) {
        infoF$multialSNP<-infoF$multialSNP.rec<-0
		infoF <- infoF[,c("chr","pos","variant", "ref", "alt",
          "freq.alt", "R2_mach.alt","INFO"),drop=F]
   } else { 		
        infoF <- infoF[,c("chr", "pos", "variant", "ref", "alt","freq.alt", 
          "R2_mach.alt", "multialSNP", "multialSNP.rec", "ID", "QUAL", 
          "FILTER", "INFO", "old.ref", "old.alt"), drop=F]
   }
    
    if(!is.null(info_score)) { 
      infoF<-infoF[infoF$INFO>=info_score,, drop=F]
      dosageF <- dosageF[,colnames(dosageF)%in%c("ID", as.character(infoF$variant)), drop=F]
    }
    
    if(!is.null(qc_machr2)) { 
      infoF<-infoI[infoF$R2_mach.alt>=qc_machr2,, drop=F]
      dosageF <- dosageF[,colnames(dosageF)%in%c("ID", as.character(infoF$variant)),drop=F]
    }
     
    if(!is.null(region)) { infoF$region<-as.character(region) }
   	if(!is.null(start.bp)) { infoF$start.bp <- start.bp }
	if(!is.null(end.bp)) { infoF$end.bp <- end.bp }
 	 
    if(isFALSE(multiallelic)) {
		message("Total No. of subjects in VCF file: ", nrow(vcfile@fix))
		message("Total No. of variants in original VCF file: ", nrow(dosageI))
		message("Total No. of biallelic SNPs kept in VCF file: ", length(infoF$pos))
     } else { 	
		message("Total No. of subjects in VCF file: ", nrow(vcfile@fix))
		message("Total No. of variants in original VCF file: ",length(infoF$POS))
		message("Total No. of multiallelic SNPs kept in VCF region file: ", length(which(infoF$multialSNP==1)))
		message("Total No. of biallelic SNPs kept in VCF region file: ", length(which(infoF$multialSNP==0)))
   	}   
	return(list(dosage=dosageF, info=infoF))
}   