#' LocusPlot : produces a pdf file including a plot of region-level results for a set of contiguous regions
#'
#' @chr: chromosome number of the locus to plot
#' @regionlist: list of contiguous regions to plot
#' @pheno: name of the phenotype 
#' @regscanout: output object from regionscan 
#' @p_region: region-level significance level (default is 5.86e-07)
#' @p_single: significance-level for single-SNP tests (default is  5e-8)
#' @real_size: based on physical size of the regions (default is TRUE)
#' @allSNPs: plot single-SNP results for pruned SNPs as well  (default is FALSE)
#' @region_tests: specify the list of region-level tests to plot (default is NULL, all tests will be plotted)
#' @outname: name of the output file produced 
#' @export

LocusPlot <- function(chr, regionlist ,pheno=pheno, regscanout , p_region = 5.86e-07,  p_single = 5e-8, 
  real_size=TRUE, allSNPs=FALSE, region_tests=NULL, outname=NULL) {
  
  region_results<-subset(regscanout$regionout, (chr==chr & region %in% regionlist))
  region_results$chr <- as.numeric(region_results$chr)
  
  if(isTRUE(allSNPs)) {
      var_results<-subset(regscanout$outsingleSNPall, (chr==chr & region %in% regionlist))
  } else {
      var_results<-subset(regscanout$snpout, (chr==chr & region %in% regionlist)) 
  }
  var_results$chr <- as.numeric(var_results$chr)
  var_results$bp <- as.numeric(var_results$bp)
    
  if(is.null(region_tests)) {  region_tests<- colnames(region_results)[grep(".p",colnames(region_results),fixed=T)] }
  region_union <- data.frame(do.call('rbind', lapply(region_tests, function(test) {
    tmp<-data.frame(cbind(chr=chr, region=regionlist, region_results[,c("start.bp","end.bp", test, "Wald.df")], test=test))
    colnames(tmp)<-c("chr","region","start.bp","end.bp", "pv", "no.snps", "test")
      return(tmp)
   })))
  region_union$start.bp<-as.integer(region_union$start.bp)
  region_union$end.bp<-as.integer(region_union$end.bp)
  region_union$pos <- (region_union$start.bp+region_union$end.bp)/2
  region_union$pv <- -log10(as.numeric(as.character(region_union$pv)))
  
  sg_snp_union<-subset(var_results, (chr==chr & var_results$region %in% regionlist))[,c("chr","region","start.bp","end.bp","sglm.pvalue","bp")]
  colnames(sg_snp_union)[1:5]<-colnames(region_union)[1:5]
  sg_snp_union$region <- as.integer(sg_snp_union$region)
  sg_snp_union$start.bp <- as.integer(sg_snp_union$start.bp)
  sg_snp_union$end.bp <- as.integer(sg_snp_union$end.bp)
  sg_snp_union$bp <- as.integer(sg_snp_union$bp)
  sg_snp_union$pv <- -log10(as.numeric(as.character(sg_snp_union$pv)))
  sg_snp_union$type <- "single SNP"
  #  pv=-log10(var_results$sg.pval[single_indi]) )

  sg_snp_union_minp <- data.frame()
  for (i in regionlist){
    n <- sg_snp_union [sg_snp_union$region==i,]
    sg_snp_union_minp <- rbind(sg_snp_union_minp, n[which.max(n$pv),])
  }

  sg_snp_union$region_pos <- 0
  for (i in regionlist){
    sg_snp_union$region_pos[sg_snp_union$region==i] <- 
      i-0.5+(sg_snp_union$bp[sg_snp_union$region==i]-unique(region_union$start.bp[region_union$region==i]))/(unique(region_union$end.bp[region_union$region==i])-unique(region_union$start.bp[region_union$region==i]))  
  }
  
  sg_snp_union_minp$region_pos <- 0
  for (i in regionlist){
    sg_snp_union_minp$region_pos[sg_snp_union_minp$region==i] <- 
      i-0.5+(sg_snp_union_minp$bp[sg_snp_union_minp$region==i]-region_union$start.bp[region_union$region==i][1])/(region_union$end.bp[region_union$region==i][1]-region_union$start.bp[region_union$region==i][1])  
  }
  
  sg_snp_union_minp$pos <- region_union$pos[1:length(regionlist)]
  sg_snp_union_minp$type <- "single SNP minP"
  #sg_snp_union$region_pos<-as.integer(subset(var_results, (chr==chr & var_results$region %in% regionlist))[,"bp"])

  #####  plot region and all single SNPs in these regions
  region_union$pv<-as.numeric(as.character(region_union$pv))
  y_limit <- max(region_union$pv) + 5
  region_union$region<-as.integer(region_union$region)
  
  if(isTRUE(real_size)) {
    ggplot2::ggplot() +
    ggplot2::geom_rect(data=NULL,aes(xmin=unique(as.numeric(region_union$start.bp[seq(1,nrow(region_union),2)])),
    xmax=unique(as.numeric(region_union$end.bp[seq(1,nrow(region_union),2)])),ymin=0.01,ymax=Inf),  fill="dark grey",alpha=0.6) +
    ggplot2::geom_segment(data=region_union,linewidth=3,aes(x=start.bp, xend=end.bp, y=pv, yend=pv,color=test)) +
    ggplot2::geom_point(data=sg_snp_union,aes(x=bp, y=pv,shape=type,size=type),color="#00BFC4")+
    ggplot2::geom_hline(yintercept = -log10(p_region), color = "#F8766D", linetype = "dashed")+
    ggplot2::geom_hline(yintercept = -log10(p_single), color = "#00BFC4", linetype = "dashed")+
    ggplot2::scale_x_continuous(name=paste("chr",chr,"(bp)", sep=""))+
    ggplot2::scale_y_continuous(name="-log10(P-value)",trans='sqrt',limits=c(0.01,y_limit),expand = c(0, 0))+
    #scale_y_continuous(name="-log10(p-v)",limits=c(0.5,120))+
    ggplot2::ggtitle(paste0(pheno," region-level and single SNP P-values")) +
    #scale_color_manual(name = "",values = cols)+ 
    ggplot2::scale_shape_manual(name = "",values = c(18, 1))+ 
    ggplot2::scale_size_manual(name = "",values = c(2.2, 0.6))+ 
    ggplot2::scale_alpha_manual(name = "",values = c(1, 0.6))+
    ggplot2::geom_text(region_union,mapping=aes(x=pos, y=0.1,label=region),size=2.5)+
    ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
    #geom_line(data=sg_snp_union_minp,aes(x=bp, y=pv,color=type),alpha = 1)+
    #geom_line(data=region_union,aes(x=pos, y=pv,color=test),alpha = 1)+   
  } else { # equal size
  
  xmin_<-unique(as.numeric(unique(region_union$region)[seq(1,length(unique(region_union$region)),2)]))
  
  ggplot2::ggplot() +  
    ggplot2::geom_rect(data=NULL,aes(xmin=xmin_-0.5, xmax=xmin_+0.5,ymin=0.01,ymax=Inf),  fill="dark grey",alpha=0.6)+
    ggplot2::geom_segment(data=region_union,linewidth=3,aes(x=region-0.5, xend=region+0.5, y=pv, yend=pv, color=test)) +
    ggplot2::geom_point(data=sg_snp_union,aes(x=region_pos, y=pv,shape=type,size=type), color="#00BFC4")+
    ggplot2::geom_line(data=sg_snp_union_minp,aes(x=region_pos, y=pv,color=type), alpha = 1)+
    ggplot2::geom_line(data=region_union,aes(x=region, y=pv,color=test),alpha = 1)+
    ggplot2::geom_hline(yintercept = -log10(p_region), color = "#F8766D", linetype = "dashed")+
    ggplot2::geom_hline(yintercept = -log10(p_single), color = "#00BFC4", linetype = "dashed")+
    ggplot2::scale_x_continuous(name=paste("chr",chr,": regions #",sep=""), breaks = region_union$region)+
    ggplot2::scale_y_continuous(name="-log10(P-value)",trans='sqrt',limits=c(0.01,y_limit),expand = c(0, 0))+
    #scale_y_continuous(name="-log10(p-v)",limits=c(0.5,120))+
    ggplot2::ggtitle(paste0(pheno," region-level and single SNP P-values")) +
    #scale_color_manual(name = "",values = cols)+ 
    ggplot2::scale_shape_manual(name = "",values = c(18, 1))+ 
    ggplot2::scale_size_manual(name = "",values = c(2.2, 0.6))+ 
    ggplot2::scale_alpha_manual(name = "",values = c(1, 0.6))+
    ggplot2::geom_text(region_union,mapping=aes(x=region, y=pv,label=no.snps),size=2.5)+
    ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
  }
  
  ggplot2::ggsave(filename=paste(outname,".pdf",sep=""), width=25, height=15, units="cm", dpi=300)
}
