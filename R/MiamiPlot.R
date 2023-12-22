#' MiamiPlot : produces a pdf file including a miami plot for a pair of tests specified 
#'
#' @regscanout : regionscan output object
#' @p_region: significance threshold for region-level tests (default is FALSE)
#' @p_single: significance threshold for single-SNP tests  (default is 2)
#' @test1: specify which test (from the regionscan output) will be plotted on the top panel; must be  Wald.p, MLCB.p, PC80.p, SKAT.p,  SKATO.p, GATES.p,  SimpleM.p or sglm.pvalue
#' @test2 : specify which test (from the regionscan output) will be plotted on the bottom panel ; must be  Wald.p, MLCB.p, PC80.p, SKAT.p,  SKATO.p, GATES.p,  SimpleM.p or sglm.pvalue
#' @allSNPs: if TRUE, includes all single-SNP results (including SNPs pruned on LD; default is FALSE)
#' @outname=name for the output file (without extension) 
#' @chr= chromosome name (to zoom on a particular region), optional
#' @zoom= region to zoom (e.g. chr22:100000:200000)
#' @test1_highlight: region to highlight (e.g. chr22:100000:200000)
#' @test2_highlight: region to highlight (e.g. chr22:100000:200000)
#' @export


MiamiPlot  <- function( regscanout , p_region = 5.86e-07,  p_single = 5e-8, 
   allSNPs=FALSE, test1=NULL, test2=NULL, outname=NULL, chr=NULL, zoom=zoom, 
   test1_highlight=NULL,test2_highlight=NULL) {  
           
  regionout<-regscanout$regionout
  if(isTRUE(allSNPs)) {
     snpout<-regscanout$outsingleSNPall
  } else {
    snpout<-regscanout$snpout 
  }
  
  #snp_GRange <- toGRanges(data.frame(chr=paste("chr",snpout$chr,sep=""), start=snpout$bp, end=snpout$bp, pval=snpout$sg.pval))
  regionout$chr_region <- paste("chr",regionout$chr,"_", regionout$region,sep="")
  regionout <-regionout[order(regionout$chr,regionout$start.bp),]
  
  RegionScan_prep <- do.call('rbind', by(regionout, 
  	list(regionout$chr_region), 
  	function(x) {
  		chr_region<-x$chr_region
  		chr<-x$chr
  		start<-x$start.bp
  		end<-x$end.bp
  		Wald.p<-x$Wald.p
  		MLCB.p<-x$MLCB.p
  		PC80.p<-x$PC80.p
  		SKAT.p<-x$SKAT.p
  		SKATO.p<-x$SKATO.p
  		GATES.p<-x$GATES.p
  		SimpleM.p<-x$SimpleM.p
  		y<-snpout[which(snpout$chr==chr & snpout$bp>=start & snpout$bp<=end),]
  		y$Wald.p<-Wald.p
  		y$MLCB.p<-MLCB.p
  		y$PC80.p<-PC80.p
  		y$SKAT.p<-SKAT.p
  		y$SKATO.p<-SKATO.p
  		y$GATES.p<-GATES.p
  		y$SimpleM.p<-SimpleM.p	
  		y$start.bp<-start
  		y$end.bp<-end
  		y$chr_region<-chr_region
  		return(y)
  	})) 
  
  RegionScan_prep2 <-RegionScan_prep[order(RegionScan_prep$chr,RegionScan_prep$region),]
  RegionScan_GRange <- toGRanges(data.frame(chr=paste("chr",RegionScan_prep2$chr,sep=""), start=RegionScan_prep2$bp, end=RegionScan_prep2$bp))
  maxy<-max(-log10(na.omit(as.numeric(as.character(c(RegionScan_prep2[,test1],RegionScan_prep2[,test2]))))))+5
  region_tests<-c("MLCB.p","MLCZ.p","LCB.p","LCZ.p","Wald.p","PC80.p","SKAT.p","SKATO.p","MinPJ.p","GATES.p","simpleM.p")
  
  if(isTRUE(test1%in%region_tests) & isTRUE(test2%in%region_tests)) {  
    pdf(paste(outname, ".pdf",sep=""),height=7, width=14)
     if(!is.null(chr)) {  kp<-plotKaryotype("hg19",plot.type=4, cytobands = GRanges(),	cex=0.9, labels.plotter = NULL, chromosomes=paste("chr",chr,sep=""))} 
     if(!is.null(zoom)) { kp<-plotKaryotype("hg19",plot.type=4, cytobands = GRanges(),	cex=0.9, labels.plotter = NULL, main=zoom, zoom=zoom)} 
  
      kp <- karyoploteR::kpPlotManhattan(kp,ymax=maxy, data=RegionScan_GRange, pval=as.numeric(as.character(RegionScan_prep2[,test1])),
    		points.col="2blues", r0=0.5, r1=1,highlight=test1_highlight_region, highlight.col="orange",
    		genomewideline=-log10(p_region),genomewide.lwd=2,genomewide.lty=2,
    		suggestiveline=-log10(p_region),suggestive.lwd=2,suggestive.lty=2,points.cex=1.5)
     karyoploteR::kpAxis(kp, ymin=0, ymax=maxy, r0=0.5, r1=1, tick.pos = seq(from=5,to=maxy,by=10), cex=1.2)
     karyoploteR::kpAddLabels(kp, labels = paste(test1,"(-log10)"), srt=90, pos=3, r0=0.5, r1=1, cex=1.2, label.margin = 0.035, side="left")
       
      kp <- karyoploteR::kpPlotManhattan(kp,ymax=maxy, data=RegionScan_GRange, pval=as.numeric(as.character(RegionScan_prep2[,test2])),
    		points.col="2grays", r0=0.5, r1=0,highlight=test2_highlight_region, highlight.col="orange",
    		genomewideline=-log10(p_region),genomewide.lwd=2,genomewide.lty=2,
    		suggestiveline=-log10(p_region),suggestive.lwd=2,suggestive.lty=2, points.cex=1.5)
      karyoploteR::kpAxis(kp, ymin=0, ymax=maxy, r0=0.5, r1=0, tick.pos = seq(from=5,to=maxy,by=10), cex=1.2)
      karyoploteR::kpAddLabels(kp, labels = paste(test2,"(-log10)"), srt=90, pos=3, r0=0, r1=0.5, cex=1.2, label.margin = 0.035, side="left")
      if(!is.null(zoom)) {  karyoploteR::kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = 1e5, cex=1.2) }
      dev.off()
    } 
   
    if(isTRUE(test1=="sglm.pvalue")) { 
        pdf(paste(outname, ".pdf",sep=""),height=7, width=14)
        if(!is.null(chr)) {  kp<-karyoploteR::plotKaryotype("hg19",plot.type=4, cytobands = GRanges(),	cex=0.9, labels.plotter = NULL, chromosomes=paste("chr",chr,sep=""))} 
        if(!is.null(zoom)) { kp<-karyoploteR::plotKaryotype("hg19",plot.type=4, cytobands = GRanges(),	cex=0.9, labels.plotter = NULL, main=zoom, zoom=zoom)} 
        kp <- karyoploteR::kpPlotManhattan(kp,ymax=maxy, data=RegionScan_GRange, pval=as.numeric(as.character(RegionScan_prep2[,test1])),
      		points.col="2blues", r0=0.5, r1=1,highlight=test1_highlight_region, highlight.col="orange",
      		genomewideline=-log10(p_single),genomewide.lwd=2,genomewide.lty=2,
      		suggestiveline=-log10(p_single),suggestive.lwd=2,suggestive.lty=2,points.cex=1.5)
        karyoploteR::kpAxis(kp, ymin=0, ymax=maxy, r0=0.5, r1=1, tick.pos = seq(from=5,to=maxy,by=10), cex=1.2)
       kpAddLabels(kp, labels = paste(test1,"(-log10)"), srt=90, pos=3, r0=0.5, r1=1, cex=1.2, label.margin = 0.035, side="left")
         
        kp <- karyoploteR::kpPlotManhattan(kp,ymax=maxy, data=RegionScan_GRange, pval=as.numeric(as.character(RegionScan_prep2[,test2])),
      		points.col="2grays", r0=0.5, r1=0,highlight=test2_highlight_region, highlight.col="orange",
      		genomewideline=-log10(p_region),genomewide.lwd=2,genomewide.lty=2,
      		suggestiveline=-log10(p_region),suggestive.lwd=2,suggestive.lty=2, points.cex=1.5)
        karyoploteR::kpAxis(kp, ymin=0, ymax=maxy, r0=0.5, r1=0, tick.pos = seq(from=5,to=maxy,by=10), cex=1.2)
        karyoploteR::kpAddLabels(kp, labels = paste(test2,"(-log10)"), srt=90, pos=3, r0=0, r1=0.5, cex=1.2, label.margin = 0.035, side="left")
        if(!is.null(zoom)) {  karyoploteR::kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = 1e5, cex=1.2) }
        dev.off()
    } 
    
      if(isTRUE(test2=="sglm.pvalue")) { 
        pdf(paste(outname, ".pdf",sep=""),height=7, width=14)
        if(!is.null(chr)) {  kp<-karyoploteR::plotKaryotype("hg19",plot.type=4, cytobands = GRanges(),	cex=0.9, labels.plotter = NULL, chromosomes=paste("chr",chr,sep=""))} 
        if(!is.null(zoom)) { kp<-karyoploteR::plotKaryotype("hg19",plot.type=4, cytobands = GRanges(),	cex=0.9, labels.plotter = NULL, main=zoom, zoom=zoom)} 
    
        kp <- karyoploteR::kpPlotManhattan(kp,ymax=maxy, data=RegionScan_GRange, pval=as.numeric(as.character(RegionScan_prep2[,test1])),
      		points.col="2blues", r0=0.5, r1=1,highlight=test1_highlight, highlight.col="orange",
      		genomewideline=-log10(p_region),genomewide.lwd=2,genomewide.lty=2,
      		suggestiveline=-log10(p_region),suggestive.lwd=2,suggestive.lty=2,points.cex=1.5)
        karyoploteR::kpAxis(kp, ymin=0, ymax=maxy, r0=0.5, r1=1, tick.pos = seq(from=5,to=maxy,by=10), cex=1.2)
        karyoploteR::kpAddLabels(kp, labels = paste(test1,"(-log10)"), srt=90, pos=3, r0=0.5, r1=1, cex=1.2, label.margin = 0.035, side="left")
         
        kp <- karyoploteR::kpPlotManhattan(kp,ymax=maxy, data=RegionScan_GRange, pval=as.numeric(as.character(RegionScan_prep2[,test2])),
      		points.col="2grays", r0=0.5, r1=0,highlight=test2_highlight, highlight.col="orange",
      		genomewideline=-log10(p_single),genomewide.lwd=2,genomewide.lty=2,
      		suggestiveline=-log10(p_single),suggestive.lwd=2,suggestive.lty=2, points.cex=1.5)
        karyoploteR::kpAxis(kp, ymin=0, ymax=maxy, r0=0.5, r1=0, tick.pos = seq(from=5,to=maxy,by=10), cex=1.2)
        karyoploteR::kpAddLabels(kp, labels = paste(test2,"(-log10)"), srt=90, pos=3, r0=0, r1=0.5, cex=1.2, label.margin = 0.035, side="left")
        if(!is.null(zoom)) {  karyoploteR::kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = 1e5, cex=1.2) }
        dev.off()
    } 
}
