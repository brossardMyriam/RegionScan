#' MLCbinsnpPlot : produces a pdf file including the position of the SNPs along the chr per LD bin for the region specified
#'
#' @rscanout: chromosome number of the locus to plot
#' @chr_: list of contiguous regions to plot
#' @region_: name of the phenotype 
#' @export

MLCbinsnpPlot <- function( rscanout, chr_, region_) 
{
  snpout<-rscanout$snpout[,c("chr","region", "bin","variant","bp")]
  filterout<-rscanout$filterout[,c("chr","region","bin","variant","bp")]
  binout<-rscanout$binout[,c("chr", "region","start.bp","end.bp","bin","binstart.afP.bp", "binend.afP.bp", 
    "binstart.bfP.bp", "binend.bfP.bp","binsize.bfP", "binsize.afP")]
  snpout<-subset(snpout, chr==chr_ & region==region_)
  filterout<-na.omit(subset(filterout, chr==chr_ & region==region_))
  binout<-subset(binout, chr==chr_ & region==region_)
  allout<- unique(cbind(rbind(snpout,filterout),kept=c(rep(1,nrow(snpout)),rep(0,nrow(filterout)))))
  nbins <- max(as.numeric(as.character(snpout$bin)))
  colvect <- scales::hue_pal()(nbins)  
       
  binout<-binout[order(as.numeric(as.character(binout$bin))),]
  binout$binstart.afP.bp <-as.numeric(as.character(binout$binstart.afP.bp )) 
  binout$binend.afP.bp <-as.numeric(as.character(binout$binend.afP.bp ))
  binout$binstart.bfP.bp <-as.numeric(as.character(binout$binstart.bfP.bp )) 
  binout$binend.bfP.bp <-as.numeric(as.character(binout$binend.bfP.bp ))
  
  binout$bin<-factor(binout$bin,level=unique(binout$bin),order=TRUE)
  nSNP.bfP<-sum(as.numeric(as.character(binout$binsize.bfP)))
  nSNP.afP<-sum(as.numeric(as.character(binout$binsize.afP)))
  snpout$bin<-factor(snpout$bin,level=unique(snpout$bin),order=TRUE)
  snpout$bp<-as.numeric(as.character(snpout$bp))
  
  outname<-paste("MLC_LDbin_SNPpos_chr",chr_,"_region",region_,sep="")
  ggplot2::ggplot() +  ggplot2::geom_segment(data=binout, aes(x=binstart.afP.bp, y=bin, xend=binend.afP.bp, yend=bin, color=bin)) + 
  ggplot2::geom_point(data=snpout, aes(x=bp, y=bin, color=bin))  +
  ggplot2::ggtitle(paste("region: ", region_, ", after pruning : ", nSNP.afP, " SNP(s)",sep="")) +  
  ggplot2::xlab("pos") + ggplot2::ylab("MLC LD bin#") +
  ggplot2::scale_color_manual(values=colvect, labels=paste("Bin", binout$bin, ": ", binout$binsize.afP , " SNP(s)", sep=""))
  ggplot2::ggsave(paste(outname, "_after_pruning.pdf",sep=""),width = 15, height = 10,  units ="cm")
  
  ggplot2::ggplot() + ggplot2::geom_segment(data=binout, aes(x=binstart.bfP.bp, y=bin, xend=binend.bfP.bp, yend=bin, color=bin)) + 
  ggplot2::geom_point(data=subset(allout, kept==0 ), aes(x=bp, y=bin),color="darkgrey") +
  ggplot2::geom_point(data=subset(allout, kept==1 ), aes(x=bp, y=bin, color=bin)) +
  ggplot2::ggtitle(paste("region", region_, ", after (before) pruning : ", nSNP.afP, " (", nSNP.bfP,") ", " SNPs",sep="")) +  
  ggplot2::xlab("pos") + ggplot2::ylab("MLC LD bin#")+
  ggplot2::scale_color_manual(values=colvect, labels=paste("Bin", binout$bin, ": ", binout$binsize.afP, " (", binout$binsize.bfP, ") SNP(s)", sep=""))
  ggplot2::ggsave(paste(outname, "_before_and_after_pruning.pdf",sep=""),width = 15, height = 10,  units ="cm")
  
}