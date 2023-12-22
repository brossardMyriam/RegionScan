#' qqregscanPlot : qq plot for a specified test
#'
#' @regionout: output from regscan
#' @outname: name of the output file
#' @test_name: name of the test to be plotted
#' @export

qqregscanPlot  <- function( regionout , outname, test_name ) {  
  pval<-as.numeric(as.character(regionout[,test_name]))
  stat <- qchisq(pval, df=1, lower.tail=F) 
  gc <- round(median(stat)/qchisq(p=0.5,df=1),3)
  np <- length(pval)
  data_p <- data.frame(log.p = sort(-log10(pval)), EP = -log10(np:1/(np+1)))
  y_lim <- max(apply(data_p,2,max)) + 5


  qqp <- ggplot(data_p, aes(EP)) +
    geom_point(aes(y =log.p)) + 
    geom_abline(slope = 1, intercept = 0,color="darkgrey")+
    xlab("-log10(Expected)")+
    ylab("-log10(Observed)")+
    xlim(0,ceiling(max(-log10(np:1/(np+1)))))+
    ylim(0,y_lim)+
    ggtitle(paste("QQ plot of ", test_name,sep=""))+
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(0,0), 
          legend.position=c(0.05, 0.65),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=15))
    qqp + guides(color=guide_legend(""))
    ggsave(paste(outname, ".pdf",sep=""))
    message(paste("GC=",gc,sep=""))
}