#' RegionHeatmap : produces heatmaps of correlation matrices within  each region (internal function)
RegionHeatmap <- function(inputmat, binvector=NULL, type_mat=type_mat, 
  ptitle, colvect=colvect, clustorder=clustorder, lim=c(-1,1) )
{
	input <- inputmat
	input <- get_upper_tri(input)
	melted_cormat <- melt(input, na.rm = TRUE)
	if(nrow(input)<=10) { txt=TRUE ; size.txt=3 ; size.X=10 ; size.Y=10}
	if(nrow(input)>10 & nrow(input)<=20 ) { txt=TRUE ; size.txt=2 ; size.X=10 ; size.Y=10  }
	if(nrow(input)>20 & nrow(input)<=30 ) { txt=TRUE ; size.txt=1.8 ; size.X=10 ; size.Y=10  }
	if(nrow(input)>30 & nrow(input)<=40 ) { txt=TRUE ; size.txt=1.5 ; size.X=10 ; size.Y=10  }
	#if(nrow(input)>40 & nrow(input)<=60 ) { txt=TRUE ; size.txt=1.5 ; size.X=8 ; size.Y=8  }
	if(nrow(input)>40 ) {  txt=FALSE ;size.txt=0 ; size.X=0; size.Y=0}
	names(colvect) <- colvect
	
	if (clustorder==FALSE ) {
	hggplot <- ggplot2::ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
		geom_tile(color = "white")+
		scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
		midpoint = 0, limit = c(-1,1), space = "Lab", 
		name=type_mat) +
		theme_minimal()+ 
		theme(axis.text.x = element_text(angle = 90, vjust = 1, 
		size = size.X, hjust = 1)) +
		theme(axis.text.y=element_text(size=size.Y)) +
		coord_fixed() +
		geom_text(aes(Var2, Var1, label = value), 
		color = c("black"), size = size.txt) +
		ggtitle(ptitle) +
		theme(
		  axis.title.x = element_blank(),
		  axis.title.y = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.border = element_blank(),
		  panel.background = element_blank(),
		  axis.ticks = element_blank(),
		  legend.justification = c(1, 0),
		  legend.position = c(0.6, 0.7),
		  legend.direction = "horizontal")
	} else {	  
		
		hggplot <- ggplot2::ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
		ggplot2::geom_tile(color = "white")+
		ggplot2::scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
		midpoint = sum(lim)/2, limit = lim, space = "Lab", 
		name=type_mat) +
		ggplot2::theme_minimal()+ 
		ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 1, 
		size = size.X, hjust = 1)) +
		ggplot2::theme(axis.text.y=element_text(size=size.Y)) +
		ggplot2::coord_fixed() +
		ggplot2::ggtitle(ptitle) +
		ggplot2::theme(
		  axis.title.x = element_blank(),
		  axis.title.y = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.border = element_blank(),
		  panel.background = element_blank(),
		  axis.ticks = element_blank(),
		  legend.justification = c(1, 0),
		  legend.position = c(0.6, 0.7),
		  legend.direction = "horizontal")
	}
		
	if(!is.null(binvector)) {
    mat<-data.frame(cbind(variant=names(binvector),bin=binvector))
		clst <- summary(as.factor(as.character(mat[,"bin"])))
		if (length(clst)>1) {
		  clust <- data.frame(cbind(names(clst),clst))
		  clust2 <-  apply(clust,2,function(x) { as.numeric(as.character(x)) })
		  clust3 <- data.frame(clust2[order(clust2[,1]),]) 
		  colnames(clust3) <- c("bin","nSNPs")
		  clust4<-do.call(rbind, lapply(unique(mat[,"bin"]),function(x) { 
		  clust3[which(clust3[,"bin"]==x),] 
		  }))
	
    cum <- NULL
		
		for(i in 1:length(clust4[,2])) {
			c <- sum(clust4[1:i,2])
			cum <- rbind(cum,c)
			}
		clust4$xend <- cum
		clust4$x <- clust4$xend-clust4$nSNPs
		clust4$y <- clust4$x
		clust4$yend <- clust4$y
		horiz <- clust4
		vertic <- clust4
		vertic$x <- clust4$xend 
		vertic$y <- clust4$yend
		vertic$xend <- vertic$x 
		vertic$yend <- vertic$xend
		seg <- rbind(horiz,vertic)
		
		hggplot2 <- hggplot  +
		ggplot2::geom_text(aes(Var2, Var1, label = value), 
		color = c("black"), size = size.txt) +
		ggplot2::geom_segment(data=seg, aes(x=x+0.5,xend=xend+0.5, 
		y=y+0.5,yend=yend+0.5), size=1, inherit.aes=F)	
		
		seg2 <- horiz
		seg2$yend <- seg2$y <- 0
		seg3 <- seg2
		seg3$x <- seg2$y
		seg3$xend <- seg2$yend		
		seg3$y <- seg2$x
		seg3$yend <- seg2$xend			
		
		names(colvect) <- colvect
		hggplot3 <- hggplot2 + ggplot2::geom_segment(data=seg2, aes(x=x+1,xend=xend+1, 
		  y=y-5,yend=yend-5 ,color=names(colvect)) , 
      size=10, inherit.aes=F,show.legend=F)	+ 
		ggplot2::scale_color_manual(values=colvect)	+
		ggplot2::geom_segment(data=seg3, aes(x=x-5,xend=xend-5, 
      y=y,yend=yend,color=names(colvect)) , 
      size=10, inherit.aes=F,show.legend=F) 
			
		# hggplot3 <- hggplot2 + geom_segment(data=seg2, aes(x=x+1,xend=xend+1, 
		# y=y-5,yend=yend-5 ) , size=10, inherit.aes=F,show.legend=F)	+ 
		# geom_segment(data=seg3, aes(x=x-5,xend=xend-5, 
		# y=y,yend=yend) , size=10, inherit.aes=F,show.legend=F) +
		# scale_color_manual(values=list(colvect))
		
		return(hggplot3) 
		} else {
			seg <- data.frame(x=c(0,nrow(mat)),y=c(0,0),
        xend=c(nrow(mat),nrow(mat)),yend=c(0,nrow(mat)))
			hggplot3 <- hggplot  +
			ggplot2::geom_text(aes(Var2, Var1, label = value), 
			color = c("black"), size = size.txt) +
			ggplot2::geom_segment(data=seg, aes(x=x+0.5,xend=xend+0.5, 
			y=y+0.5,yend=yend+0.5), size=1, inherit.aes=F)	+ 
			ggplot2::scale_color_manual(values=colvect)
			return(hggplot3)
	 }  
	}  else { return(hggplot) }
}

