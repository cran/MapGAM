colormap <-
function(obj, map=NULL, exp=F, add=F,  mapmin=NULL, mapmax=NULL, border.gray = 0.3, contours="none", contours.drawlabels=FALSE,
	contours.lty=1, contours.lwd=1, contours.levels, contours.labcex=0.7,interval.exclude=0, arrow=T, axes=F, ptsize=0.9, 
	mai, legend.name = "predicted values",legend.cex=1,legend.add.line,alpha=0.05,...) {
  op.mai <- par()$mai
  op.xpd <- par()$xpd
	if (exp) fitvals=exp(obj$fit) else fitvals=obj$fit 
	results = cbind(obj$grid,fitvals)
	if (is.null(mapmin)) mapmin=min(results[,3], na.rm=TRUE)
	if (is.null(mapmax)) mapmax=max(results[,3], na.rm=TRUE)
	if (!is.null(map)) {
		if (add==T) {add=F; warning("add=T ignored because the map argument was used")}  
		if (class(map)=="map") maprange = map$range else
		if (class(map)=="SpatialPolygonsDataFrame") {
			mapcenter = apply(bbox(map),1,mean)
			center2edges = c(-1,1)*max(apply(bbox(map),1,diff))/2
			maprange = c(mapcenter[1]+center2edges,mapcenter[2]+center2edges)
		} else
		maprange = c(Inf,-Inf,Inf,-Inf)
	}
	dataXmin=min(results[,1],if(!is.null(map)) maprange[1])
	dataXmax=max(results[,1],if(!is.null(map)) maprange[2])
	dataYmin=min(results[,2],if(!is.null(map)) maprange[3])
	dataYmax=max(results[,2],if(!is.null(map)) maprange[4])
	offsetX=(dataXmax-dataXmin)/5
	offsetY=(dataYmax-dataYmin)/5
	qu = seq(mapmin,mapmax,length=2251)
	cp = rainbow(2252,start=0,end=0.66)
	col.seq = rev(cp)
	grad = cut(results[,3],breaks=c(0,qu,Inf),labels=F)
	par(xpd=TRUE)
	if (add==T) points(obj$grid,col=col.seq[grad],pch=15,cex=ptsize,...) else
	if (axes==F) {
	  mmai = if(missing(mai))c(0,0.35,0,0.35) else mai
		par(mai=mmai)
		plot(obj$grid,col=col.seq[grad],pch=15,cex=ptsize,type="p",
			 xlim=c(dataXmin-0.4*offsetX,dataXmax+0.4*offsetX),
			 ylim=c(dataYmin-0.8*offsetY,dataYmax),ann=F,axes=F,...)
	   }  else
	if (axes==T) {
#		par(mai=c(1.52,0.82,0.82,1.42))
	  mmai = if(missing(mai))c(1,0.82,1,0.5) else mai
	  par(mai=mmai)
		plot(obj$grid,col=col.seq[grad],pch=15,cex=ptsize,type="p",
			xlab="", xaxt="n")
		axis(3)
		mtext(names(obj$grid[1]), side=3, line=2)
		#mtext(names(obj$grid[2]), side=2, line=2)
	}
	if (!is.null(map)) {
		if (class(map)=="SpatialPolygonsDataFrame") plot(map, border=gray(border.gray), add=T) else {
			if (class(map)!="map") warning("map class not recognized--attempting generic plot")
			lines(map,col=gray(border.gray),type="l")
		}
	}
	
	# CONTOUR LINES / ISOBOLES
	if (contours!="none") {		# if contour lines should be drawn
	  result_contour = NULL
	  if(contours=="interval"){
	    if(any(c(is.null(obj$conf.low),is.null(obj$conf.high))))
	      warning("interval contours omitted because no conf.low or conf.high in obj")
	    else{
	      if(!exp)
	        result_contour = apply(cbind(obj$conf.low,obj$conf.high),1,
	         function(x){if(any(is.na(x))) 0 else 
	           sign(x[1]-interval.exclude)*((x[1]-interval.exclude)*(x[2]-interval.exclude)>=0)})
	      else
	        result_contour = apply(cbind(exp(obj$conf.low),exp(obj$conf.high)),1,
	         function(x){if(any(is.na(x))) 0 else sign(x[1]-exp(interval.exclude))*
	             ((x[1]-exp(interval.exclude))*(x[2]-exp(interval.exclude))>=0)})
	      levels= c(0,1)
	      labels = c("lower","higher")
	      drawlabels=if(!missing(contours.drawlabels))contours.drawlabels else FALSE
	      lwd=if(!missing(contours.lwd)) contours.lwd else 2
	      lty=if(!missing(contours.lty)) contours.lty else 1
	    } 
	  }else{
	    if(!is.element(contours,names(obj))){
	      if(contours =="permrank"){ # for back-compatibility
	        contours = "pointwise.permt"
	        contours.levels = c(alpha/2, 1-alpha/2)
	      }else{
	        if(contours == "response")  ## for back-compatibility
	        contours = "fit"
	        else warning(paste("contours omitted because no",contours,"in obj"))
	      }
	    }
	    if(is.element(contours,names(obj))){
	      result_contour = obj[[contours]]
	      if(exp) result_contour = exp(result_contour)
	      levels= if(!missing(contours.levels)) contours.levels else pretty(range(result_contour,finite=TRUE),10) 
        labels = levels
	      drawlabels=if(!missing(contours.drawlabels))contours.drawlabels else FALSE
	      lwd=contours.lwd
	      lty=contours.lty
	    }
	  }
	  if(!is.null(result_contour)){
    	df=cbind(results[(order(results[,2])),1:2],result_contour[(order(results[,2]))])
  		cx=as.matrix(unique(df[(order(df[,1])),1]))
  		cy=as.matrix(unique(df[,2]))
  		cz=as.matrix(reshape(df,v.names=names(df)[3],idvar=names(df)[1],
  			timevar=names(df)[2],direction="wide")[,-1])
  		cz=cz[order(unique(df[,1])),]	
  		contour(x=cx, y=cy, z=cz, add=TRUE,levels=levels, labels = labels,drawlabels=drawlabels, lwd=lwd, lty=lty, labcex = contours.labcex)
	  }
	}
		
	# LEGEND
	#if (!exp) leglab = legend.name else leglab = "Predicted Values"
	leglab = legend.name
	# if (!is.null(obj$m)) leglab = paste(obj$m,leglab)
  fY = 0.5 + 0.5*(axes==T)
	ypos = dataYmin-fY*offsetY
	len = (dataXmax-dataXmin)*7/12
	points(dataXmin+(1:2252*len/2252),rep(ypos,2252),cex=1.6,col=cp[2253-1:2252],pch=15)
	text(x=dataXmin,y=ypos,pos=1,labels=format(round(mapmin,2),digits=3),cex=.8*legend.cex)
	text(x=dataXmin+len,y=ypos,pos=1,labels=format(round(mapmax,2),digits=3),cex=.8*legend.cex)
	text(x=dataXmin+len/2,y=ypos,pos=3,labels=leglab,cex=legend.cex) 
  if(!missing(legend.add.line)){  # Add one line at legend
    points(dataXmin+len*(legend.add.line-mapmin)/(mapmax-mapmin),
           ypos,cex=.8*legend.cex,col="black", pch="|", lwd=2)  
    text(x=dataXmin+len*(legend.add.line-mapmin)/(mapmax-mapmin),
         y=ypos,pos=1,labels=paste(legend.add.line),cex=.8*legend.cex)
  }
	# NORTH ARROW
	if (arrow) {
		points(dataXmax,dataYmin-0.65*offsetY,cex=1.2,col="black", pch="|", lwd=2)
		points(dataXmax-0.01*offsetX,ypos,col="black", bg="black", pch=24)
		text(x=dataXmax,y=ypos,"N",cex=legend.cex,pos=3)
	}

	# SCALE BAR
	if (is.null(map)==F & class(map)=="SpatialPolygonsDataFrame") {
		d = 0.17*(dataXmax - dataXmin)
		d = signif(d,1)
		leftedge = dataXmin+1.125*len
		lines(c(leftedge,leftedge+d),c(ypos,ypos),col="black",lwd=3)
		points(leftedge,ypos,cex=1.2,col="black", pch="|", lwd=2)
		points(leftedge+d,ypos,cex=1.2,col="black",pch="|",lwd=2)
		if (d>=1000) {scdist=d/1000; units="km"} else {scdist=d; units="m"}
		sclab = paste(as.character(format(scdist,digits=2)),units)
		text(x=leftedge+0.5*d,y=ypos,sclab,cex=0.8*legend.cex,pos=1)
	}
	par(mai=op.mai,xpd=op.xpd)
}
