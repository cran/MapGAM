colormap <-
function(modgamobj, map=NULL, add=F, mapmin=NULL, mapmax=NULL, arrow=T, ptsize=0.9, alpha=0.05) {
	if (!is.null(modgamobj$OR)) fitvals=modgamobj$OR else fitvals=modgamobj$fit 
	results = cbind(modgamobj$grid,fitvals)
	if (!is.null(modgamobj$pointwise)) results = cbind(results,modgamobj$pointwise)
	if (is.null(mapmin)) mapmin=min(results[,3])
	if (is.null(mapmax)) mapmax=max(results[,3])
	dataXmin=min(results[,1])
	dataXmax=max(results[,1])
	dataYmin=min(results[,2])
	dataYmax=max(results[,2])
	offsetX=(dataXmax-dataXmin)/5
	offsetY=(dataYmax-dataYmin)/5
	qu = seq(mapmin,mapmax,length=2251)
	cp = rainbow(2252,start=0,end=0.66)
	col.seq = rev(cp)
	grad = cut(results[,3],breaks=c(0,qu,Inf),labels=F)
	if (add==T) points(modgamobj$grid,col=col.seq[grad],pch=15,cex=ptsize) else {
		par(mai=c(0,0,0.75,0))
		plot(modgamobj$grid,col=col.seq[grad],pch=15,cex=ptsize,type="p",xlim=c(dataXmin,dataXmax),
			 ylim=c((dataYmin-0.8*offsetY),dataYmax),ann=F,axes=F)
	}
	if (!is.null(map)) {
		if (class(map)=="SpatialPolygonsDataFrame") plot(map, border=gray(0.3), add=T) else {
			if (class(map)!="map") warning("map class not recognized--attempting generic plot")
			lines(map,col=gray(0.3),type="l")
		}
	}
	# CONTOUR LINES
	if (length(results)>=4) {	# if column exists for permutation ranks
		df=cbind(results[(order(results[,2])),1:2],results[(order(results[,2])),4])
		names(df)[3]="perm"
		df = df[!is.na(df$perm),]
		cx=as.matrix(unique(df[(order(df[,1])),1]))
		cy=as.matrix(unique(df[,2]))
		cz=as.matrix(reshape(df,v.names=names(df)[3],idvar=names(df)[1],timevar=names(df)[2],direction="wide")[,-1])
		cz=cz[order(unique(df[,1])),]
		levs = c(alpha/2,1-alpha/2)
		contour(x=cx,y=cy,z=cz,levels=levs,color="black", drawlabels=FALSE, lwd=3, lty="twodash", add=TRUE)
	}

	# LEGEND
	if (!is.null(modgamobj$OR)) leglab = "Odds Ratios" else leglab = "Predicted Values"
	if (!is.null(modgamobj$m)) leglab = paste(modgamobj$m,leglab)
	len=(dataXmax-dataXmin)*2/3
	stop=trunc(len*(1-mapmin)/(mapmax-mapmin)/(len/2252))
	points(dataXmin+(1:2252*len/2252),rep(dataYmin-0.5*offsetY,2252),cex=1.6,col=cp[2253-1:2252],pch=15)
	points(dataXmin+len*(1-mapmin)/(mapmax-mapmin),(dataYmin-0.5*offsetY),cex=.8,col="black", pch="|", lwd=2)
	text(x=dataXmin,y=(dataYmin-0.5*offsetY),pos=1,labels=format(round(mapmin,2),digits=3),cex=.8)
	text(x=dataXmin+len,y=(dataYmin-0.5*offsetY),pos=1,labels=format(round(mapmax,2),digits=3),cex=.8)
	if (!is.null(modgamobj$OR)) text(x=dataXmin+len*(1-mapmin)/(mapmax-mapmin),
		y=(dataYmin-0.5*offsetY),pos=1,labels="1.00",cex=.8)	# only add line at 1 for ORs
	text(x=dataXmin+len/2,y=(dataYmin-0.5*offsetY),pos=3,labels=leglab,cex=1) 

	# NORTH ARROW
	if (arrow) {
		points(dataXmax,dataYmin-0.65*offsetY,cex=1.2,col="black", pch="|", lwd=2)
		points(dataXmax-0.01*offsetX,dataYmin-0.5*offsetY,cex=1.2,col="black", bg="black", pch=24)
		text(x=dataXmax,y=(dataYmin-0.5*offsetY),"N",cex=1,pos=3)
	}

	# SCALE BAR
	if (is.null(map)==F & class(map)=="SpatialPolygonsDataFrame") {
		d=0.17*(dataXmax - dataXmin)
		if (d>1000) {d=d/1000; units="km"} else {units="m"}
		d=ceiling(d)
		d=5*(d%/%5)
		f=(d*1000)/(dataXmax - dataXmin)	
		d=paste(as.character(format(d,digits=2)),units)
		text(x=dataXmin+1.125*len+0.5*f*(dataXmax - dataXmin),y=(dataYmin-0.5*offsetY),d,cex=1,pos=1)
		points(dataXmin+1.125*len,dataYmin-0.5*offsetY,cex=1.2,col="black", pch="|", lwd=2)
		points(dataXmin+1.125*len+f*(dataXmax - dataXmin),dataYmin-0.5*offsetY,cex=1.2,col="black",
			   pch="|", lwd=2)
		lines(c(dataXmin+1.125*len,dataXmin+1.125*len+f*(dataXmax - dataXmin)),
			  c(dataYmin-0.5*offsetY,dataYmin-0.5*offsetY),col="black",lwd=3)
	}
}
