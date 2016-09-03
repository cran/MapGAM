plot.modgam <- function(x, map = NULL, exp = F, add = F, intervals=T,
                        mapmin = NULL, mapmax = NULL, border.gray = 0.3, contours=c("none","response","permrank","interval"), 
                        contours.drawlabels=FALSE,contours.lty=1, contours.lwd=1, contours.levels, contours.labcex=0.7, 
                        arrow=T, axes=F, ptsize=0.9, alpha=0.05, mai,legend.name = "predicted values",legend.cex=1,...){
  modgamobj = x
  leglab = legend.name
  if(missing(legend.name)){
    if(!is.null(modgamobj$family)){
      if(modgamobj$family[1]=="survival") leglab = if(exp) "hazard ratio" else "log hazard ratio"
      if(modgamobj$family[1]=="binomial"&modgamobj$family[2]=="logit") leglab = if(exp) "odds ratio" else "log odds ratio"
      if(modgamobj$family[1]=="poisson"&modgamobj$family[2]=="log") leglab = if(exp) "risk ratio" else "log risk ratio"
    }
  }
  contours = contours[1]
  if(contours!="none"){
    if (contours=="permrank" && is.null(modgamobj$pointwise.permt)) {
      warning("permrank contours omitted because permute=0 in modgam")
      contours = "none"
    }
    if (contours=="interval" && any(c(is.null(modgamobj$conf.low),is.null(modgamobj$conf.high)))){
      warning("interval contours omitted because no conf.low or conf.high in modgam")
      contours = "none"
    }
    if(!is.element(contours,c("none","response","permrank","interval"))){
      warning("contours omitted because contours type not recognized")
      contours = "none"
    }
    if(contours == "response"){
      contours = "fit"
      if(missing(contours.drawlabels)) contours.drawlabels = TRUE
    }
    if(contours == "permrank"){
      contours = "pointwise.permt"
      if(missing(contours.levels)) contours.levels = c(alpha/2, 1-alpha/2)
      if(missing(contours.lwd)) contours.lwd = 2
      if(missing(contours.lty)) contours.lty = 1
    }
    if(contours == "interval"){
      if(missing(contours.levels)) contours.levels = c(-1,1)
      if(missing(contours.lwd)) contours.lwd = 2
      if(missing(contours.lty)) contours.lty = 1
    }
  }
  legend.add.line = if(exp) 1 else 0
  if(intervals & all(c(!is.null(modgamobj$conf.low),!is.null(modgamobj$conf.high)))){
    if(is.null(mapmin))
      mapmin=min(if(exp)modgamobj$exp.conf.low else modgamobj$conf.low,rm.na=TRUE)
    if(is.null(mapmax))
      mapmax=max(if(exp)modgamobj$exp.conf.high else modgamobj$conf.high,rm.na=TRUE)
    mmai = if(missing(mai)) c(0,0,0.3,0) else mai
    legend.cex = legend.cex*1.4
    op.mfrow = par()$mfrow
    tempobj1 = modgamobj
    tempobj1$fit = modgamobj$conf.low; tempobj1$exp.fit = modgamobj$exp.conf.low
    tempobj2 = modgamobj
    tempobj2$fit = modgamobj$conf.high; tempobj2$exp.fit = modgamobj$exp.conf.high
    par(mfrow = c(1,3))
    colormap(tempobj1, map, exp, add, mapmin, mapmax, border.gray,contours, contours.drawlabels, contours.lty,
             contours.lwd, contours.levels, contours.labcex, 0, arrow, axes, ptsize,mmai,leglab,legend.cex, legend.add.line,
             ...)
    title(main=paste(round((1-modgamobj$predobj$level)*100,2),"% CI (lower)"),cex.main=legend.cex)
    colormap(modgamobj, map, exp, add, mapmin, mapmax, border.gray,contours, contours.drawlabels, contours.lty,
             contours.lwd, contours.levels, contours.labcex, 0, arrow, axes, ptsize,mmai,leglab,legend.cex, legend.add.line,
             ...)
    title(main="Point Estimate",cex.main=legend.cex)
    colormap(tempobj2, map, exp, add, mapmin, mapmax, border.gray,contours, contours.drawlabels, contours.lty,
             contours.lwd, contours.levels, contours.labcex, 0, arrow, axes, ptsize,mmai,leglab,legend.cex, legend.add.line,
             ...)
    title(main=paste(round((1-modgamobj$predobj$level)*100,2),"% CI (higher)"),cex.main=legend.cex)
    par(mfrow=op.mfrow)
  }else
    colormap(modgamobj, map, exp, add, mapmin, mapmax, border.gray,contours, contours.drawlabels, contours.lty,
             contours.lwd, contours.levels, contours.labcex,0, arrow, axes, ptsize, mai,leglab,legend.cex, legend.add.line,
             ...)
}