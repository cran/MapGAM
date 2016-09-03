toformula <- function(formula,data,m="adjusted",surv = FALSE, span=0.5,degree=2,smooth = TRUE,offset){
  add = FALSE
  if(missing(formula)){
    names = names(data)
    if(surv) response = paste("Surv(",names[1],",",names[2],")")
    else response = paste(names[1])
    coords = paste("lo(",names[2+surv],",",names[3+surv],",span=",span,",degree=",degree,")")
    if(m=="adjusted"){
      adjust = paste(names[-(1:(surv+3))],collapse=" + ")
      if(smooth)
        fmla = paste(response,paste(coords,"+",adjust),sep="~")
      else
        fmla = paste(response,adjust,sep="~")
    }else{
      if(smooth)
        fmla = paste(response,coords,sep="~")
      else
        fmla = paste(response,"1",sep="~")
    }
  }else{
    fm = Reduce(paste, deparse(formula))
    parts = unlist(strsplit(fm,"\\~"))
    mf = terms(formula)
    parts2 = attr(mf,"term.labels")
    lo.index = grep("lo\\([[:print:]]+\\)",parts2)
    offset.index = grep("offset\\[[:print:]]+\\)",parts2)
    if(length(offset.index)==0)add=TRUE
    if(smooth){
      parts.2 = parts2[lo.index]
      start <- gregexpr("\\(",parts.2)[[1]]; end <-rev(gregexpr("\\)",parts.2)[[1]]) 
      psub <- substr(parts.2,start[1]+1,end[1]-1)
      split = paste(unlist(strsplit(psub,"\\([^()]*\\)(*SKIP)(*F)|\\h*,\\h*", perl=T))[1:2],collapse=",")
      parts.2= paste("lo(",split,",span=",span,",degree=",degree,")")
      parts.2=c(parts.2,parts2[-lo.index])
      parts2 <- paste(parts.2,collapse=" + ")
    }else{
      parts2 <- parts2[-lo.index]
      if(length(parts2)>0)
        parts2 <- paste(parts2,collapse=" + ")
      else parts2 = "1"
    }
    fmla = paste(parts[1],parts2,sep=" ~ ")
    
  }
  if(!missing(offset)&add)fmla =paste(fmla,"offset(",offset,")")
  as.formula(fmla)
}