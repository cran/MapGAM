optspan=function (formula,data, offset,spans=seq(0.05,0.95,by=0.05),m = "adjusted", family = binomial(), verbose = TRUE, degree=1,...)
{
  if(!missing(formula)){
    if (missing(data)) 
      data <- environment(formula)
    mf <- match.call()
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mt <- if (missing(data)) 
      terms(formula, c("Surv","lo"))
    else terms(formula, c("Surv","lo"), data = data)
    surv = !is.null(attr(mt,"specials")$Surv)
    order = attr(mt,"specials")$lo
    if(surv){
      family = "survival"
    }
    if(is.character(family))
      if(!surv & tolower(family)[1]=="survival")
        stop("formula misspecified for survival data")
    mf$formula <- mt
    mf <- eval(mf, parent.frame())
    ncols <- attr(mf[[order]],"ncols")
    if(is.null(order) | ncols<2)
      stop("Two parameters must specified in lo() in the formula")
    if(length(order)>1)
      stop("Only one smoothing function could be included in the model")
    los <- attr(mt,"term.labels")[order-1-surv]
    los <- gsub("[[:blank:]]","",los)
    start <- gregexpr("\\(",los)[[1]]; end <-rev(gregexpr("\\)",los)[[1]]) 
    losub <- substr(los,start[1]+1,end[1]-1)
    parts <- strsplit(losub,"\\([^()]*\\)(*SKIP)(*F)|\\h*,\\h*", perl=T)
    coords.name <- parts[[1]][1:2]
    if(ncols>2)
      warning(paste("Only the first two variables", paste(coords.name,collapse=" and "),"will be used for smoothing"))
    
    if(length(mf)==2) m = "unadjusted" else m="adjusted"
  }
  dim.min = 3; surv = FALSE
  if(is.character(family)){
    if(tolower(family)== "survival") {
      dim.min = 4; surv = TRUE
    }
  }
  if(dim.min==3){
    if(is.character(family)) 
      family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) family = family()
    if (is.null(family$family)) {
      print(family)
      stop("`family' not recognized")
    }
  }
  data.name = colnames(data)
  data.dim = length(data.name)
  if(missing(formula)){
    if(data.dim < dim.min) stop(paste("data must include at least", dim.min,"columns"))
    if(tolower(m)=="adjusted"){
      if(data.dim == dim.min){
        warning("cannot run adjusted model without additional covariates in data")
        m = "unadjusted"
      }
    }else if(tolower(m) == "unadjusted" | tolower(m) == "crude"){
      data = data[,1:dim.min]
      data.name=data.name[1:dim.min]
      data.dim = dim.min
    }else stop(paste("model type", m, "not recognized"))
    coords.name <- data.name[2:3+surv]
  }

  N.span = length(spans)
  span.aic = rep(0,N.span)
  for (S in 1:N.span) {
      span = spans[S]
      fmla = toformula(formula,data, m,surv, span=span, degree=degree,TRUE,offset)
      if(surv)
        fit = eval(substitute(gamcox(fmla,data,span=span,degree=degree,...)))
      else
        fit = eval(substitute(gam(fmla,family=family,data=data,...)))
      span.aic[S] = fit$aic
      if (verbose) cat(paste("The AIC for span=", format(span, nsmall = 2), 
          " is ", format(fit$aic, nsmall = 2), ".", sep = ""), 
          fill = TRUE)
  }
  sp = spans[which.min(span.aic)]
  if(verbose) cat(paste("The optimal span is",sp,"\n"))
  return(sp)
}
