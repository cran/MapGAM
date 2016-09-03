gamcox<- function(formula,data,subset,weights,span=0.5,I.span=0.2,degree=1,loess.trace="exact", Maxiter=40, tol=1e-7){
  call <- match.call()
  #if (!missing(formula)){
    if (missing(data)) 
      data <- environment(formula)
    mf <- match.call()
    m <- match(c("formula", "data","subset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mt <- if (missing(data)) 
      terms(formula, c("lo"))
    else terms(formula, c("lo"), data = data)
    mf$formula <- mt
    mf <- eval(mf, parent.frame())
    order = attr(mt,"specials")$lo
    if(is.null(order)){
      if(missing(weights))
      fit = coxph(formula,data[subset,])
      else fit = coxph(formula,data,weights,subset)
      rslt <- fit
      rslt$additive.predictors = fit$linear.predictors
      if(!is.null(fit$coefficients)){
        rslt$deviance = -2*fit$loglik[2]
        rslt$var=fit$var
        rslt$aic = AIC(fit)
        rslt$df.residual = fit$n - if(!is.null(fit$coefficients))length(fit$coefficients) else 0
      }else{
        rslt$deviance = -2*fit$loglik
        rslt$df.residual = fit$n
      }
      class(rslt)<-"coxph"
      return(rslt)
    }else{
      ncols <- attr(mf[[order]],"ncols")
      spanf = if(!missing(span)) span else NULL
      degreef = if(!missing(degree)) degree else NULL
      los <- attr(mt,"term.labels")[order-1]
      los <- gsub("[[:blank:]]","",los)
      start <- gregexpr("\\(",los)[[1]]; end <-rev(gregexpr("\\)",los)[[1]]) 
      losub <- substr(los,start[1]+1,end[1]-1)
      parts <- strsplit(losub,"\\([^()]*\\)(*SKIP)(*F)|\\h*,\\h*", perl=T)
      coords.name <- parts[[1]][1:2]
      if(is.null(dim(mf[[order]]))|dim(mf[[order]])[2]<2) stop("Two predictors must be specified in lo()")
      if(ncols>2)
        warning(paste("Only the first two variables", paste(coords.name,collapse=" and "),"will be used for smoothing"))
      mf[[order]] <- mf[[order]][,1:ncols]
      colnames(mf[[order]]) <- coords.name
      if(length(parts[[1]])>2){
        for(i in 3:length(parts[[1]]))
          eval(parse(text=parts[[1]][i]))
      }
      if(!is.null(spanf))if(spanf!=span) 
        warning(paste("span size of",span,"in the formula will be used instead of value of argument sp=", sp))
      if(!is.null(span)) sp = span
      if(!is.null(degreef))if(degreef!=degree)
        warning(paste("degree=", degree,"in the formula will be used instead of value of argument degree=", degreef))
      
      Y = list(time=mf[[1]][,1],event = mf[[1]][,2])
      X.matrix <- if(!is.empty.model(mt))model.matrix(mt,mf,contrasts)[,-1]
#      X = as.matrix(cbind(mf[[order]],X.matrix[,-grep("lo\\([[:print:]]+\\)",colnames(X.matrix))]))
      smooth.frame = mf[[order]]
#      colnames(X)[-c(1,2)] <- colnames(as.matrix(mf))[-c(1,2,order+1,order+2)]
    }
#  }else{
#    if(min(dim(data))<4)
#      stop("data must include at least 4 columns")
#    Y = list(time=data[subset,1],event=data[subset,2])
    #X = as.matrix(data[subset,-c(1,2)])
#    X = model.matrix(as.formula(paste("~",paste(names(data)[-c(1,2)],collapse="+"),"-1")),data[subset,],contrasts)
#  }
  if(!missing(weights)){
    if(any(weights<0)) stop("Negative weights not allowed")
    weights = weights[subset]
  }
  fit = gamcox.fit(Y,X.matrix,smooth.frame,weights,span,I.span,degree,loess.trace, Maxiter, tol)
  rslt <- fit
#  if(!missing(formula)){
    rslt$formula=formula
    rslt$terms=mt
#  }
  rslt$xlevels = .getXlevels(mt,mf)
  rslt$data <- data
  rslt$span <- span
  rslt$I.span <- I.span
  rslt$Maxiter <- Maxiter
  rslt$tol<- tol
  rslt$degree <- degree
  rslt$loess.trace = loess.trace
  rslt$X=X.matrix
  rslt$Y = Y
  rslt$smooth.frame = mf[[order]]
  rslt$weights = if(missing(weights))NULL else weights
  rslt$call <- call
  class(rslt)<-"gamcox"
  rslt
}