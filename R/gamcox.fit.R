gamcox.fit <- function(Y,X,smooth.frame,weights,span=0.5, I.span = 0.2,degree=1,loess.trace="exact",Maxiter=40, tol=1e-7){
  N <- length(Y$time)
  fit <- coxph(Surv(Y$time,Y$event)~X,weights=if(missing(weights))rep(1,N)else weights)
  eta <- fit$linear.predictors
  adjust = TRUE
#  if(dim(X)[2]==2)adjust=FALSE
  which = grep("lo\\([[:print:]]+\\)",colnames(X))
  ls <- dls(Y,X,which,eta,I.span,adjust)
  w <- ls$w*(ls$w>=0)
  deltaeta <- ls$deltaeta
  dev.old <- -ls$l*2
  Max.backfitting <- Maxiter
  if(!adjust) Max.backfitting = 1
  for(i in 1:Maxiter){
    Z <- eta+deltaeta
    residuals <- Z
    fit <- list(fitted.values=rep(0,N),coefficients=0)
    s <- rep(0,N)
    for (j in 1:Max.backfitting){
      R <- residuals + fit$fitted.values
#      if(adjust){
        fit <- lm(R~X,weights=w)
        residuals <- fit$residuals
        old <- s
        R <- residuals + s
#      }
      smooth.fit <- loess(R ~ smooth.frame,weights=w,span=span,degree=degree,normalize=FALSE,control=loess.control(trace.hat=loess.trace))
      s <- smooth.fit$fitted-mean(smooth.fit$fitted)
      residuals <- R-s
#      if(adjust){
        deltaf <- weighted.mean((s - old)^2,w)
        RATIO <- sqrt(deltaf/sum(s^2))
        if (RATIO <= tol) break;
#      }
    }  
    eta <- Z-residuals
    ls <- dls(Y,X,which,eta,I.span,adjust)	
    dev.new <- -2*ls$l
    if(i==1 | dev.new<=dev.old){
      rslt.smooth.fit=smooth.fit
      rslt <- list()
#      if(adjust){
        rslt$cov <- vcov(fit)[-1,-1]
        rslt$coefficients <- fit$coefficients[-1]
#      }
      rslt$weights <- w
      rslt$additive.predictors <- eta
      rslt$smooth <- s
      dev.out <- dev.new
      rslt$residuals <- residuals
      deltaeta <- ls$deltaeta
      w <- ls$w*(ls$w>=0)
    }
    if(dev.new>dev.old | (abs(dev.old - dev.new)/(dev.old + 0.1) <=tol)) break;
    dev.old <- dev.new
  }
  if(i == Maxiter)warnings(paste("The process does not converge in",Maxiter,"iterations"))
  rslt$deviance <- dev.out
  dfs = predict(rslt.smooth.fit,se=TRUE)
  rslt$var = dfs$se^2
  rslt$df.residual = dfs$df-if(!is.null(rslt$coefficients))length(rslt$coefficients) else 0
  rslt$aic = rslt$deviance + 2*(N-rslt$df.residual)
  if(!is.null(rslt$coefficients))
    names(rslt$coefficients) <- colnames(X)
  rslt
}