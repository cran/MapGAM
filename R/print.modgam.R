print.modgam <- function(x,...){
  fit = x
  cat("Call:\n")
  print(fit$call)
  cat("\nModel:\n")
  if(!is.null(fit$gamobj$formula)){
    print(fit$gamobj$formula,showEnv = F)
  }else{
    data = as.matrix(fit$gamobj$data)
    response.name=colnames(data)[1:2]
    other.name=colnames(fit$gamobj$X)
    coords = paste("lo(",other.name[1],",",other.name[2],")",sep="")
    covariates = paste(other.name[-c(1,2)],collapse="+")
    print(as.formula(paste("Surv(",response.name[1],",",response.name[2],")~",covariates,"+",coords)),showEnv = F)
  }
  if(tolower(fit$family)[1]!="survival")
    cat(paste(c("Family:", "Link:"), fit$family[1:2]), fill = T)
  else cat(paste("span: ",fit$gamobj$span,"\n"))
  cat("\nCoefficients:\n")
  #print(fit$gamobj$coefficients)
  estimate = fit$gamobj$coefficients
  if(any(class(fit$gamobj)=="gam")) sd = sqrt(diag(vcov(fit$gamobj)))
  else sd = sqrt(diag(fit$gamobj$cov))
  zstat = estimate/sd
  pvalue = 2*pnorm(-abs(zstat))
  
  mat = cbind(estimate,sd,zstat,pvalue)
  colnames(mat)=c("Estimate", "Std.Error", "z value", "P(>|z|)")
  print(mat)
  cat("\nDegrees of Residual Freedom:",fit$gamobj$df.residual,"\n")
  cat("Residual Deviance:", fit$gamobj$deviance,"\n")
  cat("AIC:", fit$gamobj$aic,"\n")
  if(fit$global.pvalue>=1e-5)
    cat("p value for testing the global spatial effect:", fit$global.pvalue,"\n")
  else 
    cat("p value for testing the global spatial effect: <1e-5 \n")
  cat("\n Spatial effect predictions:\n")
  print(summary(fit$fit))
  
  if(!is.null(fit$conf.low)){
    cat("\nSpatial effect",100*(1-fit$predobj$level), "% lower interval:\n")
    print(summary(fit$conf.low))
    cat("\nSpatial effect",100*(1-fit$predobj$level), "% higher interval:\n")
    print(summary(fit$conf.high))
  }
  
  if(!is.null(fit$N.permt))
    if(fit$global.permt>=1/fit$N.permt)
      cat(paste("\np value for testing the global spatial effect based on",fit$N.permt,"permutations:",
              fit$global.permt,"\n"))
    else 
      cat(paste("\np value for testing the global spatial effect based on",fit$N.permt,"permutations: <",
                1/fit$N.permt,"\n"))
}