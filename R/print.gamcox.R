print.gamcox <- function(x,...){
  fit = x
  cat("Call:\n")
  print(fit$call)
  cat("\nModel:\n")
  if(!is.null(fit$formula)){
    print(fit$formula,showEnv = F)
  }else{
    data = as.matrix(fit$data)
    response.name=colnames(data)[1:2]
    other.name=colnames(fit$X)
    coords = paste("lo(",other.name[1],",",other.name[2],")",sep="")
    covariates = paste(other.name[-c(1,2)],collapse="+")
    print(as.formula(paste("Surv(",response.name[1],",",response.name[2],")~",covariates,"+",coords)),showEnv = F)
  }
  cat("\nCoefficients:\n")
  print(fit$coefficients)
  
  cat("\nSummary of spatial effect(centered at median):\n")
  sp = as.vector(fit$smooth+fit$smooth.frame%*%fit$coefficients[grep("lo\\([[:print:]]+\\)",names(fit$coefficients))])
  print(summary(sp-median(sp)))
  cat("\nDegrees of Residual Freedom:",fit$df.residual,"\n")
  cat("Residual Deviance:", fit$deviance,"\n")
}