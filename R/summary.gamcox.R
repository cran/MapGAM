#***********************************************************************************
#
# Summarize the gamcox Object
# Copyright (C) 2016, The University of California, Irvine
#
# This library is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this library??? if not, <http://www.gnu.org/licenses/>.
#
#*******************************************************************************
summary.gamcox <- function(object,...){
  fit = object
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
  if(!is.null(fit$global.pvalue)) {
    if(fit$global.pvalue>=1e-5)
      cat("p value for testing the global spatial effect:", fit$global.pvalue,"\n")
    else 
      cat("p value for testing the global spatial effect: <1e-5 \n")
  }
}