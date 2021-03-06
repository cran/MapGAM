#***********************************************************************************
#
# Print modgam object
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
print.modgam <- function(x,...){
  fit = x
  cat("Call:\n")
  ## print the model
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
  print(fit$gamobj$coefficients)
}