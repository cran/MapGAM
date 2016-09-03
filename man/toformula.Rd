\name{toformula}
\alias{toformula}
\title{
Build a formula based on data for \code{\link{modgam}} function}
\description{
Based on the arguments in \code{\link{modgam}}, build a formula used for model fitting.
}
\usage{
toformula(formula, data,m="adjusted",surv = FALSE, span=0.5,degree=2,
          smooth = TRUE,offset)
}
\arguments{
  \item{formula}{
User specified formula expression (for \code{\link{modgam}})
}
  \item{data}{
A well structured data frame with the two geolocation parameters in the first and second columns. \code{data} is ingnored if \code{formula} is specified
}
  \item{m}{
Use \code{m="adjusted"} to specify a model including adjusted confounders. Use \code{m="unadjusted"} to specify a unadjusted model
}
  \item{surv}{
\code{surv=TRUE} specifies a model for censored surival dataset. 
}
  \item{span, degree}{
Smoothing parameters used for smoothing functions (see details in \code{\link{modgam}}) if \code{smooth=TRUE}.
}
  \item{smooth}{
\code{smooth=TRUE} specifies a smooth term in the model.
}
  \item{offset}{
Name for offset.
}
}

\value{
The function returns the formula for the model.
}


\author{
Lu Bai

Send bug reports to \email{sbartell@uci.edu}.
}
 

\seealso{
\code{\link{modgam}}, .
}
\examples{
\donttest{
data(CAdata)
toformula(data=CAdata, surv=TRUE)
}
}

\keyword{misc}
\keyword{smooth}
