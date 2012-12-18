\name{optspan}
\alias{optspan}
\title{
Determine the Optimal Span Size for modgam  
}
\description{
Determines the optimal span size for \code{modgam}, a spatial logistic generalized additive model (GAM) with a LOESS smooth for location, by minimizing the AIC.  Currently evaulates span sizes from 5\% to 95\% in increments of 5\%.  
}
\usage{
optspan(rdata, m = "adjusted")
}
\arguments{
  \item{rdata}{
Data set (required).  The data must be structured so that the outcome is in the 1st column and X and Y location values are in the 2nd and 3rd columns.  The outcome column should be coded as 0s and 1s.  Any additional columns will be entered as covariates with linear effects in the gam model. }
  \item{m}{
Model type.  Options are "crude" or "adjusted" (default).  If "crude", only the spatial smooth term for location is included in the model.  If "adjusted", all covariates in the data set (columns >= 4) are included in the model.  
}
}

\value{
The optimal span size, determined by minimizing the AIC
}

\author{
Veronica Vieira, Scott Bartell, and Robin Bliss

Send bug reports to \email{sbartell@uci.edu}.
}

\note{
This function does not return model predictions--only the optimal span size.  To obtain model predictions use the modgam function.   
}

\seealso{
\code{\link{predgrid}}, 
\code{\link{modgam}}, 
\code{\link{colormap}}.
}
\examples{
\donttest{
data(MAdata)
optspan(MAdata, m="crude")
}
}

\keyword{misc}
\keyword{smooth}
\keyword{optimize}