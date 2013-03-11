\name{sampcont}
\alias{sampcont}
\title{
Unmatched Control Sampling
}
\description{
Take all cases and a random sample of controls from a data frame.  Simple random sampling and stratified random sampling are available.  For statified random sampling, strata can be defined by region, or by region and time.  If no specific regions are specified then the function will create a regular grid for sampling.     
}
\usage{
sampcont(rdata, type = "stratified", regions = NULL, times = NULL, n = 1, nrow = 100, ncol = 100)
}
\arguments{
  \item{rdata}{
A data frame with the outcome (coded as 0/1) in the 1st column, and the geocoordinates (e.g., X and Y) in the 2nd and 3rd columns.  Additional columns are not used in the sampling scheme but are retained in the sampled data frame.  
}
  \item{type}{
"stratified" (default) or "simple".  If "simple" then a simple random sample of \code{n} controls (rows of \code{rdata} with outcome=0) is obtained.  If "stratified" then a stratified random sample of controls is obtained, with up to \code{n} controls per stratum.  Sampling strata are defined by the \code{regions} and \code{times} arguments.  All cases (rows with outcome=1) are taken for the sample regardless of the value supplied for \code{type}.  
}
  \item{regions}{
A vector of length equal to the number of rows in \code{rdata}, used to construct sampling strata.  If \code{regions = NULL} then the function will define \code{regions} as a vector of specific grid cells for each row in \code{rdata}, on a regular grid with \code{nrow} rows and \code{ncol} columns.  If \code{times = NULL} then the nonempty regions are used as the sampling strata.  If \code{times} is a vector, then the sampling strata are all nonempty combinations of \code{regions} and \code{times}.      
}
  \item{times}{
A vector of length equal to the number of rows in \code{rdata}, used to construct sampling strata.  If \code{times = NULL} then the sampling strata are defined only by the \code{regions} argument.  If \code{times} is a vector, then the sampling strata are all nonempty combinations of \code{regions} and \code{times}.  Continuous times should generally be binned before being passed through this argument, as there are no efficiency gains if each value in \code{times} is unique.     
}
  \item{n}{
The number of controls to sample from the eligible controls in each stratum.  All available controls will be taken for strata with fewer than \code{n} eligible controls.  
}
  \item{nrow}{
The number of rows used to create a regular grid for sampling regions.  Only used when \code{regions = NULL}.
}
  \item{ncol}{
The number of columns used to create a regular grid for sampling regions.  Only used when \code{regions = NULL}.
}
}
\value{
  \item{rdata}{
A data frame with all cases and a random sample of controls.  
}
  \item{w}{
Inverse probability weights for the rows in \code{rdata}.  Important to include as weights in subsequent analyses.    	
}
  \item{ncont}{
The total number of controls in the sample. 
}
}
\author{
Scott Bartell
\email{sbartell@uci.edu}.
}

\seealso{
\code{\link{modgam}}
}

\examples{
#### load beertweets data, which has 719 cases and 9281 controls
data(beertweets)
# take a simple random sample of 1000 controls
samp1 <- sampcont(beertweets, type="simple", n=1000)

# take a stratified random sample of controls on a 80x50 grid
# requires PBSmapping package
samp2 <- NULL
\donttest{
if(require(PBSmapping)) samp2 <- sampcont(beertweets, nrow=80, ncol=50)

# Compare locations for the two sampling designs (cases in red)
par(mfrow=c(2,1), mar=c(0,3,4,3))
plot(samp1$rdata$longitude, samp1$rdata$latitude, col=3-samp1$rdata$beer,
	cex=0.5, type="p", axes=FALSE, ann=FALSE)
# Show US base map if maps package is available
mapUS <- require(maps)
if (mapUS) map("state", add=TRUE)
title("Simple Random Sample, 1000 Controls")

if (!is.null(samp2)) {
	plot(samp2$rdata$longitude, samp2$rdata$latitude, 
		col=3-samp2$rdata$beer, cex=0.5, type="p", axes=FALSE, 
		ann=FALSE)
	if (mapUS) map("state", add=TRUE)
	title(paste("Spatially Stratified Sample,",samp2$ncont,"Controls"))
	}
}
par(mfrow=c(1,1))

## Note that weights are needed in statistical analyses
# Prevalence of cases in sample--not in source data
mean(samp1$rdata$beer)		 
# Estimated prevalence of cases in source data 
weighted.mean(samp1$rdata$beer, w=samp1$w)	
## Do beer tweet odds differ below the 36.5 degree parallel?
# Using full data
glm(beer~I(latitude<36.5), family=binomial, data=beertweets) 
# Stratified sample requires sampling weights 
if (!is.null(samp2)) glm(beer~I(latitude<36.5), family=binomial, 
	data=samp2$rdata, weights=samp2$w)
}

\keyword{misc}
