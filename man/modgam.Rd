\name{modgam}
\alias{modgam}
\title{
Fit a Generalized Additive Model with a Two-Dimensional Smooth  (GAM)}
\description{
Fits a crude or adjusted regression on a user-supplied grid for spatial analysis using a generalized additive model with a two-dimensional LOESS smooth for geolocation (or any other two-dimensional predictor).  Includes optional permutation tests for global and local tests of the null hypothesis that the two-dimensional predictor (e.g., geolocation) is not associated with the outcome.  Most applications will pass the output of this function to \code{colormap} to map the resulting odds ratios or other effect estimates.}
\usage{
modgam(rdata, rgrid, family=binomial, permute = 0, conditional = TRUE, 
       m = "adjusted", sp = NULL, keep = FALSE, verbose = TRUE, \dots)
}
\arguments{
  \item{rdata}{
Data set (required). The data must be structured so that the outcome is in the 1st column and the X and Y coordinates for two-dimensional predictor (e.g., geolocation) are in the 2nd and 3rd columns, respectively.  If more than three columns are provided and \code{m = "adjusted"}, the additional columns will be entered as linear predictors in the GAM model. 
}
  \item{rgrid}{
A data frame containing the values at which to generate predictions (required).  X and Y coordinates for the two-dimensional predictor must be in the 1st and 2nd columns.  Additional covariates values for predictions may be provided in additional columns of rgrid.  If \code{m = "adjusted"} and rdata includes covariates not present in rgrid, then for each missing covariate the median value will be taken from rdata.  The \code{predgrid} function can be used to supply an appropriate rgrid for rdata, with the grid clipped according to any specified base map (see examples below).   
}
  \item{family}{
A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.  Note that unlike other packages, MapGAM defaults to the binomial family and logit link, i.e., a logistic model.)
}
  \item{permute}{
The number of permutations of the data set for testing the significance of the two-dimensional predictor. \code{permute = 0} (default) produces no permutation tests.  If \code{permute} > 0, the paired coordinates for the two-dimensional predictor are randomly permuted in order to simulate the distribution of results under the null hypothesis that location is not associated with the outcome.  1000 permutations are recommended for reasonable accuracy of p-values.  WARNING: Because each permutation requires refitting the GAM, permutation tests can be quite slow.    
}
  \item{conditional}{
Logical value indicating whether to run a conditional or unconditional permutation test; this argument is used only when \code{permute} > 0.  The default is a conditional permutation test: the span size is held fixed throughout all permutations even if \code{sp = NULL} to optimize the span size for the original data (see Young et al., 2012).  The unconditional permutation test repeats the span size optimization for each permutation, which is more conservative but takes about 20 times longer to compute.  If \code{conditional = FALSE} the sp argument is ignored when fitting both the original and permuted data sets.               
}
  \item{m}{
Model type for the GAM.  Options are "crude" or "adjusted" (default).  If "crude", only the smooth for the two-dimensional predictor is included in the model (i.e., "crude" is a synonym for "unadjusted").  If "adjusted", all covariates in the data set (columns >= 4) are included in the model.     
}
  \item{sp}{
Span size for the LOESS smooth of the two-dimensional predictor.  If \code{sp = NULL} (default) then an optimal span size will be determined using the optspan function.  
}
  \item{keep}{
Logical value indicating whether to store and return the pointwise odds ratios, and if conditional = FALSE the span size, for every permuted data set.  These values aren't necessary for mapping or cluster identification, and storing them slows the permutation tests, so the default is \code{keep = FALSE}.
}  
  \item{verbose}{
Logical value indicating whether to print the GAM model statement, the percentile rank for the global deviance statistic in the permutation test, and the progress of the permutation test (report completion of every 10 permutations).  The default is \code{verbose = TRUE}.      
}
  \item{\dots}{
Further arguments to be passed to the \code{\link{gam}} function (e.g., weights).  	
}
}

\details{
The model used to fit the data is a generalized additive model with a LOESS smooth for a two-dimensional predictor such as geolocation (Hastie and Tibshirani, 1990; Kelsall and Diggle, 1998; Webster et al., 2006).  Although any family and link function supported by the \code{gam} function is supported, the default binomial family with logit link yields the following model:
\deqn{\ln \left(\frac{\pi_{i}}{1-\pi_{i}}\right) = S(x_{i},y_{i}) + \mathbf{Z_{i}} \boldsymbol{\beta}}{
	ln(\pi / (1-\pi)) = S(x,y) + Z \beta}	
where \eqn{\pi_{i}}{\pi} is the probability that the outcome is 1 for participant \emph{i}, \eqn{x_{i}}{x} and \eqn{y_{i}}{y} are predictor coordinates for participant \emph{i} (i.e., projected distance east and north, respectively, from an arbitrarily defined origin location), \eqn{S(.,.)} is a 2-dimensional smoothing function (currently LOESS), \eqn{\mathbf{Z_{i}}}{Z} is a row vector of covariate values for participant \emph{i}, and \eqn{\boldsymbol{\beta}}{\beta} is a vector of unknown regression parameters (including an intercept).  
When a permutation test is requested, for each permutation the paired X and Y coordinates in the data set are randomly reassigned to participants, consistent with the null hypothesis that the geolocation (or another two-dimensional predictor entered in place of X and Y) is not associated with the outcome.  Note that this procedure intentionally preserves associations between other covariates and the outcome variable so the permutation test reflects the significance of geolocation.  See the references for more details.              
}

\value{
\item{grid}{A data frame with X and Y coordinates from rgrid.}
\item{m}{Whether the GAM was unadjusted (only spatial location as a predictor) or adjusted (included other covariates in addition to spatial location).  "Crude" is a synonym for "unadjusted."}
\item{span}{The span size used for the LOESS smooth.  If \code{keep = TRUE} and \code{conditional = FALSE}, a vector with optimized span sizes for the original data set and each of the permuted data sets.}
\item{gamobj}{The GAM model object from the fit to the data.}
\item{family}{The family and link function used for the model.}
\item{fit}{Predicted values (on the linear predictor scale) for each point on the grid, from the smoothed GAM model.  For example, for a binomial family logit link model the fit values are the log odds of the outcome.  Predicted responses can be obtained from the fit values by applying the inverse link function.}
\item{OR}{Odds ratios for each grid point, comparing fits from the smoothed GAM model with a referent model without the two-dimensional predictor but otherwise identical to the GAM.  These are only provided when the binomial family and logit link are selected.}
If \code{permute > 0} then the following values are also provided:
\item{global}{The p-value based on the deviance statistic comparing models with and without geolocation (or other two-dimensional predictor), using the distribution of deviance statistics from the permuted data sets to represent the null distribution.  For a test of H0: geolocation is unassociated with the outcome variable (adjusting for any covariates if \code{m = "adjusted"}), reject H0 if the percentile rank is below alpha.   WARNING: by default \code{modgam} uses a conditional permutation test which produces inflated type I error rates; Young et al. (2012) recommend using alpha=0.025 to limit the type I error rate to approximately 5\%.} 
\item{pointwise}{For each point on the grid, the percentile rank of the local linear predictor for the model compared to the local linear predictor distribution from the permuted data sets.  This result is needed to define clusters of the outcome (e.g., spatial regions with statistically significant outcomes--e.g., unusually high or low risks.)}
If \code{permute} > 0 and \code{keep = T} then the following values are also provided:
\item{permutations}{A matrix containing null permuted values for each point on the grid, with the results for each permutation in a separate column.  For the binomial family and logit link these are provided as odds ratios, otherwise they are reported as linear predictors.}
}

\references{
Hastie TJ, Tibshirani RJ.  Generalized Additive Models.  (Chapman & Hall/CRC Monographs on Statistics & Applied Probability, Boca Raton, Florida, 1990).  

Kelsall J, Diggle P. Spatial variation in risk of disease: a nonparametric binary regression approach. J Roy Stat Soc C-App 1998, 47:559-573. 	
	
Vieira V, Webster T, Weinberg J, Aschengrau A, Ozonoff D. \href{http://www.ehjournal.net/content/4/1/11}{Spatial analysis of lung, colorectal, and breast cancer on Cape Cod: An application of generalized additive models to case-control data}. Environmental Health 2005, 4:11. 

Webster T, Vieira V, Weinberg J, Aschengrau A. \href{http://www.ij-healthgeographics.com/content/5/1/26}{Method for mapping population-based case-control studies using Generalized Additive Models}. International Journal of Health Geographics 2006, 5:26.

Young RL, Weinberg J, Vieira V, Ozonoff A, Webster TF. \href{http://www.ij-healthgeographics.com/content/9/1/37}{A power comparison of generalized additive models and the spatial scan statistic in a case-control setting}. International Journal of Health Geographics 2010, 9:37.

\url{http://www.busrp.org/projects/project2.html}

}

\author{
Veronica Vieira, Scott Bartell, and Robin Bliss

Send bug reports to \email{sbartell@uci.edu}.
}

\section{Warning }{Permutation tests are computationally intensive, often requiring several hours or more.} 

\seealso{
\code{\link{predgrid}}, 
\code{\link{optspan}}, 
\code{\link{colormap}}, 
\code{\link{readShapePoly}}.
}
\examples{
\donttest{
# Load base map in SpatialPolygonsDataFrame format
# This map was read from ESRI shapefiles using the readShapePoly function
data(MAmap)	
# Load data and create grid on base map
data(MAdata)							
gamgrid <- predgrid(MAdata, MAmap)    # requires PBSmapping package		
# Fit crude logistic GAM to the MA data using span size of 50\%
# and predict odds ratios for every point on gamgrid
fit1 <- modgam(MAdata, gamgrid, m="crude", sp=0.5)  
# Summary statistics for pointwise crude odds ratios 
summary(fit1$OR)			
# Summary stats for pointwise crude log odds (linear predictor) 
summary(fit1$fit)			

# fit adjusted GAM using span size of 50\%, 
# including a (too small) conditional permutation test
fit2 <- modgam(MAdata, gamgrid, permute=25, m="adjusted", sp=0.5)  

# Detailed example with a continuous outcome variable, map projections, 
# and data trimming:  investigating tweet times by geolocation
# NOTE: this example requires the maps, mapproj, and PBSmapping packages
# Convert base map and beer tweet data locations to US Albers projection 
# for better distance estimates than using (lat,long) as (X,Y) coords 
if(require(maps) & require(mapproj) & require(PBSmapping)) {
	USmap <- map("state",projection="albers",parameters=c(29.5,45.5),
		plot=FALSE,fill=TRUE,col="transparent")
	data(beertweets)
	# Reuse last map projection to convert data coordinates	
	XY <- mapproject(beertweets$longitude,beertweets$latitude)[1:2]  
	jtime <- julian(beertweets$time)
	# Convert tweet dates and times to time of day (24-hour)
	tweettime <- as.numeric(jtime-trunc(jtime))*24
	beerproj <- data.frame(tweettime, XY[1], XY[2], beertweets$beer)			 
	# Generate grid on the US map, trimmed to range of beer data
	USgrid <- predgrid(beerproj, USmap)						
	# Fit adjusted model--adjusting for beer indicator variable 
	fit3 <- modgam(beerproj, USgrid, family=gaussian, m="adjusted", sp=0.2)	
	# Get summary statistics for predicted tweet times across grid points 
	summary(fit3$fit)	
}

# Smoothing for two-dimensional chemical exposure instead of geolocation
# case status in 1st column, mercury and selenium in 2nd and 3rd columns   
ma2 <- MAdata[,c(1,5:6)]  
expgrid <- predgrid(ma2)
fit4 <- modgam(ma2,expgrid,sp=.5,m="crude") 
summary(fit4$OR)
# plot the results, with mercury on the X axis and selenium on the Y axis
colormap(fit4, contours="response", arrow=FALSE) 
}
}

\keyword{misc}
\keyword{smooth}
