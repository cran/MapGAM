\name{modgam}
\alias{modgam}
\title{
Fit a Spatial Logistic Generalized Additive Model (GAM)}
\description{
Calculates crude or adjusted odds ratios on a user-supplied grid for spatial analysis using generalized additive model with a LOESS smooth for geolocation.  Includes an optional permutation test for global and local tests of the null hypothesis that spatial location is not associated with the outcome.}
\usage{
modgam(rdata, rgrid, permute = 0, m = "adjusted", sp = NULL, keep = F, verbose = T, \dots)
}
\arguments{
  \item{rdata}{
Data set (required). The data must be structured so that the outcome (coded as 0s and 1s) is in the 1st column and X and Y location values are in the 2nd and 3rd columns.  If additional columns are provided and m = "adjusted", then the additional columns will be entered as linear predictors in the GAM model. 
}
  \item{rgrid}{
A data frame containing the values at which to generate odds ratio predictions (required).  X and Y location values must be in the 1st and 2nd columns.  Additional covariates values for predictions may be provided in additional columns of rgrid.  If m = "adjusted" and rdata includes covariates not present in rgrid, then for each missing covariate the median value will be taken from rdata.  The \code{predgrid} function can be used to supply an appropriate rgrid for rdata, with the grid clipped according to any specified base map (see examples below).   
}
  \item{permute}{
The number of permutations of the data set for testing the significance of location. permute = 0 (default) produces no permutation tests.  If permute > 0, locations are randomly permuted in order to simulate the distribution of results under the null hypothesis that location is not associated with the outcome.  This is a conditional permutation test: the span size is held fixed throughout all permutations even when sp = NULL (see Young et al., 2012).  1000 permutations are recommended for reasonable accuracy of p-values.  WARNING: Because each permutation requires refitting the GAM, permutation tests can be quite slow.    
}
  \item{m}{
Model type for the GAM.  Options are "crude" or "adjusted" (default).  If "crude", only the spatial smooth term for location is included as a predictor in the model (i.e., "crude" is a synonym for "unadjusted").  If "adjusted", all covariates in the data set (columns >= 4) are included in the model.     
}
  \item{sp}{
Span size for the LOESS smooth of spatial location.  If sp = NULL (default) then an optimal span size will be determined using the optspan function.  
}
  \item{keep}{
Logical value indicating whether to store and return the pointwise odds ratios for every permuted data set.  These values aren't necessary for mapping or cluster identification, and storing them slows the permutation tests, so the default is keep = FALSE.
}  
  \item{verbose}{
Logical value indicating whether to print the GAM model statement, the percentile rank for the global deviance statistic in the permutation test, and the progress of the permutation test (report completion of every 10 permutations).  The default is verbose = TRUE.      
}
  \item{\dots}{
Further arguments to be passed to the \code{gam} function (e.g., weights).  	
}
}

\details{
The model used to fit the data is a generalized additive model with a LOESS smooth for geolocation (Hastie and Tibshirani, 1990; Kelsall and Diggle, 1998; Webster et al., 2006):
\deqn{\ln \left(\frac{\pi_{i}}{1-\pi_{i}}\right) = S(x_{i},y_{i}) + \mathbf{Z_{i}} \boldsymbol{\beta}}{
	ln(\pi / (1-\pi)) = S(x,y) + Z \beta}	
where \eqn{\pi_{i}}{\pi} is the probability that the outcome is 1 for participant \emph{i}, \eqn{x_{i}}{x} and \eqn{y_{i}}{y} are spatial coordinates for participant \emph{i} (i.e., projected distance east and north, respectively, from an arbitrarily defined origin location), \eqn{S(.,.)} is a 2-dimensional smoothing function (currently LOESS), \eqn{\mathbf{Z_{i}}}{Z} is a row vector of covariate values for participant \emph{i}, and \eqn{\boldsymbol{\beta}}{\beta} is a vector of unknown regression parameters (including an intercept).  
When a permutation test is requested, for each permutation the geolocations (paired x and y coordinates) in the data set are randomly reassigned to participants, consistent with the null hypothesis that geolocation is not associated with risk for the outcome.  Note that this procedure preserves associations between other covariates and risk of the outcome, so the permutations only test the signficance of the smoothed geolocation.  See the references for more details.              
}

\value{
\item{grid}{A data frame with X and Y locations from rgrid}
\item{OR}{The estimated odds ratio for each point on the grid, from the spatially smoothed GAM model.  The referent is logistic model without spatial location.}
\item{m}{Whether the GAM was unadjusted (only spatial location as a predictor) or adjusted (included other covariates in addition to spatial location).  "Crude" is a synonym for "unadjusted."}
\item{span}{The span size used for the LOESS smooth of the effects of location}
\item{gamobj}{The GAM model object from the fit to the data}
If permute > 0 then the following values are also provided:
\item{global}{The p-value based on the deviance statistic comparing models with and without spatial locations, using the distribution of deviance statistics from the permuted data sets to represent the null distribution.  For a test of H0: spatial location is unassociated with risk of the outcome (adjusting for any covariates if m = "adjusted"), reject H0 if the percentile rank is below alpha.   WARNING: \code{modgam} uses a conditional permutation test which produces inflated type I error rates; Young et al. (2012) recommend using alpha=0.025 to limit the type I error rate to approximately 5\%.} 
\item{pointwise}{For each point on the grid, the percentile rank of the local log odds for the model compared to the local log odds distribution from the permuted data sets.  This result is needed to define clusters of the outcome (spatial regions with statistically significant risks--unusually high or low risks.)}
If permute > 0 and keep = T then the following values are also provided:
\item{permutations}{A matrix containing null permuted odds ratio estimates for each point on the grid, with the results for each permutation in a separate column.}
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

\section{Warning }{Permutation tests are computationally intensive, often requiring several hours for a large data set.} 

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
gamgrid <- predgrid(MAdata, MAmap)		
# Fit crude GAM model to the MA data using span size of 50\%
# and predict odds ratios for every point on gamgrid
fit1 <- modgam(MAdata, gamgrid, m="crude", sp=0.5)  
# Get summary statistics for pointwise crude odds ratios 
summary(fit1$OR)				

# fit adjusted GAM model to the MA data using span size of 50\%, and run a (too small) permutation test
fit2 <- modgam(MAdata, gamgrid, permute=25, m="adjusted", sp=0.5)  
}
}

\keyword{misc}
\keyword{smooth}
