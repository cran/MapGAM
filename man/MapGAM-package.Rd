\name{MapGAM-package}
\alias{MapGAM-package}
\alias{MapGAM}
\docType{package}
\title{
Mapping Smoothed Odds Ratios from Individual-Level Spatial Data
}
\description{
Contains functions for mapping risk surfaces using individual-level data such as case-control study data, using generalized additive models (GAMs) to smooth over two-dimensional location while adjusting for any confounding variables using methods described by Kelsall and Diggle (1998) and Webster at al. (2006).  Includes convenient functions for mapping, efficient control sampling, and permutation tests for the null hypothesis that spatial location is not associated with risk (adjusting for any confounders).     
}
\details{
\tabular{ll}{
Package: \tab MapGAM\cr
Type: \tab Package\cr
Version: \tab 0.6-1\cr
Date: \tab 2012-12-23\cr
License: \tab GPL-3\cr
}
Typical applications will start with the predgrid function to create a regular grid of points within the study area, with optional map boundaries (e.g., a country, state, or regional map).  Crude or adjusted odds ratios are then estimated at each grid point using the modgam function to smooth by geolocation.  Finally, the odds ratios (and optionally, "clusters"--areas of signficantly increased or decreased risk determined via permutation tests) are plotted using the colormap function.  The trimdata and sampcont functions can be used to restrict data to those within map boundaries and to conduct simple or spatiotemporal stratified sampling from eligible controls.  The optspan function can be used to find an optimal span size for the LOESS smoother used by the modgam function; it is automatically used within the modgam function when the span size is not provided by the user.}

\author{
Veronica Vieira, Scott Bartell, and Robin Bliss

Send bug reports to \email{sbartell@uci.edu}.
}

\references{
Hastie TJ, Tibshirani RJ.  Generalized Additive Models.  (Chapman & Hall/CRC Monographs on Statistics & Applied Probability, Boca Raton, Florida, 1990).  

Kelsall J, Diggle P. Spatial variation in risk of disease: a nonparametric binary regression approach. J Roy Stat Soc C-App 1998, 47:559-573. 	
	
Vieira V, Webster T, Weinberg J, Aschengrau A, Ozonoff D. \href{http://www.ehjournal.net/content/4/1/11}{Spatial analysis of lung, colorectal, and breast cancer on Cape Cod: An application of generalized additive models to case-control data}. Environmental Health 2005, 4:11. 

Webster T, Vieira V, Weinberg J, Aschengrau A. \href{http://www.ij-healthgeographics.com/content/5/1/26}{Method for mapping population-based case-control studies using Generalized Additive Models}. International Journal of Health Geographics 2006, 5:26.

Young RL, Weinberg J, Vieira V, Ozonoff A, Webster TF. \href{http://www.ij-healthgeographics.com/content/9/1/37}{A power comparison of generalized additive models and the spatial scan statistic in a case-control setting}. International Journal of Health Geographics 2010, 9:37.

\url{http://www.busrp.org/projects/project2.html}
}

\keyword{ package }
\keyword{ smooth }
\keyword{ datasets }
\keyword{ hplot }
\keyword{ misc }

\seealso{
\code{\link{trimdata}}, 
\code{\link{sampcont}}, 
\code{\link{predgrid}}, 
\code{\link{optspan}}, 
\code{\link{modgam}}, 
\code{\link{colormap}}.
}

\examples{
# Load synthetic data and a preformatted base map
data(MAmap)		
data(MAdata)
# Create a grid on the base map (PBSmapping package recommended)	
if(require(PBSmapping)) gamgrid <- predgrid(MAdata, MAmap) else
	 gamgrid <- predgrid(MAdata) 		
# Fit a GAM with a smooth term for spatial location
fit1 <- modgam(MAdata, gamgrid, m="crude", sp=0.5)   
\donttest{
# Display odds ratio estimates on the base map
colormap(fit1, MAmap)								
}
#### See colormap and modgam help files for more examples 
}
