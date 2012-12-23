\name{colormap}
\alias{colormap}
\title{
Maps Odds Ratios and Clusters for modgam Objects
\description{
Displays a color image map, including a legend, scale bar, and North arrow, showing crude or adjusted odds ratios for a grid of points.  Irregular grids are allowed.  Also draws contour lines for regions of increased or decreased odds of the outcome ("clusters"), if permutation ranks are provided.  Designed to display \code{modgam} objects but can be used with other model results if the modgamobj list contains the correct elements.      
}
\usage{
colormap(modgamobj, map = NULL, add=F, mapmin = NULL, mapmax = NULL, ptsize = 0.9, alpha = 0.05)
}
\arguments{
  \item{modgamobj}{
(Required) A list containing at least these two elements:  grid (a 2 column data frame with X and Y coordinates) and OR (a vector of odds ratios, with length equal to the number of grid points).  If the list contains the element pointwise (a vector of percentile ranks generated by permutation test) then those values will be used to generate contours.  The correct list format is provided by the modgam function (see examples).   
  \item{map}{
Can be used to map odds ratios on a base map from the \code{map} function in the \pkg{maps} package, or on a base map produced by the \code{readShapePoly} function in \pkg{maptools} package.  ReadShapePoly reads maps from external files in the ESRI shapefile format.  \code{map=NULL} produces a color image without any base map.  To add the color image to an existing plot use \code{add=T}.        
}
  \item{add}{
Use add=T to add the color map to an existing plot.  	
}
  \item{mapmin}{
The minimum odds ratio for the color scale legend 
}
  \item{mapmax}{
The maximum odds ratio for the color scale legend 
}
  \item{ptsize}{
The size of the points used to fill in colors on the map.  Increase to remove white space inside the map or decrease to remove color outside the map boundaries.  NOTE: white space can also be eliminated by increasing the grid size in \code{predgrid}, which is often preferable as it results in a higher resolution map.   
}
  \item{alpha}{
The nominal pointwise type I error rate.  If the modgam object includes pointwise permutation percentile ranks, contour lines are added at the percentiles alpha/2 and 1-alpha/2 from the null permutations, denoting clusters with statistically significant decreases or increases in odds.    
}
}

\value{
Produces an image map showing crude or adjusted odds ratios for the outcome.  If the base map is in \code{readShapePoly} format a scale bar is included.  The scale bar assumes that the X and Y coordinates are provided in meters.      
}

\author{
Veronica Vieira, Scott Bartell, and Robin Bliss

Send bug reports to \email{sbartell@uci.edu}.
}

\section{Warning}{Note that the contour lines use a pointwise nominal type I error rate of alpha; the chance of a type I error occurring for \strong{at least} one point on the map is much higher, typically approaching 100\% at alpha=0.05, because a spatial prediction grid generally contains many points.} 

\seealso{
  \code{\link{trimdata}}, 
  \code{\link{predgrid}}, 
  \code{\link{optspan}}, 
  \code{\link{modgam}}.
}

\examples{
\donttest{
### Load readShapePoly base map and data 
data(MAmap)		
data(MAdata) 
# create grid on base map
MAgrid <- predgrid(MAdata, MAmap)		
# fit crude GAM model to the MA data using span size of 50\%
fit1 <- modgam(MAdata, MAgrid, m="crude", sp=0.5)  
# Plot a map showing crude odds ratios
colormap(fit1, MAmap)					

#### A detailed example including map projections and data trimming
# NOTE: this example requires the maps package and mapproj package
# Convert base map and beer tweet data locations to US Albers projection 
# This results in better distance estimates than using (lat,long) as (X,Y) coords 
if(require(maps) & require(mapproj)) {
	USmap <- map("state",projection="albers",parameters=c(29.5,45.5),
				plot=FALSE,fill=TRUE,col="transparent")
	data(beertweets)
	case <- beertweets$beer
    # Reuse last map projection to convert data coordinates	
	XY <- mapproject(beertweets$longitude,beertweets$latitude)[1:2]  
	beerproj <- data.frame(case, XY[1], XY[2])			 
	# Generate grid on the US map, trimmed to range of beer data
	USgrid <- predgrid(beerproj, USmap)						
    # Fit unadjusted model--geolocation only
	fit2 <- modgam(beerproj, USgrid, m="unadjusted", sp=0.2)	
	colormap(fit2, USmap)								
	title(main="Beer Tweet Odds Ratios")	
}
}
}

\keyword{hplot}
\keyword{misc}
\keyword{smooth}