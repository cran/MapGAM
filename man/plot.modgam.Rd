\name{plot.modgam}
\alias{plot.modgam}
\title{
Maps Predicted Values and Clusters for modgam Objects
}
\description{
Displays a color image map, including a legend, scale bar, and optional North arrow, showing crude or adjusted odds ratios (or linear predictors) for a grid of points.  Irregular grids are allowed.  Also draws contour lines for regions of signficantly increased or decreased values of the outcome variable ("clusters"), if permutation ranks are provided.  Designed to display \code{modgam} objects but can be used with other model results if the modgamobj list contains the correct elements.      
}
\usage{
\method{plot}{modgam}(x, map = NULL, exp = F, add = F, intervals=T, mapmin = NULL, 
       mapmax = NULL, border.gray = 0.3,
       contours=c("none","response","permrank","interval"),
       contours.drawlabels=FALSE, contours.lty=1, contours.lwd=1, 
       contours.levels, contours.labcex=0.7, arrow=T, axes=F, ptsize=0.9,
       alpha=0.05, mai, legend.name = "predicted values", legend.cex=1, \dots)
}
\arguments{
  \item{x}{
(Required) An object of class \code{modgam}, usually returned by the \code{modgam} function.  
}
  \item{map}{
Can be used to map predicted values on a base map from the \code{map} function in the \pkg{maps} package, or on a base map produced by the \code{readShapePoly} function in \pkg{maptools} package.  ReadShapePoly reads maps from external files in the ESRI shapefile format.  \code{map=NULL} produces a color image without any base map.
}
  \item{exp}{
If \code{exp=T}, \code{fit} will be displayed in exponential scale.  
}
  \item{add}{
Use \code{add=T} to add the color map to an existing plot.  This will often result in loss of the legend and scale, which are added outside of the normal
    map boundaries.  \code{add} is ignored when a map is provided using the 
    \code{map} argument.
}
  \item{intervals}{
If the \code{modgam} object contains confidence intervals, the lower and higher intervals will be plotted if \code{intervals=T}  
}
  \item{mapmin}{
The minimum value for the color scale legend. 
}
  \item{mapmax}{
The maximum value for the color scale legend. 
}
  \item{border.gray}{
Gray scale for the border of \code{map}, ranges from 0 to 1. It is white with \code{border.gray=1}, and black with \code{border.gray=0}. The defalt scale is 0.3.
}
  \item{contours}{
Use \code{contours="response"} to add contour lines for the predicted response, for example to draw isoboles for mixtures of exposures.  Use \code{contours="permrank"} to add contour lines for pointwise p-values computed from the permutation ranks, at \code{alpha/2} and \code{(1-alpha)/2}. use \code{contours = "intervals"} to add contour lines for areas where the confidence intervals do not include 0. The default is "none" which produces no contour lines.
}
  \item{contours.drawlabels}{
Use \code{contours.drawlabels = TRUE} add labels on contour lines.   
}
  \item{contours.lty}{
The line type for contour lines
}
  \item{contours.lwd}{
The width for contour lines.
} 
  \item{contours.levels}{
The levels for contour lines. When \code{exp=T}, the levels will also be used in an exponential scale.
} 
  \item{contours.labcex}{
cex for contour labelling. This is an absolute size, not a multiple of \code{par("cex")} 
}
  \item{arrow}{
Use arrow=T to add a North arrow to the map.  	
}
  \item{axes}{
Use axes=T to add axes to the map (useful for chemical mixture isoboles).
}
  \item{ptsize}{
The size of the points used to fill in colors on the map.  Increase to remove white space inside the map or decrease to remove color outside the map boundaries.  NOTE: white space can also be eliminated by increasing the grid size in \code{predgrid}, which is often preferable as it results in a higher resolution map.   
}
  \item{alpha}{
The nominal pointwise type I error rate; only used when \code{contours="permrank"}.  
}
  \item{mai}{
The margins of the plot. Details see \link{par},}
  \item{legend.name}{
The name of displayed for the legend bar.}
  \item{legend.cex}{
A numerical value giving the amount by which legend text should be magnified relative to the default.  
}
  \item{...}{
Other arguments to be passed to \link{colormap} function
  }
}

\value{
Produces an image map showing crude or adjusted linear predictors (or odds ratios).  If the base map is in \code{readShapePoly} format a scale bar is included.  The scale bar assumes that the X and Y coordinates are provided in meters.      
}

\author{
Scott Bartell, Lu Bai, Robin Bliss, and Veronica Vieira

Send bug reports to \email{sbartell@uci.edu}.
}

\section{Warning}{Note that the contour lines use a pointwise nominal type I error rate of alpha; the chance of a type I error occurring for \strong{at least} one point on the map is much higher, typically approaching 100\% at alpha=0.05, because a spatial prediction grid generally contains many points.} 

\seealso{
  \code{\link{trimdata}}, 
  \code{\link{predgrid}}, 
  \code{\link{optspan}}, 
  \code{\link{modgam}},
  \code{\link{colormap}}.
}

\examples{
\donttest{
data(MAmap)		
data(MAdata) 
# Create a grid on the base map (PBSmapping package recommended)
if(require(PBSmapping)) MAgrid <- predgrid(MAdata, map=MAmap) else 
MAgrid <- predgrid(MAdata)		
# fit crude GAM model to the MA data using span size of 50\%
fit1 <- modgam(data=MAdata, rgrid=MAgrid, m="crude", sp=0.5)  
# Plot a map showing crude odds ratios
plot(fit1, MAmap)					

#### A detailed example including map projections and data trimming
# NOTE: this example requires the maps, mapproj, and PBSmapping packages
# Convert base map and beer tweet data locations to US Albers projection 
# projected coords yield better distance estimates than (lat,long)
if(require(maps) & require(mapproj) & require(PBSmapping)) {
  USmap <- map("state",projection="albers",parameters=c(29.5,45.5),
               plot=FALSE,fill=TRUE,col="transparent")
  data(beertweets)
  case <- beertweets$beer
  # Reuse last map projection to convert data coordinates	
  XY <- mapproject(beertweets$longitude,beertweets$latitude)[1:2]  
  beerproj <- data.frame(case, XY[1], XY[2])			 
  # Generate grid on the US map, trimmed to range of beer data
  USgrid <- predgrid(beerproj, map=USmap)						
  # Fit unadjusted model--geolocation only
  fit2 <- modgam(data=beerproj, rgrid=USgrid, m="unadjusted", sp=0.2)	
  dev.new(width=7,height=5)
  plot(fit2, USmap)
  title(main="Beer Tweet Odds Ratios")	
}
}
}

\keyword{hplot}
\keyword{misc}
\keyword{smooth}
