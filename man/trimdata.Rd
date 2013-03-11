\name{trimdata}
\alias{trimdata}
\title{
Trim a Data Set To Map Boundaries 
}
\description{
Takes a subset of a data frame: returns rows with geolocations inside the boundaries determined by a map.  If \code{rectangle=FALSE}, strict map boundaries are used.  If \code{rectangle=TRUE}, a rectangular boundary is determined based on the range of X and Y coordinates for the map, along with an optional buffer.   
}
\usage{
trimdata(rdata, map, Xcol=2, Ycol=3, rectangle=F, buffer=0.05)
}
\arguments{
  \item{rdata}{
Data set (required).  If \code{rdata} has only 2 columns then they are assumed to be the X and Y coordinates, in that order.  If \code{rdata} has more than 2 columns then identify the positions of the X and Y coordinates with \code{Xcol} and \code{Ycol}, respectively.
}
  \item{map}{
A map for clipping the grid, provided in a recognized format.  Supported classes include "map", the format produced by the maps package, and "SpatialPolygonsDataFrame", which can be produced from an ESRI shapefile using the readShapePoly function in the maptools package.
}
  \item{Xcol}{
The column number of rdata for the X coordinates of geolocation.  Only used if rdata has >2 columns.
  }
  \item{Ycol}{
The column number of rdata for the Y coordinates of geolocation.  Only used if rdata has >2 columns.
  }
  \item{rectangle}{
 If \code{rectangle=FALSE} (default), only data with geolocations strictly within the map boundaries are retained.  If \code{rectangle=TRUE}, only data with geolocations within a rectangular boundary are retained.  The rectangular boundary coordinates are determined using the ranges of X and Y coordinates for the map, along with an optional buffer.   
}
  \item{buffer}{
A fraction of the map range to add to either side of the rectangular boundary before trimming the data (default 5\%).  The buffer argument is ignored if \code{rectangle=FALSE}.  If a vector of length 2 is supplied, the first value as the buffer fraction for the X coordinate range and the second value is used as the buffer fraction for the Y coordinate range. 
}
}
\details{
Various functions from the PBSmapping package and maptools package are used to convert map formats, check whether the map and data are aligned, and clip the grid.  Be sure to use \code{fill=T} when using the \code{map} function to produce maps for \code{trimdata}: maps produced using the default \code{fill=F} argument do not work properly with this function.  If the map centroid is not in the range of the data a warning message is printed; this might indicate differing projections but can occur naturally when the data were not sampled from the entire extent of the map or when map boundaries are concave.      
}

\value{
A subset of the rdata data frame containing only those rows with geolocations within the specified boundaries.      
}

\author{
Veronica Vieira, Scott Bartell, and Robin Bliss

Send bug reports to \email{sbartell@uci.edu}.
}

\seealso{
\code{\link{predgrid}},
\code{\link{optspan}}, 
\code{\link{modgam}}, 
\code{\link{colormap}}.
}

\examples{
# These examples require the "PBSmapping" and "maps" packages
if (require(maps) & require(PBSmapping)) {
data(beertweets)
dim(beertweets)

### Trim data to US base map, and plot them
basemap1 <- map("usa", fill=TRUE, col="transparent")	
dUS <- trimdata(beertweets, basemap1)	 				 
# Plot tweet locations (beer tweets in red)
points(dUS$longitude, dUS$latitude, col=dUS$beer+1, cex=0.5)		

### Trim data to Texas base map, and plot them
basemap2 <- map("state", regions="texas", fill=TRUE, col="transparent") 
dTX <- trimdata(beertweets, basemap2) 								
# Plot tweet locations (beer tweets in red)
points(dTX$longitude, dTX$latitude, col=dTX$beer+1, cex=0.5) 
}
}
\keyword{misc}
