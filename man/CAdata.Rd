\name{CAdata}
\alias{CAdata}
\docType{data}
\title{
Deidentified survival data for California
}
\description{
Survival time, censoring status, and geolocations (jittered to preserve anonymity) for 5000 ovarian cancer patients in California, using a Lambert projection (in meters).  \code{\link{CAmap}} is a map of California using the same projection.}

\usage{data(CAdata)}
\format{
  A data frame with 5000 observations on the following 6 variables.
  \describe{
    \item{\code{time}}{Survival time}
    \item{\code{event}}{Censoring status}
    \item{\code{X}}{Projected X coordinate}
    \item{\code{Y}}{Projected Y coordinate}
    \item{\code{AGE}}{Patient age}
    \item{\code{INS}}{Insurance type: \code{Mng} (Managed Care: managed care, HMO, PPO and other private insurance), \code{Mcr} (Medicare), \code{Mcd} (Medicaid), \code{Oth}: (Other Insurance: military, county-funded), \code{Uni} (Not Insured, self-pay) and \code{Unk} (Unknown)}
  }
}

\examples{
data(CAdata)
summary(CAdata)
plot(CAdata$X,CAdata$Y) 
}

\keyword{datasets}
