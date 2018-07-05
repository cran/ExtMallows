\name{simu1}
\alias{simu1}
\docType{data}
\title{
Simulation data 1
}
\description{
This data set is simulated as described in the Simulation Study 1 of the reference. It is a 30 by 6 data frame, representing 6 independent top-30 partial rankings.
}
\usage{data("simu1")}
\format{
  A data frame with 30 observations on the following 6 variables.
  \describe{
    \item{\code{V1}}{a numeric vector}
    \item{\code{V2}}{a numeric vector}
    \item{\code{V3}}{a numeric vector}
    \item{\code{V4}}{a numeric vector}
    \item{\code{V5}}{a numeric vector}
    \item{\code{V6}}{a numeric vector}
  }
}
\references{
An extended Mallows model for ranked data aggregation
}
\examples{
data(simu1)
dim(simu1)
}
\keyword{datasets}
