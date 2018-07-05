\name{simu2}
\alias{simu2}
\docType{data}
\title{
Simulation data 2
}
\description{
This data set is simulated as described in the Simulation Study 2 of the reference. It is a 40 by 6 data frame, representing 6 independent top-40 partial rankings.
}
\usage{data("simu2")}
\format{
  A data frame with 40 observations on the following 6 variables.
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
data(simu2)
dim(simu2)
}
\keyword{datasets}
