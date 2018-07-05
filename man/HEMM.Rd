\name{HEMM}
\alias{HEMM}
\title{ A hierarchical extended Mallows model for aggregating multiple ranking lists
}
\description{
It uses the hierarchical extended Mallows model to aggregate multiple full/partial ranking lists.
}
\usage{
HEMM(rankings, num.kappa, is.kappa.ranker, initial.method, it.max)
}
\arguments{
  \item{rankings}{
A n by m matrix, with each column representing a ranking list, which ranks the items from the most preferred to the least preferred. For missing items, use 0 to denote them.
}
  \item{num.kappa}{the number of over-correlated ranking groups}
  \item{is.kappa.ranker}{ a list of over-correlated ranking groups, with the k-th element denoting the column numbers of the rankings that belong to the k-th group}
  \item{initial.method}{
the method for initializing the value of pi0, with four options: mean, median, geometric and random (the mean of three randomly sampled ranking lists). By default, initial.method="mean".
}
  \item{it.max}{
the maximum number of iterations. By default, it.max=20.
}
}
\value{
  \item{op.phi}{optimal value of phi}
  \item{op.phi1}{optimal value of phi1, the phi value in over-correlated ranking groups}
  \item{op.omega}{optimal value of omega}
  \item{op.alpha}{optimal value of alpha}
  \item{op.pi0}{optimal value of pi0, ranking the items from the most preferred to the least preferred}
  \item{op.kappa}{optimal value of kappa, denoting the items from the most preferred to the least preferred}
  \item{max.logL}{maximum value of log-likelihood}
}
\references{
An extended Mallows model for ranked data aggregation
}
\author{
Han Li, Minxuan Xu, Jun S. Liu and Xiaodan Fan
}
\examples{
data(simu3)
res=HEMM(rankings = simu3, num.kappa = 2, is.kappa.ranker = list(1:5, 6:10),
    initial.method = "mean", it.max = 20)
res$op.phi
res$op.phi1
res$op.omega
res$op.pi0

data(NBArankings)
res=HEMM(rankings = NBArankings, num.kappa = 1, is.kappa.ranker = list(1:6),
    initial.method = "mean", it.max = 20)
res$op.omega
res$op.pi0
res$op.kappa

}
