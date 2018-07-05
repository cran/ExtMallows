\name{EMM}
\alias{EMM}
\title{An extended Mallows model for aggregating multiple ranking lists
}
\description{
It uses the extended Mallows model to aggregate multiple full/partial ranking lists.
}
\usage{
EMM(rankings, initial.method, it.max)
}
\arguments{
  \item{rankings}{
A n by m matrix, with each column representing a ranking list, which ranks the items from the most preferred to the least preferred. For missing items, use 0 to denote them.
}
  \item{initial.method}{
the method for initializing the value of pi0, with four options: mean, median, geometric and random (the mean of three randomly sampled ranking lists). By default, initial.method="mean".
}
  \item{it.max}{
the maximum number of iterations. By default, it.max=20.
}
}
\value{
  \item{op.phi}{optimal value of phi}
  \item{op.omega}{optimal value of omega}
  \item{op.alpha}{optimal value of alpha}
  \item{op.pi0}{optimal value of pi0, ranking the items from the most preferred to the least preferred}
  \item{max.logL}{maximum value of log-likelihood}
}
\references{
An extended Mallows model for ranked data aggregation
}
\author{
Han Li, Minxuan Xu, Jun S. Liu and Xiaodan Fan
}
\examples{
data(simu1)
res=EMM(rankings = simu1, initial.method = "mean", it.max = 20)
res$op.phi
res$op.omega
res$op.pi0

}
