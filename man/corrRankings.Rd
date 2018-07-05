\name{corrRankings}
\alias{corrRankings}
\title{p value for measuring the correlation of pairwise rankings
}
\description{
It caclulates the p values that measure the correlation of pariwise rankings.
}
\usage{
corrRankings(rankings)
}
\arguments{
  \item{rankings}{
A n by m data frame, with each column representing a ranking list, which ranks the items from the most preferred to the least preferred. For missing items, use 0 to denote them.
}
}

\value{
  \item{pair.pvalue}{a symmetric matrix of p values, with the (i,j)-th element denoting the p value of the i,j-th rankings.}
}
\references{
An extended Mallows model for ranked data aggregation
}
\author{
Han Li, Minxuan Xu, Jun S. Liu and Xiaodan Fan
}
\note{
Note that the input rankings should have at least 8 rankings. When constructing the samples of rescaled V distance for a given rank position, the number of samples should at least be 28 and the number of rankings that have items up to this position should account for at least 2/3 of the total number of rankings, otherwise the p value calculation stops at this position.
}
\examples{

data(simu3)
pvalue=corrRankings(rankings = simu3)

#threshold the p values

threshold=0.05
pvalue.trunc=ifelse(pvalue<=0.05, pvalue, 1)

#plot the p values

x=y=1:ncol(pvalue)
par(mfrow=c(1,2))
image(x, y, pvalue, xlab = NA, ylab = NA, sub = "rank coefficient")
image(x, y, pvalue.trunc, xlab = NA, ylab = NA, sub = "rank coefficient < 0.05")

}
