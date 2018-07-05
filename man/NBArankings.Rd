\name{NBArankings}
\alias{NBArankings}
\docType{data}
\title{
A real example of rankings of NBA teams
}
\description{
This example is about aggregating the multiple rankings of NBA teams and was studied by Deng et al. (2014). They collected 34 rankings, including 6 professional rankings and 28 amateur rankings, for the 30 NBA teams in the 2011-2012 season. For the missing items in the partial rankings, we use number 0 to denote them.
}
\usage{data("NBArankings")}
\format{
  A data frame with 30 observations on the following 34 variables.
  \describe{
    \item{\code{V1}}{a factor with levels \code{76ers} \code{Bobcats} \code{Bucks} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Clippers} \code{Grizzlies} \code{Hawks} \code{Heat} \code{Hornets} \code{Jazz} \code{Kings} \code{Knicks} \code{Lakers} \code{Magic} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Pacers} \code{Pistons} \code{Raptors} \code{Rockets} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves} \code{TrailBlazers} \code{Warriors} \code{Wizards}}
    \item{\code{V2}}{a factor with levels \code{76ers} \code{Bobcats} \code{Bucks} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Clippers} \code{Grizzlies} \code{Hawks} \code{Heat} \code{Hornets} \code{Jazz} \code{Kings} \code{Knicks} \code{Lakers} \code{Magic} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Pacers} \code{Pistons} \code{Raptors} \code{Rockets} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves} \code{TrailBlazers} \code{Warriors} \code{Wizards}}
    \item{\code{V3}}{a factor with levels \code{76ers} \code{Bobcats} \code{Bucks} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Clippers} \code{Grizzlies} \code{Hawks} \code{Heat} \code{Hornets} \code{Jazz} \code{Kings} \code{Knicks} \code{Lakers} \code{Magic} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Pacers} \code{Pistons} \code{Raptors} \code{Rockets} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves} \code{TrailBlazers} \code{Warriors} \code{Wizards}}
    \item{\code{V4}}{a factor with levels \code{76ers} \code{Bobcats} \code{Bucks} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Clippers} \code{Grizzlies} \code{Hawks} \code{Heat} \code{Hornets} \code{Jazz} \code{Kings} \code{Knicks} \code{Lakers} \code{Magic} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Pacers} \code{Pistons} \code{Raptors} \code{Rockets} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves} \code{TrailBlazers} \code{Warriors} \code{Wizards}}
    \item{\code{V5}}{a factor with levels \code{76ers} \code{Bobcats} \code{Bucks} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Clippers} \code{Grizzlies} \code{Hawks} \code{Heat} \code{Hornets} \code{Jazz} \code{Kings} \code{Knicks} \code{Lakers} \code{Magic} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Pacers} \code{Pistons} \code{Raptors} \code{Rockets} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves} \code{TrailBlazers} \code{Warriors} \code{Wizards}}
    \item{\code{V6}}{a factor with levels \code{76ers} \code{Bobcats} \code{Bucks} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Clippers} \code{Grizzlies} \code{Hawks} \code{Heat} \code{Hornets} \code{Jazz} \code{Kings} \code{Knicks} \code{Lakers} \code{Magic} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Pacers} \code{Pistons} \code{Raptors} \code{Rockets} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves} \code{TrailBlazers} \code{Warriors} \code{Wizards}}
    \item{\code{V7}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Hawks} \code{Heat} \code{Lakers} \code{Pacers} \code{Spurs} \code{Thunder}}
    \item{\code{V8}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Clippers} \code{Heat} \code{Knicks} \code{Lakers} \code{Spurs} \code{Thunder}}
    \item{\code{V9}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Heat} \code{Knicks} \code{Lakers} \code{Mavericks} \code{Spurs} \code{Thunder}}
    \item{\code{V10}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Clippers} \code{Heat} \code{Lakers} \code{Mavericks} \code{Spurs} \code{Thunder}}
    \item{\code{V11}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Heat} \code{Knicks} \code{Lakers} \code{Nuggets} \code{Warriors} \code{Wizards}}
    \item{\code{V12}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Clippers} \code{Heat} \code{Lakers} \code{Mavericks} \code{Spurs} \code{Thunder}}
    \item{\code{V13}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Hornets} \code{Jazz} \code{Kings} \code{Lakers} \code{Magic} \code{Rockets}}
    \item{\code{V14}}{a factor with levels \code{0} \code{76ers} \code{Celtics} \code{Heat} \code{Kings} \code{Lakers} \code{Rockets} \code{Spurs} \code{Suns}}
    \item{\code{V15}}{a factor with levels \code{0} \code{Bulls} \code{Celtics} \code{Heat} \code{Lakers} \code{Mavericks} \code{Rockets} \code{Spurs} \code{Thunder}}
    \item{\code{V16}}{a factor with levels \code{0} \code{Celtics} \code{Hawks} \code{Heat} \code{Lakers} \code{Mavericks} \code{Raptors} \code{Spurs} \code{Thunder}}
    \item{\code{V17}}{a factor with levels \code{0} \code{76ers} \code{Celtics} \code{Heat} \code{Knicks} \code{Lakers} \code{Mavericks} \code{Nets} \code{Thunder}}
    \item{\code{V18}}{a factor with levels \code{0} \code{76ers} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Heat} \code{Lakers} \code{Mavericks} \code{Thunder}}
    \item{\code{V19}}{a factor with levels \code{0} \code{Bulls} \code{Heat} \code{Kings} \code{Lakers} \code{Rockets} \code{Spurs} \code{Suns} \code{Warriors}}
    \item{\code{V20}}{a factor with levels \code{0} \code{Bucks} \code{Celtics} \code{Heat} \code{Lakers} \code{Magic} \code{Mavericks} \code{Rockets} \code{Suns}}
    \item{\code{V21}}{a factor with levels \code{0} \code{Celtics} \code{Heat} \code{Kings} \code{Lakers} \code{Mavericks} \code{Spurs} \code{Suns} \code{Timberwolves}}
    \item{\code{V22}}{a factor with levels \code{0} \code{Celtics} \code{Heat} \code{Kings} \code{Lakers} \code{Spurs} \code{Suns} \code{Thunder} \code{Timberwolves}}
    \item{\code{V23}}{a factor with levels \code{0} \code{Bobcats} \code{Celtics} \code{Heat} \code{Lakers} \code{Mavericks} \code{Nuggets} \code{Spurs} \code{Suns}}
    \item{\code{V24}}{a factor with levels \code{0} \code{76ers} \code{Heat} \code{Knicks} \code{Lakers} \code{Pistons} \code{Rockets} \code{Spurs} \code{Wizards}}
    \item{\code{V25}}{a factor with levels \code{0} \code{76ers} \code{Celtics} \code{Hawks} \code{Heat} \code{Knicks} \code{Lakers} \code{Magic} \code{Thunder}}
    \item{\code{V26}}{a factor with levels \code{0} \code{Bulls} \code{Cavaliers} \code{Celtics} \code{Hawks} \code{Heat} \code{Knicks} \code{Lakers} \code{Rockets}}
    \item{\code{V27}}{a factor with levels \code{0} \code{76ers} \code{Clippers} \code{Lakers} \code{Magic} \code{Mavericks} \code{Pacers} \code{Raptors} \code{Warriors}}
    \item{\code{V28}}{a factor with levels \code{0} \code{76ers} \code{Bulls} \code{Celtics} \code{Heat} \code{Lakers} \code{Pistons} \code{Rockets} \code{Wizards}}
    \item{\code{V29}}{a factor with levels \code{0} \code{76ers} \code{Bulls} \code{Grizzlies} \code{Hawks} \code{Kings} \code{Knicks} \code{Nets} \code{Timberwolves}}
    \item{\code{V30}}{a factor with levels \code{0} \code{76ers} \code{Bucks} \code{Bulls} \code{Knicks} \code{Raptors} \code{Rockets} \code{Thunder} \code{Timberwolves}}
    \item{\code{V31}}{a factor with levels \code{0} \code{76ers} \code{Heat} \code{Lakers} \code{Magic} \code{Mavericks} \code{Pacers} \code{Pistons} \code{Suns}}
    \item{\code{V32}}{a factor with levels \code{0} \code{76ers} \code{Bulls} \code{Celtics} \code{Heat} \code{Knicks} \code{Lakers} \code{Magic} \code{Pacers}}
    \item{\code{V33}}{a factor with levels \code{0} \code{Clippers} \code{Heat} \code{Knicks} \code{Lakers} \code{Mavericks} \code{Nets} \code{Nuggets} \code{Wizards}}
    \item{\code{V34}}{a factor with levels \code{0} \code{Bulls} \code{Hawks} \code{Heat} \code{Jazz} \code{Knicks} \code{Nets} \code{Rockets} \code{Timberwolves}}
  }
}
\references{
Deng, K., Han, S., Li, K. J. and Liu, J. S. (2014). Bayesian aggregation of order-based rank data, Journal of the American Statistical Association 109(507): 1023-1039.
}
\examples{
data(NBArankings)
dim(NBArankings)
}
\keyword{datasets}
