\name{plot.AMIAS.1d}
\alias{plot.AMIAS.1d}
\title{plot observe data with fitted value for AMIAS model}
\description{
  Similar to other plot methods, this functions plot fitted lines with the given sequence from a fitted \code{"AMIAS.1d"} object.
}
\usage{
  \method{plot}{AMIAS.1d}(x, s=length(x$df), betatype="l", betacolor="blue", betalwd=2, ...)
}
\arguments{
  \item{x}{A list with class attribute 'AMIAS.1d'.}
  \item{s}{If the AMIAS object contain a T sequence to select adatively, the s-th model would be plot, default is the selected model.}
  \item{betatype}{Character indicating the type of plotting; actually any of the types as in \code{\link[graphics]{plot.default}}.}
  \item{betacolor}{Color of the fitted line.}
  \item{betalwd}{Wide of the fitted line.}
  \item{\dots}{Other arguments,see \code{\link[graphics]{plot.default}} for details.}
}

\seealso{
  \code{\link{AMIAS}} \code{\link[graphics]{plot.default}} \code{\link[graphics]{lines.default}}.
}

\examples{
  library(AMIAS)

  set.seed(12580)
  n <- 100
  sigma <- 0.3
  y0 <- rep(0,n)
  y0[10:15] <- 2
  y0[40:60] <- -1
  y0[80:82] <- 4
  y <- y0 + sigma*rnorm(n)
  y[80:82] <- y0[80:82] + sigma*rnorm(3)
  object <- AMIAS(y, h=1, outer_itermax=5)
  plot(object)
}
