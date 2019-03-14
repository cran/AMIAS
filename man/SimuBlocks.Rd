\name{SimuBlocks}
\alias{SimuBlocks}
\alias{SimuWave}
\alias{SimuDoppler}
\title{Generate simulated data}
\description{
  Generate data for simulations.
}
\usage{
SimuBlocks(n, sigma = 0.1, seed = NA)
SimuWave(n, sigma = 0.1, seed = NA)
SimuDoppler(n, sigma = 0.1, seed = NA)

}
\arguments{
  \item{n}{The number of observations. }
  \item{sigma}{A positive value to specify the standard deviation of the Gaussian noise.}
  \item{seed}{a single value to specify the seed used for generating data. See \code{\link{set.seed}} for details.}
}

\details{
Let \eqn{x = i/n, i=1,\dots,n}. Assume that \eqn{y_x} follows an underlying trend \eqn{y0_x} with additive noise \eqn{\epsilon_x}, i.e., 
\deqn{y_x = y0_x + \epsilon_x,}
where \eqn{y0_x} is a smooth curve representing relatively long-term movements and \eqn{\epsilon_x} is often assumed to follow i.i.d. \eqn{N(0, \sigma^2)}.
}
\value{
A list with the following components:
  \item{y}{The simulated data of length \eqn{n}.}
  \item{x}{The locations of \code{y}, of length \eqn{n}.}
  \item{y0}{The true signal of length \eqn{n}.}
  \item{tau}{Locations of the true knots.}
  \item{SetA}{Indexs of the true knots.}
}


\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.
  
  Maintainer: Canhong Wen <wench@ustc.edu.cn>
}
\references{
  Wen, C., Zhu, J., Wang, X., and Zhang, A. (2019) \emph{L0 trend filtering}, technique report.
}

\examples{
  ##------ Piecewise constant trend filtering example in Wen et al.(2018).-----
  set.seed(1)
  data <- SimuBlocks(2048)
  fit <- samias(data$y, kmax = 15)           # With default input argument
  plot(fit, type="coef", main = "Blocks")    # Plot the optimal estimate
  lines(data$x, data$y0, type="s")           # Add the true signal for reference
  plot(fit, type= "vpath", main = "Blocks")  # Plot the solution path

}
