\name{print.amias}
\alias{print.amias}
\alias{print.samias}
\title{Print an amias or a samias object}
\description{
  Print a summary of the fitted \code{"amias"} or \code{"samias"} object.
}
\usage{
\method{print}{amias}(x, ...)
\method{print}{samias}(x, digits = max(3, getOption("digits") - 3), ...)
  
}
\arguments{
  \item{x}{A list with class attribute 'amias' or 'samias'.}
  \item{digits}{significant digits in printout}
  \item{\dots}{Additional print arguments.}
}

\seealso{
  \code{\link{amias}}, \code{\link{samias}}.
}
\details{
 For an \code{"amias"} object, the call that produced the object \code{x} is printed, followed by a summary of the output.

  For a \code{"samias"} object, the call that produced the object \code{x} is printed, followed by a three-column matrix with columns \code{k}, \code{Df} and \code{\%BIC}. The \code{k} column is the number of knots and the \code{Df} column is the degree of freedom defined as \eqn{Df = k+q+1}. \code{BIC} is the BIC value.
}
\value{
A summary or a three-column matrix is silently returned.
}

\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.
  
  Maintainer: Canhong Wen <wench@ustc.edu.cn>
}
\references{
  Wen, C., Zhu, J., Wang, X., and Zhang, A. (2018) \emph{L0 trend filtering}, technique report.
}


\examples{
  ##----- A toy example -------
  set.seed(0)
  n <- 100
  x = seq(1/n, 1,length.out = n)
  y0 = 0*x; y0[x>0.5] = 1
  y = y0 + rnorm(n, sd = 0.1)
  
  # For 'amias' object
  fit <- amias(y, k = 1) 
  print(fit)
  
  # For 'samias' object
  fit <- samias(y, kmax = 5)
  print(fit)
  
}
