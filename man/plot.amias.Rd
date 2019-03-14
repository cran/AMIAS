\name{plot.amias}
\alias{plot.amias}
\alias{plot.samias}
\title{Plotting method for amias and samias objects}
\description{
   The function \code{plot.amias} plots the coefficient from an \code{"amias"} object with a user-specified number of knots as a function of the input positions (which are assumed to be evenly spaced if not specified). The function \code{plot.samias} produces either a profile plot of the solution path for the primal variable or the optimal fitted coefficients (which, recall, are determined by minimizing the Bayesian Information Criterion) from a fitted \code{"samias"} object.
}
\usage{
\method{plot}{amias}(x, add.knots = TRUE,...)
\method{plot}{samias}(x, type = c("coef", "vpath"), k, add.label = TRUE, add.knots = TRUE,...)

}
\arguments{
  \item{x}{A list with class attribute 'amias' or 'samias'.}
  \item{add.knots}{Whether to add the locations of the detected knots at the top of the figure. It only works when \code{type = "coef"} for the function \code{plot.samias}.}
  \item{type}{Types of the plotting figure, either "coef" or "vpath". For \code{type = "coef"} and missing \code{k}, we plot the optimal \eqn{l_0} fitted coefficients as well as the raw data.  For \code{type = "coef"} and a user-specified \code{k}, we plot the \eqn{l_0} fitted coefficients with number of knots \code{k} as well as the raw data. For \code{type = "vpath"}, we plot the solution path of the primal variable \eqn{v} versus the cardinality \eqn{k}, i.e., number of the detected knots. }
  \item{k}{The number of knots. Users can choose k from \eqn{1,2,\dots, kmax}; if none are specified, then it plots the optimal coefficient determined by miniming the BIC value.}
  \item{add.label}{Whether to add labels at the right side to indicate the locations of the detected knots. It only works when \code{type = "vpath"} for the function \code{plot.samias}.}
  \item{\dots}{Other arguments,see \code{\link[graphics]{plot.default}} for details.}
}

\seealso{
  \code{\link{amias}}, \code{\link{samias}}.
}

\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.
  
  Maintainer: Canhong Wen <wench@ustc.edu.cn>
}
\references{
  Wen, C., Zhu, J., Wang, X., and Zhang, A. (2019) \emph{L0 trend filtering}, technique report.
}


\examples{

  ##----- A toy example of piecewise constant signal -------
  set.seed(0)
  n <- 100
  x = seq(1/n, 1,length.out = n)
  y0 = 0*x; y0[x>0.5] = 1
  y = y0 + rnorm(n, sd = 0.1)
  
  # For 'amias' object
  fit <- amias(y, k = 1) 
  plot(fit)
  
  # For 'samias' object
  fit <- samias(y, kmax = 5)
  plot(fit, type = "coef", main = "Piecewise Constant")
  
  op <- par(mfrow=c(1,2))
  plot(fit, type= "coef", add.knots = FALSE, main = "Piecewise Constant")
  plot(fit, type = "vpath", main = "Piecewise Constant")
  par(op)
  
  ##----- A toy example of piecewise linear trend -------
  set.seed(0)
  y0 = 2*(0.5-x); y0[x>0.5] = 2*(x[x>0.5]-0.5)
  y = y0 + rnorm(n, sd = 0.1)
  
  # For 'amias' object
  fit <- amias(y, D_type = "tfq", q = 1, k = 1) 
  plot(fit, main = "Piecewise Linear")
  
  # For 'samias' object
  fit <- samias(y, D_type = "tfq", q = 1, kmax = 4)
  op <- par(mfrow=c(1,2))
  plot(fit, type = "coef", main = "Piecewise Linear")
  plot(fit, type = "vpath", main = "Piecewise Linear")
  par(op)
  
}
