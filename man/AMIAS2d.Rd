\name{AMIAS2d}
\alias{AMIAS2d}
\title{Alternating Minimization Induced Active Set Algorithms(2d)}
\description{
  Extremely Efficient Procedures for Composite L0 Penalized Estimation(2d).
}
\usage{
  AMIAS2d(Mat, D_type = "fused.2d", composite = FALSE, W_type = "identity",
  k1=NULL, k2=1, D=NULL, W=NULL, T1=NULL, T2=min(10,nmw[1]-1), rho1=(dim1*dim2)^6, 
  rho2=(dim1*dim2)^6, h=5, tao=1, outer_itermax=20, select_max=min(20,nm[1]-1), 
  eps=0.1, iter_max=10, ...)
}
\arguments{
  \item{Mat}{a grapha matrix.}
  \item{D_type}{Specifcation of \eqn{D}, either one of "fused.2d" or "fused.tfk". when \eqn{D} is set, it would be "user".}
  \item{composite}{Whether to use composite L0 penalty.}
  \item{W_type}{Specifcation of \eqn{W}, either one of "identity", "fused.1d" or "fused.tfk". when \eqn{W} is set, it would be "user".}
  \item{k1}{Integer: used to generate \eqn{D} for polynomial trend filtering of order k(see \code{\link{genDtf1d}} for Details).}
  \item{k2}{Integer: used to generate \eqn{W} for polynomial trend filtering of order k(see \code{\link{genDtf1d}} for Details).}
  \item{D}{Matrix in the penalty item(see Details).}
  \item{W}{Matrix in the penalty item for composite L0 panalty(see Details).}
  \item{T1}{The number of change point to detect in \eqn{D\beta} penalty,default is NULL and it would be selected adatively.}
  \item{T2}{The number of change point to detect in \eqn{W\beta} penalty,default is 10.}
  \item{rho1}{The lagrange operator of the argumented lagrange form for split item in \eqn{D\beta}.}
  \item{rho2}{The lagrange operator of the argumented lagrange form for split item in \eqn{W\beta}.}
  \item{h}{Arguments for \code{\link{my.rollmean}}.}
  \item{tao}{Number: increment of the T sequence in adative function algorithm.}
  \item{outer_itermax}{The maximum iterations to select T adatively.}
  \item{select_max}{The maximum number of change point to detect.}
  \item{eps}{The threshold to stop the adative AMIAS algorithm, when \eqn{MSE} is smaller than \eqn{eps}.}
  \item{iter_max}{The maximum iterations to setect change point for given T.}
  \item{\dots}{Other arguments.}
}
\details{
  The Parameter of models implied by T is fit by primal dual active set 
  and this is the number of change point to detect in AMIAS. Note that the 
  model to fit is \deqn{ \beta_i + \sigma_i,} and the objective function 
  is \deqn{1/2*RSS + \lambda_1*penalty_1 + \lambda_2*penalty_2(\lambda_2 = 0 for 
  non-composite type),} where \deqn{panalty_1 = |D\beta|_0 , panalty_2 = |W\beta|_0.}

  If the parameter \code{T1} is not set, a sequence of \code{T1} would be fitted 
  by increment \code{tao} adatively. For details, see the AMIAS paper.
}
\value{
A list with class attribute 'AMIAS' and named components:
  \item{call}{The call that produced this object.}
  \item{Mat}{Observe grapha matrix.}
  \item{Mathat}{The fitting coefficients in the noparameter model.}
  \item{composite}{Whether to use composite L0 penalty.}
  \item{k1}{Integer: used to generate D for polynomial trend filtering of order k(see Details).}
  \item{k2}{Integer: used to generate W for polynomial trend filtering of order k(see Details).}
  \item{alpha}{The split variable of the argumented lagrange form in \eqn{D\beta}}
  \item{u}{The lagrange operator of the argumented lagrange form in \eqn{D\beta} for linear item.}
  \item{gamma}{The split variable of the argumented lagrange form in \eqn{W\beta}}
  \item{v}{The lagrange operator of the argumented lagrange form in \eqn{W\beta} for linear item.}
  \item{df}{Degree of freedom of the seleted model.}
  \item{T2}{The number of change point to detect in \eqn{W\beta} penalty.}
  \item{D_type}{Specifcation of \eqn{D}, default is  "fused.1d". when D is set, it would be "user".}
  \item{W_type}{Specifcation of \eqn{D}, default is  "identity". when W is set, it would be "user".}
  \item{rho1}{The lagrange operator of the argumented lagrange form in \eqn{D\beta} for split item.}
  \item{rho2}{The lagrange operator of the argumented lagrange form in \eqn{W\beta} for split item.}
  \item{tao}{Number: increment of the T sequence in adative function algorithm.}
  \item{eps}{The threshold to stop the adative AMIAS algorithm, when \code{MSE} is smaller than \code{eps}.}
  \item{iter}{The iterations to setect change point for given \code{T}}
  \item{iter_max}{The maximum iterations to setect change point for given \code{T}}
  \item{outer_itermax}{The maximum iterations to select \code{T} adatively.}
  \item{select_max}{The maximum number of change point to detect.}
}

\references{
  Wen, C., Zhu, J., Wang, X., and Zhang, A. (2017) L0 trend filtering for piecewise smooth modeling.
}

\seealso{
  \code{\link{genDtf2d}} \code{\link{genD2d}} \code{\link{my.rollmean}}.
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
  y <- matrix(y,10,10)

  AMIAS2d(y, D_type = "fused.2d")
}
