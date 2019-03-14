\name{my.rollmean}
\alias{my.rollmean}
\title{Rolling Means}
\description{
  Generic functions for computing rolling means of ordered observations.
}
\usage{
my.rollmean(y, h = 5 ,...)
}
\arguments{
  \item{y}{Observe sequence, with length n.}
  \item{h}{integer width of the rolling window, \code{\link[zoo]{rollmean}} would be used with \code{k = h*2+1}.}
  \item{\dots}{Other arguments, see \code{\link[zoo]{rollmean}} for details.}
}
\value{
  \item{vector}{Vector: smooth observe sequence, with length n.}
}

\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.
  
  Maintainer: Canhong Wen <wench@ustc.edu.cn>
}
\references{
  Wen, C., Zhu, J., Wang, X., and Zhang, A. (2019) \emph{L0 trend filtering}, technique report.
}
