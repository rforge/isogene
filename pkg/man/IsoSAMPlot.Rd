\name{IsoSAMPlot}
\alias{IsoSAMPlot}
\title{Plots produced using the SAM procedure}
\description{
  The function produces four plots using the SAM procedure for
  one of the five test statistics (the likelihood ratio test, Williams,
  Marcus, the M and modified M tests): FDR vs. delta, number
  of significant genes vs. delta, number of false positives vs. delta, and the observed
  vs. expected SAM test statistics obtained from permutations.
}
\usage{

IsoSAMPlot(qqstat, allfdr, FDR, stat)
}
\arguments{
  \item{qqstat}{output from function Isoqqstat containing the test statistics of permutations}
  \item{allfdr}{the delta table obtained from function Isoallfdr}
  \item{FDR}{choose the desired FDR to control}
  \item{stat}{choose one of the five test statistics to use}
}

\value{
returns four plots produced using the SAM procedure.

}

\author{Lin et al.}
\note{
  This function produces four plots using the SAM procedure for the five test statistics. To use the SAM procedure, the number of genes in the dataset is preferably larger than 500.
}
\seealso{\code{\link{isoreg}}, \code{\link{Isoqqstat}},
\code{\link{Isoallfdr}},\code{\link{Isoqval}},
\code{\link{IsoTestSAM}}}

\examples{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(4500, 1,1),500,9) ## 500 genes with no trends
  y2 <- matrix(c(rnorm(1500, 1,1),rnorm(1500,2,1),rnorm(1500,3,1)),500,9) ## 500 genes with increasing trends
  y <- data.frame(rbind(y1, y2)) ##y needs to be a data frame
  qqstat <- Isoqqstat(x, y, fudge="pooled", niter=50)
  allfdr <- Isoallfdr(qqstat, , stat = "E2")
  IsoSAMPlot(qqstat, allfdr, FDR = 0.1, stat = "E2")
}
\keyword{hplot}
