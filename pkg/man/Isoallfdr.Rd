\name{Isoallfdr}
\alias{Isoallfdr}
\title{Obtaining the delta table in the SAM procedure}
\description{
  The function obtains the delta table in the SAM procedure for the five test statistics (the
  global likelihood test, Williams, Marcus, M, and the modified M).
}
\usage{
Isoallfdr(qqstat, ddelta, stat)
}
\arguments{
  \item{qqstat}{output from function Isoqqstat containing the test statistics of permutations}
  \item{ddelta}{give a list of values as cut-off to find the number of significant genes in the SAM procedure. 
  If unspecified, the default value is assigned using the centiles of the absolute difference between the observed and expected test statistics.}
  \item{stat}{choose one of the five test statistics to use}
}


\value{
   dtable: the delta table in the SAM procedure containing
   six columns. The first column is the cut-off value to find
   the number of significant genes, the second column is the
   meidan number of false positives, the third column is the 90\%
   percentile number of false positives, the fourth column is
   the number of significant genes, the fifth column is the median
   FDR, and the last column is the 90\% FDR.
}

\author{Lin et al.}
\note{
  This function calculates the delta table in the SAM procedure for the five test statistics. To use the SAM procedure, the number of genes in the dataset is preferably larger than 500.
}
\seealso{\code{\link{isoreg}}, \code{\link{Isoqqstat}},
\code{\link{Isoqval}},\code{\link{IsoTestSAM}},
\code{\link{IsoSAMPlot}}}

\examples{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(4500, 1,1),500,9) ## 500 genes with no trends
  y2 <- matrix(c(rnorm(1500, 1,1),rnorm(1500,2,1),rnorm(1500,3,1)),500,9) ## 500 genes with increasing trends
  y <- data.frame(rbind(y1, y2)) ##y needs to be a data frame
  qqstat <- Isoqqstat(x, y, fudge="pooled",niter=50)
  allfdr <- Isoallfdr(qqstat,,stat="E2")
 }
\keyword{htest}

