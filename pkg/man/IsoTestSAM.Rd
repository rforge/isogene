\name{IsoTestSAM}
\alias{IsoTestSAM}
\title{Obtaining the list of significant genes using the SAM procedure}
\description{
  The function obtains the list of significant genes using the SAM
  procedure for the five test statistics (the global likelihood
  test, Williams, Marcus, M, and the modified M).
}
\usage{
IsoTestSAM(x, y, fudge, niter, FDR, stat)
}
\arguments{
  \item{x}{numeric vector containing the dose levels }
  \item{y}{data frame of the gene expression with Probe ID as row names }
  \item{fudge}{option used for calculating the fudge factor in the SAM test 
    statistic, either \code{"pooled"} (fudge factor will be automatically computed in the function), or \code{"none"} if no fudge factor is used}
  \item{niter}{number of permutations to use}
  \item{FDR}{choose the desired FDR to control}
  \item{stat}{choose one of the five test statistics to use}
}

\value{
  A list with components 
  \item{sign.genes1}{a list of genes declared significant using
  the SAM procedure in a matrix of 5 columns. The first colomn
  is the probe id, the second column is the corresponding row number
  of the probe in the dataset, and the third column is the ordered test statistic values,
  and the fourth column is the q-values of the SAM procedure. The last two columns
  are raw p-values based on permutations and BH adjusted p-values.}
  \item{qqstat}{output of Isoqqstat}
  \item{allfdr}{output of Isoallfdr}
}

\author{Lin et al.}
\note{
  This function obtains the list of significant genes using the SAM procedure for the five test statistics. To use the SAM procedure, the number of genes in the dataset is preferably larger than 500.
}
\seealso{\code{\link{isoreg}}, \code{\link{Isofudge}}, \code{\link{IsoGenemSAM}}, \code{\link{Isoqqstat}},
  \code{\link{Isoallfdr}},\code{\link{Isoqval}}, \code{\link{IsoSAMPlot}}}
\examples{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(4500, 1,1),500,9) ## 500 genes with no trends
  y2 <- matrix(c(rnorm(1500, 1,1),rnorm(1500,2,1),rnorm(1500,3,1)),500,9) ## 500 genes with increasing trends
  y <- data.frame(rbind(y1, y2)) ##y needs to be a data frame
  SAM.obj <- IsoTestSAM(x, y, fudge="pooled", niter=50, FDR=0.05, stat="E2") 
}
\keyword{htest}
