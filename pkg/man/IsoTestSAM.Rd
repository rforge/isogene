\name{IsoTestSAM}
\alias{IsoTestSAM}
\title{Obtaining the list of significant genes using the SAM procedure}
\description{
  The function obtains the list of significant genes using the SAM
  procedure for the five test statistics (the global likelihood
  test, Williams, Marcus, M, and the modified M).
}
\usage{

IsoTestSAM(x, y, fudge, niter, seed, FDR, stat)
}
\arguments{
  \item{x}{is the dose levels }
  \item{y}{is the data frame of the gene expression with Probe ID as row names }
  \item{fudge}{is the option used for calculating the fudge factor in the SAM test statistic, either ="pooled" the fudge factor will be obtained in the program, or =0 then no fudge factor is used}
  \item{niter}{is the number of permutations to use }
  \item{seed}{random seed to generate the permutations }
  \item{FDR}{choose the desired FDR to control}
  \item{stat}{choose one of the five test statistics to use}
}

\value{
 sign.genes1: a list of genes declared significant using
 the SAM procedure in a matrix of 3 columns. The first colomn
 is the probe id, the second column is the corresponding row number
 of the probe in the dataset, and the third column is the ordered test statistic values,
and the last column is the q-values of the SAM procedure.
}


\author{Lin et al.}
\note{
  This function obtains the list of significant genes using the SAM procedure for the five test statistics. To use the SAM procedure, the number of genes in the dataset is preferably larger than 500.
}
\seealso{\code{\link{isoreg}}, \code{\link{Isofudge}}, \code{\link{IsoGenemSAM}}, \code{\link{Isoqqstat}},
\code{\link{Isoallfdr}},\code{\link{Isoqval}},
\code{\link{IsoSAMPlot}}}

\examples{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(90, 1,1),10,9) ## 10 genes with no trends
  y2 <- matrix(c(rnorm(30, 1,1),rnorm(30,2,1),rnorm(30,3,1)),10,9) ## 10 genes with increasing trends
  y <- data.frame(rbind(y1, y2)) ##y needs to be a data frame
  SAM.obj <- IsoTestSAM(x, y, fudge=0, niter=100, seed=123, FDR=0.05, stat="E2") ##no fudge factor is used
}
\keyword{htest}