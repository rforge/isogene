\name{IsoTestSAM}
\alias{IsoTestSAM}
\title{Obtaining the list of significant genes using the SAM procedure}
\description{
  The function obtains the list of significant genes using the SAM procedure for the five test statistics (the
  global likelihood test, Williams, Marcus, M, and the modified M).
}
\usage{

IsoTestSAM(x, y, niter, seed, FDR, stat)
}
\arguments{
  \item{x}{is the dose levels }
  \item{y}{is the data frame of the gene expression with Probe ID as row names }
  \item{niter}{is the number of permutations to use }
  \item{seed}{random seed to generate the permutations }
  \item{FDR}{choose the desired FDR to control}
  \item{stat}{choose one of the five test statistics to use}
}

\value{
 sign.genes1: a list of genes declared significant using
 the SAM procedure in a matrix of 3 columns. The first colomn
 is the probe id, the second column is the corresponding row number
 of the probe in the dataset, and the third column is the q-values of
 the SAM procedure.
}

\references{
Testing for Trend in Dose-Response Microarray Experiments: a Comparison
of Testing Procedures, Multiplicity, and Resampling-Based Inference, Lin
et al. 2007, Stat. App. in Gen. & Mol. Bio., 6(1), article 26.
  }
\author{Lin et al.}
\note{
  This function obtains the list of significant genes using the SAM procedure for the five test statistics. To use the SAM procedure, the number of genes in the dataset is preferably larger than 500.
}
\seealso{\code{\link{isoreg}}, \code{\link{Isoqqstat}},
\code{\link{Isoallfdr}},\code{\link{Isoqval}},
\code{\link{IsoSAMPlot}}}

\examples{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(9000, 1,1),1000,9) ## 1000 genes with no trends
  y2 <- matrix(c(rnorm(3000, 1,1),rnorm(3000,2,1),rnorm(3000,3,1)),1000,9) ## 1000 genes with increasing trends
  y <- data.frame(rbind(y1, y2)) ##y needs to be a data frame
  SAM.obj <- IsoTestSAM(x, y, niter=100, seed=123, FDR=0.09, stat="E2")

  }
\keyword{htest}