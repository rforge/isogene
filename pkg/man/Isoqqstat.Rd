\name{Isoqqstat}
\alias{Isoqqstat}
\title{Implementation of five SAM test statistics in the SAM procedure}
\description{
  The function calculates SAM test statistics from permutations for the
	five test statistics (the global likelihood test, Williams, Marcus, M,
	and the modified M).
}
\usage{
  Isoqqstat(x, y, fudge, niter)
}
\arguments{
  \item{x}{indicates the dose levels}
  \item{y}{gene expression for all genes}
  \item{fudge}{the fudge factor value to be used in the SAM test statistics: either fudge="pooled" then it is calculated by the function, or fudge="none" then no fudge factor is used}
  \item{niter}{number of permutations used in the SAM procedure}
}
\value{
  A list with components 
  
\item{aa1}{the matrix of the observed test statistic values using the
	likelihood ratio test with 4 columns: the first column contains the 
	observed test statistic values sorted in ascending order, the second
	contains the mean expected test statistic values obtained from permutations,
	the third column contains the difference between the first and the second
	column, and the last column gives the ranking of the genes in 
	ascending order.}
\item{to1}{the matrix of the test statistic values from permutations
	using the likelihood ratio test: each column of the matrix corresponds
	to the sorted test statistic from each permutation in an ascending order.}
\item{aa2}{the matrix of the observed test statistic values using
	Williams' test with 4 columns: the first column is the sorted observed
	test statistic values in an ascending order, the second is the mean
	expected test statistic values obtained from permutations, the third
	column is the difference between the first and the second column, and
	the last column is the rankings of the genes in an ascending order.}
\item{to2}{the matrix of the test statistic values from permutations
	using Williams' test: each column of the matrix corresponds to the
	sorted test statistic from each permutation in an ascending order.}
\item{aa3}{the matrix of the observed test statistic values using
	Marcus' test with 4 columns: the first column is the sorted observed
	test statistic values in an ascending order, the second is the mean
	expected test statistic values obtained from permutations, the third
	column is the difference between the first and the second column, and
	the last column is the rankings of the genes in an ascending order.}
\item{to3}{the matrix of the test statistic values from permutations
	using Marcus' test: each column of the matrix corresponds to the
	sorted test statistic from each permutation in an ascending order.}

\item{aa4}{the matrix of the observed test statistic values using the M
	test with 4 columns: the first column is the sorted observed test
	statistic values in an ascending order, the second is the mean
	expected test statistic values obtained from permutations, the third
	column is the difference between the first and the second column, and
	the last column is the rankings of the genes in an ascending order.}
\item{to4}{the matrix of the test statistic values from permutations
	using the M test: each column of the matrix corresponds to the sorted
	test statistic from each permutation in an ascending order.}
\item{aa5}{the matrix of the observed test statistic values using the
	modified M test with 4 columns: the first column is the sorted
	observed test statistic values in an ascending order, the second is
	the mean expected test statistic values obtained from permutations,
	the third column is the difference between the first and the second
	column, and the last column is the rankings of the genes in an
	ascending order.}
\item{to5}{the matrix of the test statistic values from permutations
  using the modified M test: each column of the matrix corresponds to
	the sorted test statistic from each permutation in an ascending order.}
}

\references{
Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors). (2012) Modeling Dose-response Microarray
Data in Early   Drug Development Experiments Using 
R. Springer.

IsoGene: An R Package for Analyzing Dose-response Studies in Microarray Experiments, Pramana S., Lin D., Haldermans P., Shkedy Z., Verbeke T., De Bondt A., Talloen W., Goehlmann H., Bijnens L. 2010, R Journal 2/1. 

}


\author{Lin et al.}
\note{
  This function calculates the SAM test statistics to be used in the SAM
	procedure for the five test statistics. To use the SAM procedure, the
	number of genes in the data set is preferably larger than 500. 
}
\seealso{\code{\link{isoreg}}, \code{\link{Isoallfdr}}, \code{\link{IsoGenemSAM}}
	       \code{\link{Isoqval}},\code{\link{IsoTestSAM}},\code{\link{IsoSAMPlot}}}
\examples{
  set.seed(1234)
 \donttest{
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(4500, 1,1),500,9) ## 500 genes with no trends
  y2 <- matrix(c(rnorm(1500, 1,1),rnorm(1500,2,1),
    rnorm(1500,3,1)),500,9) ## 500 genes with increasing trends
  y <- data.frame(rbind(y1, y2)) ##y needs to be a data frame
  qqstat <- Isoqqstat(x, y, fudge="pooled", niter = 50)
}
}
\keyword{htest}

