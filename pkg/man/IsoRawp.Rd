\name{IsoRawp}
\alias{IsoRawp}
\title{IsoRawp}
\description{
The function calculates the raw one-sided and two-sided p-values for each
test statistic using permutations.
}
\usage{
IsoRawp(x, y, niter, progressBar=TRUE)
}
\arguments{
  \item{x}{numeric vector containing the dose levels}
  \item{y}{a data frame of the gene expression with Probe IDs as row names}
  \item{niter}{number of permutations to use}
  \item{progressBar}{to enable or disable progress bar; default is TRUE and setting it as FALSE, when problems with tcltk package 	dependecy occur, solves the problem}
}
\details{
The number of permutations to use can be chosen based on the number of
possible permutations of samples. If the possible number is too big,
usually >5000 permutations can be sufficient. 
}
\value{
A list of components
  \item{raw.p.one}{returns the one-sided p-value matrix for the five
    test statisticsin in 6 columns: the first column is the probe ID, the
    second to the last columns contain the raw p-values for each test
    statistic}
  \item{raw.p.two}{returns the two-sided p-value matrix for the five
    test statistics in 6 columns: the first column is the probe ID, the
    second to the last columns contain the raw p-values for each test
    statistic}
  \item{rawp.up}{returns the one-sided p-value matrix testing increasing
    alternative for the five test statistics in 6 columns: the first
    column is the probe ID, the second to the last columns contain the
    raw p-values for each test statistic}
 \item{rawp.dn}{returns the one-sided p-value matrix testing decreasing
   alternative for the five test statistics in 6 columns: the first
   column is the probe ID, the second to the last columns contain the
   raw p-values for each test statistic}
}

\references{ 
Modeling Dose-response Microarray
Data in Early Drug Development Experiments Using 
R, Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors), (2012), Springer.

Testing for Trend in Dose-Response Microarray Experiments: a Comparison
of Testing Procedures, Multiplicity, and Resampling-Based Inference, Lin
et al. 2007, Stat. App. in Gen. & Mol. Bio., 6(1), article 26.
}

\author{Lin et al.}
\note{ 
For each gene, the one-sided p-values are calculated from 
\eqn{\min(p^{Up}, p^{Down}}{min(p^Up, p^Down)} and the 
two sided p-values are calculated from 
\eqn{\min\{2 * \min(p^{Up},p^{Down}), 1\}}{min{2 * min(p^Up, p^Down), 1}}, 
where \eqn{p^{Up}}{p^Up} and \eqn{p^{Down}}{p^Down} are the
p-values calculated for each ordered alternative.
}
\seealso{\code{\link{IsoTestBH}}}
\examples{
\dontrun{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(90, 1,1),10,9)  # 10 genes with no trends
  y2 <- matrix(c(rnorm(30, 1,1), rnorm(30,2,1),
                 rnorm(30,3,1)), 10, 9)  # 10 genes with increasing trends
  y <- data.frame(rbind(y1, y2))  # y needs to be a data frame 
  rp <- IsoRawp(x, y, niter = 1000)
  rp
}
}
\keyword{htest}
