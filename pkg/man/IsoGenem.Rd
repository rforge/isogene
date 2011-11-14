\name{IsoGenem}
\alias{IsoGenem}
\title{The five test statistics calculated for both the increasing and
  decreasing trends}
\description{
  The function calculates the values for the five test statistics (the
  global likelihood test, Williams, Marcus, M, and the modified M) for
  testing increasing and decreasing alternatives.
}
\usage{
IsoGenem(x, y)
}
\arguments{
  \item{x}{indicates the dose levels}
  \item{y}{gene expression for all genes}
}
\value{
   A list with components 
  \item{E2.up}{the test statistic of global likelihood test for testing
    increasing alternative.}
  \item{Williams.up}{the test statistic of Williams for testing
    increasing alternative.}
  \item{Marcus.up }{the test statistic of Marcus for testing increasing
    alternative.}
  \item{M.up }{the M test statistic for testing increasing alternative.}
\item{ModM.up }{the test statistic of the modified M for testing
  increasing alternative.}
  \item{E2.dn }{the test statistic of Williams for testing increasing alternative.}
\item{Williams.dn }{the test statistic of global likelihood test for
  testing increasing alternative.}
  \item{Marcus.dn }{the test statistic of Williams for testing
    increasing alternative.}
\item{M.dn}{the test statistic of global likelihood test for testing
  increasing alternative.}
  \item{ModM.dn}{the test statistic of Williams for testing increasing alternative.}
\item{direction}{the direction with the higher likelihood of increasing
  (indicated by "u") or decreasing (indicated by "d") trend.}

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
\author{Lin et al. }
\note{
  This function calculates the five test statistics for both increasing
  and decreasing ordered alternatives for all the genes (rows in the data set). 
}
\seealso{\code{\link{isoreg}}, \code{\link{IsoGene1}}}
\examples{
\dontrun{
  set.seed(1234)
  x <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(90, 1,1),10,9)  # 10 genes with no trends
  y2 <- matrix(c(rnorm(30, 1,1), rnorm(30,2,1),
                 rnorm(30,3,1)), 10, 9)  # 10 genes with increasing trends
  y <- data.frame(rbind(y1, y2))  # y needs to be a data frame 
  stat <- IsoGenem(x,y) 
  stat
}
}

\keyword{htest}

