\name{IsoGenemSAM}
\alias{IsoGenemSAM}
\title{The five SAM test statistics calculated for both the increasing and
  decreasing trends}
\description{
  The function calculates the values for the five SAM test statistics (the
  global likelihood test, Williams, Marcus, M, and the modified M) for
  the most likely direction.
}
\usage{
IsoGenemSAM(x, y, fudge.factor)
}
\arguments{
  \item{x}{indicates the dose levels}
  \item{y}{gene expression for all genes}
  \item{fudge.factor}{the fudge factor values to be used in the SAM test statistics}
}
\value{
   A list with components 
  \item{E2}{the SAM test statistic of global likelihood test for the likely direction of each gene.}
  \item{Williams}{the test statistic of Williams for the likely direction of each gene.}
  \item{Marcus }{the test statistic of Marcus for the likely direction of each gene.}
  \item{M }{the M test statistic for the likely direction of each gene.}
  \item{ModM }{the test statistic of the modified M for the likely direction of each gene.}
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
\seealso{\code{\link{isoreg}}, \code{\link{IsoGene1}}, \code{\link{Isofudge}}}
\keyword{htest}

