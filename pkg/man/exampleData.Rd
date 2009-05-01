\name{exampleData}
\alias{exampleData}
\docType{data}
\title{Dose-response microarray example data}
\description{
This dose-response microarray data contains 1000 genes and 4 doses 
(one control dose (zero dose) and three increasing dose) with 3 arrays at each dose level.
}
\usage{data(exampleData)}
\format{
  A data frame with 1000 observations on the following 12 variables.
  \describe{
    \item{\code{X1}}{Sample one with zero dose}
    \item{\code{X1.1}}{Sample two with zero dose}
    \item{\code{X1.2}}{Sample three with zero dose}
    \item{\code{X2}}{Sample one with second dose}
    \item{\code{X2.1}}{Sample two with second dose}
    \item{\code{X2.2}}{Sample three with second dose}
    \item{\code{X3}}{Sample one with third dose}
    \item{\code{X3.1}}{Sample two with third dose}
    \item{\code{X3.2}}{Sample three with third dose}
    \item{\code{X4}}{Sample one with fourth dose}
    \item{\code{X4.1}}{Sample two with fourth dose}
    \item{\code{X4.2}}{Sample three with fourth dose}
  }
}



\references{
Testing for Trend in Dose-Response Microarray Experiments: a Comparison
of Testing Procedures, Multiplicity, and Resampling-Based Inference, Lin
et al. 2007, Stat. App. in Gen. & Mol. Bio., 6(1), article 26.
}
\examples{
data(exampleData)
x <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3))
gene1 <- as.numeric(exampleData[1,])
IsoPlot(x, gene1)
}
\keyword{datasets}
