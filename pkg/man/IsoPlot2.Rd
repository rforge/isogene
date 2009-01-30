\name{IsoPlot2}
\alias{IsoPlot2}
\title{IsoPlot2}
\description{
Plot of the data points and the isotonic means in a (likely) increasing
or decreasing trend
}
\usage{
IsoPlot2(x, y, mode=c("Continuous","Ordinal"))
}
\arguments{
  \item{x}{indicates the dose levels }
  \item{y}{is the gene expression for one gene}
  \item{mode}{ specifies the dose levels to continuous or ordinal. The default is continous.}
}
\value{
Plot of the data points and the isotonic means is produced.
}
\author{Lin et al.}
\note{This function produces a plot for a single gene.
}
\seealso{ \code{\link{IsoPlot1}} }
\examples{
  x <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
  y <- c(rnorm(3,1,1), rnorm(3,2,1), rnorm(3,3,1), rnorm(3,4,1))
  IsoPlot2(x, y)
}
\keyword{hplot}
