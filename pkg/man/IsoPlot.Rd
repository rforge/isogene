\name{IsoPlot}
\alias{IsoPlot}
\title{IsoPlot}
\description{
Plot of the data points and the sample means at each dose
}
\usage{
IsoPlot1(x, y, type=c("Continuous","Ordinal"), add.curve=TRUE/FALSE)
}
\arguments{
  \item{x}{ indicates the dose levels }
  \item{y}{ is the gene expression for one gene }
  \item{type}{ specifies the dose levels to continuous or ordinal. The default is continous.}
  \item{add.curve}{ specifies whether a fitted isotonic regression curve with a likely direction is added or not. The default is FALSE.}
}
\value{
Plot of the data points, the sample means for each dose (either as continuous or ordinal), and a fitted isotonic regression curve (optional) is produced.
}
\author{ Lin et al.}
\note{This function produces a plot for a single gene.
}
\seealso{\code{\link{IsoPlot2}}}
\examples{
  x <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
  y <- c(rnorm(3,1,1), rnorm(3,2,1), rnorm(3,3,1), rnorm(3,4,1))
  IsoPlot1(x,y)
  IsoPlot1(x,y, type="ordinal", add.curve=TRUE)
}
\keyword{hplot}
