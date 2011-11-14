\name{IsoPlot}
\alias{IsoPlot}
\title{IsoPlot}
\description{
Plot of the data points and the sample means at each dose
}
\usage{
IsoPlot(x, y, type=c("continuous", "ordinal"), add.curve = FALSE)
}
\arguments{
  \item{x}{ indicates the dose levels }
  \item{y}{ is the gene expression for one gene }
  \item{type}{ specifies the dose levels to \code{"continuous"} or \code{"ordinal"}. The default is \code{"continous"}.}
  \item{add.curve}{ specifies whether a fitted isotonic regression curve with a likely direction is added or not. The default is \code{FALSE}.}
}
\value{
Plot of the data points, the sample means for each dose (either as continuous or ordinal), and a fitted isotonic regression curve (optional) is produced.
}

\references{
Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors). (2012) Modeling Dose-response Microarray
Data in Early   Drug Development Experiments Using 
R. Springer.
}


\author{ Lin et al.}
\note{This function produces a plot for a single gene.
}
\examples{
  x <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
  y <- c(rnorm(3,1,1), rnorm(3,2,1), rnorm(3,3,1), rnorm(3,4,1))
  IsoPlot(x, y)
  IsoPlot(x, y, type="ordinal", add.curve=TRUE)
}
\keyword{hplot}
