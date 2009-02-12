\name{IsopvaluePlot}
\alias{IsopvaluePlot}
\title{Plot of p-values from permutations under increasing or decreasing alternatives }
\description{
The function calculates the p-values using permutations under increasing
and decreasing ordered alternatives for one gene. The p-values (\eqn{p^{up}}{p^up} and
\eqn{p^{down}}{p^down}) are obtained from the plot of null distribution and observed statistics.
}
\usage{
IsopvaluePlot(x, y, niter, stat = c("E2", "Williams", "Marcus", "M", "ModifM"))
}
\arguments{
  \item{x}{the dose levels }
  \item{y}{the gene expressions}
  \item{niter}{the number of permutations to use}
  \item{stat}{choose one of the five test statistics to use}
}
\value{
 Plots of the null distribution and the observed test statistic under
 increasing and decreasing ordered alternatives.
}
\author{ Lin et al. }
\note{The function obtains the p-values under increasing and decreasing
  ordered alternatives for a single gene.
}
\seealso{\code{\link{IsoGene1}}}
\examples{
x <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
y <- c(rnorm(3,1,1), rnorm(3,2,1), rnorm(3,3,1), rnorm(3,4,1)) 

IsopvaluePlot(x, y, niter = 1000, stat = "Williams")
}
\keyword{hplot}
