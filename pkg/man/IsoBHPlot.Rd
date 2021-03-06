\name{IsoBHPlot}
\alias{IsoBHPlot}
\title{Plot of adjusted p-values using BH or BY adjustment}
\description{
The function produces a plot with adjusted p-values using BH (Benjamini
and Hochberg 1995) and BY (Benjamini and Yekutieli 2004) procedures
controlling for FDR. The raw p-values and adjusted BH and BY p-values
are  plotted. 
}

\usage{
IsoBHPlot(rp, FDR, stat = c("E2", "Williams", "Marcus",
"M", "ModifM"))
}
\arguments{
  \item{rp}{raw p-value matrix with each row for one gene and 6 columns,
    the first column contains the Probe.ID, the second to the sixth
    columns are raw p-values for the five test statistics}
  \item{FDR}{the desired FDR to control}
  \item{stat}{choose one of the five test statistic to use}
}
\value{
A plot of adjusted p-values using BH and BY procedures will be produced.
}
\references{
Modeling Dose-response Microarray
Data in Early Drug Development Experiments Using 
R, Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors), (2012), Springer.

Testing for Trend in Dose-Response Microarray Experiments: a Comparison
of Testing Procedures, Multiplicity, and Resampling-Based Inference, Lin
et al. 2007, Stat. App. in Gen. & Mol. Bio., 6(1), article 26.

IsoGene: An R Package for Analyzing Dose-response Studies in Microarray Experiments, Pramana S., Lin D., Haldermans P., Shkedy Z., Verbeke T., De Bondt A., Talloen W., Goehlmann H., Bijnens L. 2010, R Journal 2/1. 

    }
\author{Lin et al.}
\seealso{ \code{\link{IsoTestBH}},  \code{\link{IsoRawp}} }
\examples{
  rp <- data.frame(paste("g", 1:100), matrix(runif(500,0,1), 100, 5))
  IsoBHPlot(rp, FDR = 0.05, stat = "E2")
}
\keyword{hplot}
