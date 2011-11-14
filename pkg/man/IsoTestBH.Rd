\name{IsoTestBH}
\alias{IsoTestBH}
\title{Test of monotonic trends using the five test statistics with BH or BY adjustment}
\description{
The function adjusts for the raw p-values of the five test statistics
using BH or BY procedure.
}
\usage{
IsoTestBH(rp, FDR, type = c("BH", "BY"), stat = c("E2",
"Williams", "Marcus", "M", "ModifM"))
}
\arguments{
  \item{rp}{raw p-value matrix with each row for one gene and 6 columns,
    the first column contains the Probe.ID, the second to the sixth
    columns are raw p-values for the five test statistics}
  \item{FDR}{the desired FDR to control}
  \item{type}{choose BH or BY procedure to control FDR}
  \item{stat}{choose one of the five test statistics to use}
}
\details{
The input raw p-values to this function can be the one sided or the two
sided ones which are obtained using function raw.p. The results using
one sided p-values and FDR controlling at alpha/2 is equivalent to that
using two sided p-values and FDR controlling at alpha.
}
\value{
  
  \item{sign.genes}{A list of significant genes while controlling FDR is
    obtained, with 4 columns: the first column is the probe ID, the
    second column is the row id, the third column is the raw p-values of
    the significant genes and the last column is the adjusted p-values
    of significant genes using BH or BY procedure}

}
\references{ package{multtest} 

Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors). (2012) Modeling Dose-response Microarray
Data in Early   Drug Development Experiments Using 
R. Springer.
}

\author{ Lin et al. }
\note{
This function only allows one type of FDR adjustment, either BH or
BY. For other type of adjustment, see function \code{mt.rawp2adjp}
in package \code{multtest}.  
}
\seealso{ 'mt.rawp2adjp', \code{\link{IsoRawp}}}
\examples{
  set.seed(1234)
  rp <- data.frame(paste("g",1:100), matrix(runif(500,0,0.1), 100, 5))
  sign <- IsoTestBH(rp, FDR = 0.05, type = "BH", stat = "E2")
}
\keyword{htest}
