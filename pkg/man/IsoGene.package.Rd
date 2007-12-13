\name{IsoGene-package}
\alias{IsoGene-package}
\alias{IsoGene}
\docType{package}
\title{
IsoGene 
}
\description{
Library IsoGene aims to identify for genes with a monotonic trend in the
expression levels with respect to the increasing doses using several
test statistics. They include the global likelihood ratio test
(Bartholowmew 1961, Barlow et al. 1972 and Robertson et al. 1988),
Williams (1971, 1972), Marcus (1976), M (Hu et al. 2005) and the
modified M (Lin et al. 2007). The p-values of the five test statistics
are obtained using permutation and they are adjusted using BH (Benjamini
and Hochberg 1995) and BY (Benjamini and Yekutilie 2004) procedures are
used for controlling FDR. 
}
\details{
\tabular{ll}{
Package: \tab IsoGene\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-05-02\cr
License: \tab Free\cr
}

The package includes the following functions:

 \item{\code{\link{IsoGene1}}}{calculates the five test statistics in
   testing both increasing and decreasing alternatives for a single gene}
\item{\code{\link{IsoGenem}}}{calculates the five test statistics in
   testing both increasing and decreasing alternatives for all the genes in the data set}   
 \item{\code{\link{IsoRawp}}}{obtains the raw (one-sided and two-sided)
   p-values using permuations}
 \item{\code{\link{IsoTestBH}}}{BH or BY procedure to adjust p-values
   while controlling FDR}
 \item{\code{\link{IsoBHPlot}}}{plot of adjusted BH and BY p-values}
 \item{\code{\link{IsoPlot1}}}{plot of data and sample means}
 \item{\code{\link{IsoPlot2}}}{plot of data and isotonic means in an
   increasing or decreasing trend}
\item{\code{\link{IsopvaluePlot}}}{plot of p-values obtained using
  permutation under increasing or decreasing alternatives}
}
\author{
Lin et al.

Maintainer: Dan Lin <dan.lin@uhasselt.be>
}
\references{
Testing for Trend in Dose-Response Microarray Experiments: a Comparison
of Testing Procedures, Multiplicity, and Resampling-Based Inference, Lin
et al. 2007, submitted to Stat. App. in Gen. & Mol. Bio.
}
\keyword{package}
\seealso{
  \code{\link[<multtest>:<multtest>-package]{<multtest>}},
  \code{\link{IsoGene1}},   \code{\link{IsoGenem}},\code{\link{IsoRawp}}, \code{\link{IsoTestBH}},
  \code{\link{IsoBHPlot}}, \code{\link{IsoPlot1}}, \code{\link{IsoPlot2}}.
}
