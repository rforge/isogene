\name{dopamine}
\alias{dopamine}
\docType{ExpressionSet data}
\title{Dose-response microarray example data}
\description{ This dose-response microarray data contains 1000 genes and 6 doses 
(0 (control), 0.01, 0.04, 0.16, 0.63, 2.5mg/kg) with 4-5 arrays at each dose level.
}
\usage{data(dopamine)}
\format{
  An ExpressionSet object, the assayData has 1000 features and 26 samples, and in phenoData contains information of sample
names and dose levels.
  \describe{
   For the gene expression matrix obtained using the expres function, the column names are ( X1, X2, ..., X26). 
 These correspond to the dose levels (obtained using   pData function):  0, 0, 0.01, 0.01, 0.04, 0.04, 0.16, 0.16, 0.63,
 0.63, 2.50, 2.50, 0, 0, 0, 0.01, 0.01, 0.01, 0.04, 0.04, 0.16, 0.16, 0.63, 0.63, 2.50, 2.50.
  }
}



\references{
Testing for Trend in Dose-Response Microarray Experiments: a Comparison
of Testing Procedures, Multiplicity, and Resampling-Based Inference, Lin
et al. 2007, Stat. App. in Gen. & Mol. Bio., 6(1), article 26.

Gene Expression Studies Using Affymetrix Microarrays, Göhlmann, H. and Talloen, W., Chapman &
Hall/CRC, 2009

}
\examples{

data(dopamine)
express <- data.frame(exprs(dopamine))
dose <- pData(dopamine)$dose
IsoPlot(dose,express[56,],type="continuous", add.curve=TRUE)

}
\keyword{datasets}
