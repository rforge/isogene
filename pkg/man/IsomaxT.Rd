\name{IsomaxT}
\alias{IsomaxT}
\title{The maxT procedure for order restricted inference}
\description{
  The function calculates the adjusted p-values for the five test statistics (the
  global likelihood test, Williams, Marcus, M, and the modified M) using the maxT procedure.
}
\usage{
IsomaxT(x, y, niter) 
}
\arguments{
  \item{x}{indicates the dose levels}
  \item{y}{a data frame of the gene expression}
  \item{niter}{number of permutations to use}
}
\value{
   A matrix with adjusted p-values for the five test statistics.
}

\references{
Resampling based multiple testing, Westfall, P.H. and Young, S.S. 1993, Willy.

Modeling Dose-response Microarray
Data in Early Drug Development Experiments Using 
R, Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors), (2012), Springer.
}


\author{Lin et al. }

\note{
This function calculates the five test statistics using the maxT procedure 
that is controlling the Family Wise Error Rate. 
}
\seealso{\code{\link{IsoTestBH}}}

\examples{
 \donttest{
  x.res <- c(rep(1,3),rep(2,3),rep(3,3))
  y1 <- matrix(rnorm(90, 1,1),10,9)  # 10 genes with no trends
  y2 <- matrix(c(rnorm(30, 1,1), rnorm(30,2,1),
                 rnorm(30,3,1)), 10, 9)  # 10 genes with increasing trends
  dat.mat <- data.frame(rbind(y1, y2))  # y needs to be a data frame  
  niter=1000 
  
  set.seed(1234)
  pval.maxT <- IsomaxT(x.res, dat.mat,niter)
}
}
\keyword{htest}

