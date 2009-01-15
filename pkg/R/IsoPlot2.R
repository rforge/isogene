IsoPlot2 <- function(x, y) {
  ordx <- order(x)
  x <- x[ordx]
  unx <- sort(unique(x))
  
  y1 <- as.numeric(y)[x]
  y.m <- tapply(y1, as.factor(x), mean)
  y.m.tot <- rep(mean(y), length(unx))
  
  
  n.p <- table(x)
  n.g <- length(n.p) 
 # y.is.u <- isoreg(unx, m.y)$yf
#  y.is.d <- rev(isoreg(rev(unx), m.y)$yf)

   y.is.u <- pava(y.m, wt=n.p )
   y.is.d <- rev(pava(rev(y.m), wt=rev(n.p)))

  dire <- IsoGene1(x,as.numeric(y))[[11]]

  plot(sort(x), y1, xlab = "dose", ylab = "gene expression")
  points(sort(unique(x)), y.m, pch = "+", cex = 0.85)

  if (dire == "u"){ 
    points(unx, y.is.u, pch = "*", cex = 0.85)
    lines(unx, y.is.u, lty = 1,col = 4)
  } else {
    points(unx, y.is.d, pch = "*", cex = 0.85)
    lines(unx, y.is.d, lty = 1,col = 4)
  }
   title(paste("Gene: ", row.names(y), sep = ""))
}
