IsoPlot2 <-
function(x,y) {
    y1 <- as.numeric(y)[order(x)]
    m.y <- tapply(y1,as.factor(sort(x)),mean)

    unx <- sort(unique(x)) # compute once
    
    y.is.u <- isoreg(unx, m.y)$yf
    y.is.d <- rev(isoreg(rev(unx), m.y)$yf)

    dire <- IsoGene1(x,as.numeric(y))[[11]]

    plot(sort(x), y1, xlab = "dose", ylab = "gene expression") # typo # y1 already coerced to numeric..
    points(sort(unique(x)), m.y, pch = "+", cex = 0.85)

    if (dire == "u"){ # if else..
      points(unx, y.is.u, pch = "*", cex = 0.85)
      lines(unx, y.is.u, lty = 1,col = 4)
    } else {
      points(unx, y.is.d, pch = "*", cex = 0.85)
      lines(unx, y.is.d, lty = 1,col = 4)
    }
     title(paste("Gene: ", row.names(y), sep = ""))
}
