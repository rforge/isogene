IsoPlot1 <- function(x, y) {
  y1 <- as.numeric(y)[order(x)]
  m.y <- tapply(y1, as.factor(sort(x)), mean)

  plot(sort(x), y1, xlab = "dose", ylab = "gene expression")
  points(unique(x), m.y, pch = "+", cex = 0.85)
  
  title(paste("Gene: ", row.names(y), sep = ""))
}
