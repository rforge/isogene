IsoPlot1 <-
function(x,y) {
    y1 <- as.numeric(y)[order(x)]
    m.y <- tapply(y1, as.factor(x), mean)
    plot(x, y1, xlab = "dose", ylab = "gene expression") # typo ## no as.numeric needed
    points(unique(x), m.y, pch = "+", cex = 0.85)
    
    title(paste("Gene: ",row.names(y), sep = ""))
}

x <- c(1,1,1,2,2,2,3,3,3,4,4,4)
y <- rnorm(12, 1,1)
m.y <- tapply(y, as.factor(x), mean)

IsoPlot1(x,y)
