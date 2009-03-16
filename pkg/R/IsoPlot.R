
IsoPlot <- function (x, y, type = c("continuous", "ordinal"),add.curve=c("FALSE","TRUE")) 

{
            options(warn = -1)
            
            if (is.na(match(type, c("ordinal", "continuous")))) {
                        options(warn = 0)
                        print("Error: The dose can be only continuous or ordinal")
            }
            else {
                        miny <- min(y)
                        maxy <- max(y)
                        ordx <- order(x)
                        unx <- sort(unique(x))
                        y1 <- as.numeric(y)[ordx]
                        y.m <- tapply(y1, as.factor(x[order(x)]), mean)
                        y.m.tot <- rep(mean(y), length(unx))
                        n.p <- table(x)
                        n.g <- length(n.p)
                        y.is.u <- pava(y.m, wt = n.p)
                        y.is.d <- rev(pava(rev(y.m), wt = rev(n.p)))
                        dire <- IsoGene1(x, as.numeric(y))[[11]]
                        if (missing(type) | type == "continuous") {
                                    plot(sort(x), y1, lwd = 2, xlab = "Doses", ylab = "Gene Expression")
                                    points(sort(unique(x)), y.m, pch = "+", lwd = 3,col=2)
                                    if (dire == "u") {
                                                points(unx, y.is.u, pch = "*", lwd = 2)
                                                if (add.curve=="TRUE") 
                                                            lines(unx, y.is.u, lty = 1, col = 4, lwd = 2)
                                    }
                                    else {
                                                points(unx, y.is.d, pch = "*", lwd = 2)
                                                if (add.curve=="TRUE")
                                                            lines(unx, y.is.d, lty = 1, col = 4, lwd = 2) 
                                    }
                        }
                        if (type == "ordinal") {
                                    catx <- factor(x, levels = sort(unique.default(x)), 
                                                            labels = unx, ordered = F)
                                    a <- c(1:length(unx))
                                    plot(a, ylim = c(miny, maxy), pch = "", ylab = "Gene Expression", 
                                                            xlab = "Doses", axes = FALSE)
                                    axis(1, sort(unique(a)), as.character(unx))
                                    axis(2)
                                    points(catx, y, lwd = 2)
                                    points(sort(unique(catx)), y.m, pch = "+", lwd = 3,col=2)
                                    if (dire == "u") {
                                                points(a, y.is.u, pch = "*", lwd = 2)
                                                if (add.curve=="TRUE") {lines(a, y.is.u, lty = 1, col = 4, lwd = 2)}
                                    }
                                    else {
                                                points(a, y.is.d, pch = "*", lwd = 2)
                                                if (add.curve=="TRUE") { lines(a, y.is.d, lty = 1, col = 4, lwd = 2) }
                                    }
                        }
                        title(paste("Gene: ", row.names(y), sep = ""))
            }
}

