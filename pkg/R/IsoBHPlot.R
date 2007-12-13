"IsoBHPlot" <-
function(rp, FDR, stat = c("E2", "Williams", "Marcus", "M", "ModifM")){

    stat <- match.arg(stat)
    
    Probe.ID <- rp[,1]
    if (stat == "E2") rpraw <- rp[,2]
    else if (stat == "Williams") rpraw <- rp[,3]
    else if (stat == "Marcus") rpraw <- rp[,4]
    else if (stat == "M") rpraw <- rp[,5]
    else rpraw <- rp[,6] # (stat == "ModifM") 

    # library(multtest)
    procs <- c("Bonferroni", "Holm", "BH", "BY")
    res <- mt.rawp2adjp(rpraw, procs)
    adjp <- res$adjp[order(res$index), ]

    plot(1:nrow(rp), sort(adjp[,1]), #no c()
         col = 4, pch = ".", lty = 1, xlab = "index",
         ylab = "Adjusted Aymptotic P values")

    lines(1:nrow(rp), sort(adjp[,1]), lty = 1, col = 1) # no c()
    lines(1:nrow(rp), sort(adjp[,4]), lty = 4, col = 2)
    lines(1:nrow(rp), sort(adjp[,5]), lty = 5, col = 3)
    abline(FDR, 0, lty = 6)

    legend(nrow(rp) / 2,
           0.3, col = c(1,2,3), c("Raw P","BH(FDR)","BY(FDR)"), lty = c(1,4:5))
    title(paste(stat,": Adjusted p values by BH and BY", sep = ""))
}

