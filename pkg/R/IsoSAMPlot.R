###########################################################
############### Iso SAM Plot ##############################
###########################################################

IsoSAMPlot <- function(qqstat, allfdr, FDR, stat){
   par(mfrow = c(2,2))
   FDR50 <- allfdr[,5]
   FDR90 <- allfdr[,6]
   delta <- allfdr[,1]
   fp50 <- allfdr[,2]
   fp90 <- allfdr[,3]
   signnum <- allfdr[,4]

   # FDR vs. Delta
   plot(delta, FDR50, pch = ".", ylab = "FDR")
   lines(delta, FDR90, lty = 1)
   lines(delta, FDR50, lty = 2)
   abline(0.05, 0)
   abline(0.1, 0)
   legend(0.25, 0.5, c("FDR90%","FDR50%"), lty = 1:2)
   title("a: plot of FDR vs. Delta")

   # sign vs. Delta
   plot(delta, signnum, pch = ".", ylab = "# of significant genes")
   lines(delta,signnum)
   title("b: plot of # of significant genes vs. Delta")

   # Fp vs. Delta
   plot(delta, fp50, pch = ".", ylab = "# of false positives")
   lines(delta, fp90, lty = 1)
   lines(delta, fp50, lty = 2)
   legend(0.25, nrow(qqstat[[1]])/2, c("FP90%", "FP50%"), lty = 1:2)
   title("c: plot of # of false positives vs. Delta")

   # obs vs. exp
   if (stat == "E2") {
      observed <- qqstat[[1]][,1]
      expected <- qqstat[[1]][,2]
   }
   if (stat == "Williams") {
      observed=qqstat[[3]][,1]
      expected=qqstat[[3]][,2]
   }
   if (stat == "Marcus") {
      observed=qqstat[[5]][,1]
      expected=qqstat[[5]][,2]
   }
   if (stat == "M") {
      observed=qqstat[[7]][,1]
      expected=qqstat[[7]][,2]
   }
   if (stat == "ModifM") {
      observed <- qqstat[[9]][,1]
      expected <- qqstat[[9]][,2]
   }
   del.table <- data.frame(allfdr)
   delta1 <- min(del.table[del.table[,5] <= FDR, 1])
   plot(expected, observed)
   abline(0, 1, col = "blue")
   abline(delta1, 1, lty = 5, col = "red")
   q.mat <- qqstat[[1]][order(expected),]
   x.exp <- min(q.mat[q.mat[,3] >= delta1,2])
   y.obs <- q.mat[q.mat[,2] == x.exp, 1]
   points(q.mat[q.mat[,2] >= x.exp,2], q.mat[q.mat[,1] >= y.obs,1], col = "red" )
   legend(0.45, 0.2, paste("delta=", delta1, sep=""), lty = 5)
   title("d: plot of expected vs. observed statistics")
}
