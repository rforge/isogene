###########################################################
############### Iso SAM Function ##########################
###########################################################

# output:
# list of significant genes
#    col 1: Probe.ID
#    col 2: row number
#    col 3: q value

IsoTestSAM <- function(x, y, fudge, niter, seed, FDR, stat) {
   qqstat <- Isoqqstat(x, y, fudge, niter, seed)
   allfdr <- Isoallfdr(qqstat, , stat)
   del.table <- data.frame(allfdr)
   min_fdr <- min(na.exclude(del.table[, 5]))
   if (min_fdr > FDR) {
      FDR <- min_fdr
      delta <- min(na.exclude(del.table[del.table[,5] <= FDR,1]))
      print("FDR cannot be obtained in this dataset")
   } else {
      delta <- min(na.exclude(del.table[del.table[,5] <= FDR,1]))
   }
   qval <- Isoqval(delta,allfdr,qqstat,stat)
   q.value <- qval[[1]]
   sign.list <- q.value[q.value[,3] <= FDR,]
   sign.genes <- cbind(row.names(y[sign.list[,1],]), sign.list)
   sign.genes1 <- data.frame(sign.genes[order(sign.list[,2]),])
   row.names(sign.genes1) <- 1:nrow(sign.genes1)
   names(sign.genes1) <- c("Probe.ID", "row.number","stat.val","qvalue")

   return(sign.genes1)
}