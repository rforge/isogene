# input :
#    gid : doses
#    x : gene matrix
#    niter : number of permutation, >100 is good, =500 is sufficient
#    seed : random seed
# output:
#    length of 10
#    [[1]] observed stat of E2
#    [[2]] permutation stat matrix of E2 
#    [[3]] observed stat of Williams
#    [[4]] permutation stat matrix of Williams 
#    [[5]] observed stat of Marcus
#    [[6]] permutation stat matrix of Marcus 
#    [[7]] observed stat of M
#    [[8]] permutation stat matrix of M 
#    [[9]] observed stat of ModifM
#    [[10]] permutation stat matrix of ModifM 

Isoqqstat <- function(gid, x, niter, seed){
   ## permutations
   set.seed(seed)
   xiter.index <- t(sapply(1:niter, function(i) sample(gid)))  # TV: sample stuff to remove, cf. other code
   to1 <- to2 <- to3 <- to4 <- to5 <- matrix(0, nrow(x), niter) 
   for (i in 1:niter){
      yyy0 <- IsoGenem(xiter.index[i,], x)
      
      yyy <- apply(cbind(yyy0[[1]], yyy0[[6]]), 1, max)
      to1[, i] <- sort(yyy)
      
      yyy <- apply(cbind(yyy0[[2]], yyy0[[7]]), 1, max)
      to2[, i] <- sort(yyy)
      
      yyy <- apply(cbind(yyy0[[3]], yyy0[[8]]), 1, max)
      to3[, i] <- sort(yyy)
      
      yyy <- apply(cbind(yyy0[[4]], yyy0[[9]]), 1, max)
      to4[, i] <- sort(yyy)
      
      yyy <- apply(cbind(yyy0[[5]], yyy0[[10]]), 1, max)
      to5[, i] <- sort(yyy)
      # print(i)
   }
   L <- IsoGenem(gid, x)

   ## E2
   d <- apply(cbind(L[[1]], L[[6]]), 1, max)
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to1, 1, mean)
   aa1 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)

   ## Williams
   d <- apply(cbind(L[[2]], L[[7]]), 1, max)
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to2, 1, mean)
   aa2 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
 
   ## Marcus
   d <- apply(cbind(L[[3]], L[[8]]), 1, max)   
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]   
   #Calculate the expected SAM score
   perm.mean <- apply(to3,1,mean)
   aa3 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
       
   ## M
   d <- apply(cbind(L[[4]], L[[9]]), 1, max)   
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to4, 1, mean)
   aa4 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
       
   ## MM
   d <- apply(cbind(L[[5]], L[[10]]), 1, max)      
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   
   # Calculate the expected SAM score
   perm.mean <- apply(to5, 1, mean)
   aa5 <- cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)

   res <- list(aa1 = aa1, to1 = to1, aa2 = aa2, to2 = to2, aa3 = aa3,
               to3 = to3, aa4 = aa4, to4 = to4, aa5 = aa5, to5 = to5)
   return(res)
}
     
# Arguments :  qqstat, the output from Isoqqstat
#              ddelta contains various values for delta's
# Return FDR for various values of delta

# input :
#    qqstat: the output from function Isoqqstat
#    ddelta: can be missing
#    stat: "E2", "Williams", "Marcus", "M","ModifM" 
# output :
#    delta table:
#    col 1: delta values
#    col 2: FP 50%
#    col 3: FP 90%
#    col 4: Number of significant genes
#    col 5: FDR 50%
#    col 6: FDR 90%

Isoallfdr <- function(qqstat, ddelta, stat) {
  
   if (stat == "E2"){
      qstat <- qqstat[[1]]
      dperm <- qqstat[[2]]
   }
   if (stat == "Williams") {
      qstat <- qqstat[[3]]
      dperm <- qqstat[[4]]
   }     
   if (stat == "Marcus") {
      qstat <- qqstat[[5]]
      dperm <- qqstat[[6]]
   }
   if (stat == "M") {
      qstat <- qqstat[[7]]
      dperm <- qqstat[[8]]
   }
   if (stat == "ModifM") {
      qstat <- qqstat[[9]]
      dperm <- qqstat[[10]]
   }
   
   k1 <- nrow(qstat)
   # Use percentiles if ddelta is not specified
   if (missing(ddelta)) {
      qd <- round(quantile(abs(qstat[,3]), c(0.01, 0.999)), 2)
      bb <- (qd[2] - qd[1]) / 100
      b0 <- 10^-floor(log(bb, 10))
      bb <- trunc(bb * b0) / b0
      ddelta <- seq(from = qd[1], to = qd[2], by = bb)
   }
   
   k2 <- length(ddelta)
   # Suppress warning messages if uupcut is empty
   # Determine the number of significant genes for each value in ddelta
   low.point <- up.point <- clow <- cup <- sn <- NULL
   nsn <- array(0, c(length(ddelta), 2))
   for (i in 1:k2){
      low <- max(qstat[qstat[,1] < 0, 1][qstat[qstat[,1] < 0, 3] < - ddelta[i]])
      if (low == "-Inf") { # use  !is.finite(low) combined with negative ??
         which.low <- which(qstat[qstat[,1] < 0, 1] == low)
         low.point[i] <- which.low[length(which.low)]
         # max(sort(as.numeric(grep("TRUE",as.character(test.low)))))
      }
      if (low == "-Inf"){
         low <- min(qstat[,1])
         low.point[i] <- 0
      }
      clow[i] <- low
      up <- min(qstat[qstat[,1] > 0, 1][qstat[qstat[,1] > 0,3] > ddelta[i] ])
      if (up != "Inf") {
        which.up <- which(qstat[qstat[,1]>0,1] == up)
        up.point[i] <- which.up[1]
      }
      if (up == "Inf") {
         up <- max(qstat[qstat[,1] > 0,1])
         up.point[i] <- sum(qstat[,1] > 0) + 1
      }
      cup[i] <- up
      sn[i] <- sum(qstat[,1] > 0) + 1 - sum(up.point[i]) + sum(low.point[i])
      qq <- c(rep(1, nrow(dperm)) %*% ((dperm > cup[i]) | (dperm < clow[i])))
      nsn[i, ] <- quantile(qq, c(0.5, 0.9))
    }
    prun <- quantile(dperm, c(0.25, 0.75))
    p <- min(sum(qstat[,1] < prun[2] & qstat[,1] > prun[1]) / (nrow(qstat) *0.5), 1)
    nsn <- nsn * p
    dtable = round(cbind(Ddelta = ddelta, "FalsePositive50%" = nsn[,1],
       "FalsePositive90%" = nsn[,2], Called = sn, "FDR50%" = nsn[,1]/sn,
       "FDR90%" = nsn[, 2]/sn), 4)
    return(dtable)
}

# Arguments:  qqstat, the output from Isoqqstat
#             delta, the cutoff (critical value) for the rejection region
# Return the 50th and 90th percentile false discovery rate for a given delta
# output :
#    length of 2
#    [[1]] : 2 columns for all the genes
#       col1: row.number
#       col2: q value
#    [[2]] : 2 columns for significant genes only  
#       col1: row.number
#       col2: q value

Isoqval <- function (delta, allfdr, qqstat, stat) {
  
   if (stat == "E2") qstat <- qqstat[[1]]
   if (stat == "Williams") qstat <- qqstat[[3]]    
   if (stat == "Marcus") qstat <- qqstat[[5]]    
   if (stat == "M") qstat <- qqstat[[7]]
   if (stat == "ModifM") qstat <- qqstat[[9]]
   
   low <- max(qstat[qstat[,1] < 0, 1][qstat[qstat[,1] < 0, 3] < - delta])
   if (low != "-Inf") {
      which.low <- which(qstat[qstat[,1]<0,1] == low)
      low.point <- which.low[length(which.low)]
   }
   if (low == "-Inf") {
      low <- min(qstat[,1])
      low.point <- 0
   }
   up <- min(qstat[qstat[,1] > 0, 1][qstat[qstat[,1] > 0, 3] > delta ])
   if (up == "Inf") {
      up <- max(qstat[qstat[,1] > 0, 1])
      up.point <- nrow(qstat) + 1
   }
   if (up != "Inf") {
      which.up <- which(qstat[qstat[,1] > 0, 1] == up)
      up.point <- sum(qstat[,1] <= 0) + which.up[1]
   }
   dtable <- allfdr
   k <- nrow(dtable) - 1
   nnew <- matrix(0, nrow(qstat), 2)
   test <- abs(qstat[,3]) < dtable[1, 1]
   nnew[test, 1] <- qstat[test, 4]
   nnew[test, 2] <- dtable[1, 5]
   test2 <- abs(qstat[,3]) >= dtable[k + 1, 1]
   nnew[test2, 1] <- qstat[test2, 4]
   nnew[test2, 2] <- dtable[k + 1, 5]
   for (i in k:1) {
      j <- i + 1 
      test1 <- dtable[i, 1] <= abs(qstat[, 3]) & abs(qstat[, 3]) < dtable[j, 1]
      nnew[test1, 1] <- qstat[test1, 4]
      nnew[test1, 2] <- dtable[i, 5]
   }
   nnnew <- cbind(nnew,qstat[,3], nnew[,2])
   fdr <- NULL
   ## match delta with FDR in the delta table
   for (i in k:1) {
      j <- i + 1
      test2 <- dtable[i,1] <= delta & delta < dtable[j,1]
      if (sum(test2) > 0){
         fdr <- dtable[i, 5]
      }
      test3 <- dtable[i, 1] < delta & delta <= dtable[j,1]
      if (sum(test3) > 0) {
         fdr <- dtable[j, 5]
      }
   }
   n1 <- n2 <- NULL
   if (low.point == 0 & is.numeric(fdr)) {
      n1 <- NULL
   }
   if (low.point == 1 & is.numeric(fdr)) {
      n1<- matrix(nnnew[1,], 1, 4)
   }
   if (low.point > 1 & is.numeric(fdr)) {
      n1 <- nnnew[1:low.point,]
      n1[n1[,3]> -delta, 4] <- fdr
   }
   if (up.point == nrow(qstat)+1 & is.numeric(fdr)) {
      n2 <- NULL
   }
   if (up.point == nrow(qstat) & is.numeric(fdr)) {
      n2 <- matrix(nnnew[nrow(qstat),], 1, 4)
   }
   if (up.point < nrow(qstat) & is.numeric(fdr)) {
      n2 <- nnnew[up.point:nrow(qstat),]
      n2[n2[, 3] < delta,4] <- fdr
   }
   res <- sign.list <- NULL
   m1 <- low.point + 1
   m2 <- up.point - 1
   cols <- c(1,4)
   if (is.numeric(fdr) & is.numeric(n1)) {
      n1[n1[,4] > fdr,4] <- fdr
   }
   if (is.numeric(fdr) & is.numeric(n2)) {
      n2[n2[,4] > fdr,4] <- fdr
   }
   if (is.numeric(fdr) & is.numeric(n1) & is.numeric(n2)) {
      res <- rbind(n1[, cols], nnnew[m1:m2, cols], n2[,cols])
      sign.list <- rbind(n1[, cols], n2[, cols])
   }
   if (is.numeric(fdr) & is.numeric(n1) & is.numeric(n2) == F) { ## TV: == F ? what is this supposed to be ?
      res <- rbind(n1[,cols], nnnew[m1:m2, cols])
      sign.list <- rbind(n1[,cols])
   } 
   if (is.numeric(fdr) & is.numeric(n1) == F & is.numeric(n2)) { ## TV: == F ?
      res <- rbind(nnnew[m1:m2, cols],n2[, cols])
      sign.list <- rbind(n2[, cols])
   } 
   return(res, sign.list)
}

###########################################################
############### Iso SAM Function ##########################
###########################################################
 
# output:
# list of significant genes
#    col 1: Probe.ID
#    col 2: row number
#    col 3: q value

IsoTestSAM <- function(x.res, dat.mat, niter, seed, FDR, stat) {
   qqstat <- Isoqqstat(x.res, dat.mat, niter, seed)
   allfdr <- f.allfdr(qqstat, , stat)
   del.table <- data.frame(allfdr)
   min_fdr <- min(del.table[, 5])
   if (min_fdr > FDR) {
      FDR <- min_fdr
      delta <- min(del.table[del.table[,5] <= FDR,1])
   } else {
      delta <- min(del.table[del.table[,5] <= FDR,1])
   }
   qval <- f.qval(delta,allfdr,qqstat,stat)
   q.value <- qval[[1]]
   sign.list <- q.value[q.value[,2] <= FDR,]
   sign.genes <- cbind(row.names(dat.mat[sign.list[,1],]), sign.list)
   sign.genes1 <- data.frame(sign.genes[order(sign.list[,2]),])
   row.names(sign.genes1) <- 1:nrow(sign.genes1)
   names(sign.genes1) <- c("Probe.ID", "row.number", "qvalue")
   
   return(sign.genes1)
}

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
