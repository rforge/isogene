# input :
#    x : doses
#    y : gene expression matrix
#    fudge: fudge factor value
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

Isoqqstat <- function(x, y, fudge, niter, seed){
   ## permutations
   set.seed(seed)
   xiter.index <- t(sapply(1:niter, function(i) sample(x)))  # TODO TV: sample stuff to remove, cf. other code
   to1 <- to2 <- to3 <- to4 <- to5 <- matrix(0, nrow(y), niter)
   
   if (fudge=="pooled") {fudge.factor <- Isofudge(x,y)}
   if (fudge==0) {fudge.factor <- c(rep(0,5))}

   for (i in 1:niter){
     yyy0 <- IsoGenemSAM(xiter.index[i,], as.matrix(y), fudge.factor)

     to1[, i] <- sort(yyy0[[1]])
     to2[, i] <- sort(yyy0[[2]])
     to3[, i] <- sort(yyy0[[3]])
     to4[, i] <- sort(yyy0[[4]])
     to5[, i] <- sort(yyy0[[5]])
   }
   L <- IsoGenemSAM(x, as.matrix(y), fudge.factor)

   ## E2
   d <- L[[1]]
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to1, 1, mean)
   aa1 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
   row.names(aa1) <- row.names(y)[d.sort.list]


   ## Williams
   d <- L[[2]]
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to2, 1, mean)
   aa2 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
   row.names(aa2) <- row.names(y)[d.sort.list]

   ## Marcus
   d <- L[[3]]
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to3,1,mean)
   aa3 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
   row.names(aa3) <- row.names(y)[d.sort.list]

   ## M
   d <- L[[4]]
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]
   # Calculate the expected SAM score
   perm.mean <- apply(to4, 1, mean)
   aa4 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
   row.names(aa4) <- row.names(y)[d.sort.list]

   ## MM
   d <- L[[5]]
   d.sort.list <- sort.list(d)
   d.sort <- d[d.sort.list]

   # Calculate the expected SAM score
   perm.mean <- apply(to5, 1, mean)
   aa5 <- cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
   row.names(aa5) <- row.names(y)[d.sort.list]

     
   res <- list(aa1 = aa1, to1 = to1, aa2 = aa2, to2 = to2, aa3 = aa3,
               to3 = to3, aa4 = aa4, to4 = to4, aa5 = aa5, to5 = to5)
   
   return(res)
}