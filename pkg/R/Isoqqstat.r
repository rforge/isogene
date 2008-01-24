# input :
#    x : doses
#    y : gene matrix
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

Isoqqstat <- function(x, y, niter, seed){
   ## permutations
   set.seed(seed)
   xiter.index <- t(sapply(1:niter, function(i) sample(x)))  # TV: sample stuff to remove, cf. other code
   to1 <- to2 <- to3 <- to4 <- to5 <- matrix(0, nrow(y), niter)
   for (i in 1:niter){
      yyy0 <- IsoGenem(xiter.index[i,], as.matrix(y))

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
   L <- IsoGenem(x, y)

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