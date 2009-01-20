
IsoRawp <- function(x, y, niter, seed) {

  if (niter > 2500) {
    jmatsize <- 2500 # TODO turn into function argument
    nmats <- floor(niter / jmatsize)
    endmat <- c(1:nmats * jmatsize, niter)
    begmat <- c(1, endmat[-length(endmat)] + 1)
  } else {
    nmats <- 1
    begmat <- 1
    endmat <- niter
  }

  y <- data.matrix(y, rownames.force = TRUE)

  E <- IsoGenem(x, y)

  obs.E.up  <- matrix(E[[1]], nrow(y), 1)
  obs.W.up  <- matrix(E[[2]], nrow(y), 1)
  obs.WC.up <- matrix(E[[3]], nrow(y), 1)
  obs.M.up  <- matrix( E[[4]] , nrow(y), 1)
  obs.I.up  <- matrix(E[[5]] , nrow(y), 1)

  obs.E.dn  <- matrix(E[[6]], nrow(y), 1)
  obs.W.dn  <- matrix(E[[7]] , nrow(y), 1)
  obs.WC.dn <- matrix(E[[8]] , nrow(y), 1)
  obs.M.dn  <- matrix( E[[9]] , nrow(y), 1)
  obs.I.dn  <- matrix(E[[10]] , nrow(y), 1)
  dire <- E[[11]]

  rm(E)
  # simplified version (number of chunks)
  nchunks <- 10  ## still memory problem
  chunklength <- floor(nrow(y) / nchunks)

  endpos <- c(1: (nchunks-1) * chunklength, nrow(y))   ##changed
  begpos <- c(1, endpos[-length(endpos)] + 1)

  exp.E.up  <- ff("exp.E.up",  vmode = "double", dim = c(nrow(y), niter))
  exp.W.up  <- ff("exp.W.up",  vmode = "double", dim = c(nrow(y), niter))
  exp.WC.up <- ff("exp.WC.up", vmode = "double", dim = c(nrow(y), niter))
  exp.M.up  <- ff("exp.M.up",  vmode = "double", dim = c(nrow(y), niter))
  exp.I.up  <- ff("exp.I.up",  vmode = "double", dim = c(nrow(y), niter))
  exp.E.dn  <- ff("exp.E.dn",  vmode = "double", dim = c(nrow(y), niter))
  exp.W.dn  <- ff("exp.W.dn",  vmode = "double", dim = c(nrow(y), niter))
  exp.WC.dn <- ff("exp.WC.dn", vmode = "double", dim = c(nrow(y), niter))
  exp.M.dn  <- ff("exp.M.dn",  vmode = "double", dim = c(nrow(y), niter))
  exp.I.dn  <- ff("exp.I.dn",  vmode = "double", dim = c(nrow(y), niter))

  set.seed(seed)
  x.niter <- t(replicate(niter, sample(x)))

  ffmatrices <- c("exp.E.up", "exp.W.up", "exp.WC.up", "exp.M.up", "exp.I.up",
                  "exp.E.dn", "exp.W.dn", "exp.WC.dn", "exp.M.dn", "exp.I.dn")

  obsvecs <- c("obs.E.up", "obs.W.up", "obs.WC.up", "obs.M.up", "obs.I.up",
               "obs.E.dn", "obs.W.dn", "obs.WC.dn", "obs.M.dn", "obs.I.dn")

#  compvec <- c("<", "<", "<", "<", "<", "<", ">", ">", ">", ">")

#  outnames <- c("raw1.up", "raw2.up", "raw3.up", "raw4.up", "raw5.up",
#                "raw1.dn", "raw2.dn", "raw3.dn", "raw4.dn", "raw5.dn")

#  raw1.up <- raw2.up <- raw3.up <- raw4.up <- raw5.up <-
#    raw1.dn <- raw2.dn <- raw3.dn <- raw4.dn <- raw5.dn <-
#    1:nrow(y) # one element for each row

raw.count.up <- raw.count.dn <- matrix(0, nrow(y),5)

  for (ichunk in seq(along = begpos)){
    begchunk <- begpos[ichunk]
    endchunk <- endpos[ichunk]
    suby <- y[begchunk:endchunk,]


    for (jmat in 1:nmats){
      # chunk the matrix of permutations
      jbegmat <- begmat[jmat]
      jendmat <- endmat[jmat]
      ncolmat <- jendmat - jbegmat + 1
      subx.niter <- x.niter[jbegmat:jendmat,]

      res <- apply(subx.niter, 1, function(x) IsoGenem(x = factor(x), y = suby))

      # write results to ff matrices
      exp.E.up[begchunk:endchunk,jbegmat:jendmat] <-
        matrix(sapply(res, function(x) x[[1]]), length(begchunk:endchunk), length(jbegmat:jendmat))      ## changed
      exp.W.up[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[2]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.WC.up[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[3]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.M.up[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[4]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.I.up[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[5]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.E.dn[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[6]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.W.dn[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[7]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.WC.dn[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[8]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.M.dn[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[9]]) , length(begchunk:endchunk), length(jbegmat:jendmat))
      exp.I.dn[begchunk:endchunk,jbegmat:jendmat] <-
         matrix(sapply(res, function(x) x[[10]]), length(begchunk:endchunk), length(jbegmat:jendmat))
    }


##problem of numeric for matrix exp

#    for (i in seq(along = ffmatrices)){
#      lhs <- as.list(quote(OBSVEC[begchunk:endchunk]))
#      lhs[[2]] <- as.name(obsvecs[i])
#      rhs <- as.list(quote(FFMAT[begchunk:endchunk,]))
#      rhs[[2]] <- as.name(ffmatrices[i])
#      rsarg <- as.call(list(as.name(compvec[i]), as.call(lhs), as.call(rhs)))
#      RHS <- as.call(list(as.name("rowSums"), rsarg))

#      LHS <- as.list(quote(OUTNAME[begchunk:endchunk]))
#      LHS[[2]] <- as.name(outnames[i])

#      finalcall <- as.call(list(as.name("<-"), as.call(LHS), as.call(RHS)))

#      eval(finalcall)


##changed to
      for (i in 1:6){
      x1<-obsvecs[i]
      x2<-ffmatrices[i]
      a1<-get(x1)[begchunk:endchunk]
      a2<-matrix(as.numeric(get(x2)[begchunk:endchunk,]),byrow=FALSE,nrow=length(begchunk:endchunk),ncol=niter)
      tr<-rowSums(a1<a2)
      if (i <=5) {raw.count.up[begchunk:endchunk,i] <- tr} else {
      raw.count.dn[begchunk:endchunk,1] <- tr }
       }

      for (i in 7:10) {
      x1<-obsvecs[i]
      x2<-ffmatrices[i]
      a1<-get(x1)[begchunk:endchunk]
      a2<-matrix(as.numeric(get(x2)[begchunk:endchunk,]),byrow=FALSE,nrow=length(begchunk:endchunk),ncol=niter)
      tr<-rowSums(a1>a2)
      raw.count.dn[begchunk:endchunk,i-5] <- tr
      }


  } ##loop of ichunk

  rm(suby)
  rm(res)
  rm(subx.niter)

  raw.p.up <- data.frame(raw.count.up/niter)

  raw.p.dn <- data.frame(raw.count.dn/niter)

  rawp.up <- data.frame(row.names(y), raw.p.up)
  rawp.dn <- data.frame(row.names(y), raw.p.dn)

  raw.p.one <- data.frame(row.names(y),
                          apply(cbind(raw.p.up[,1], raw.p.dn[,1]), 1, min),
                          apply(cbind(raw.p.up[,2], raw.p.dn[,2]), 1, min),
                          apply(cbind(raw.p.up[,3], raw.p.dn[,3]), 1, min),
                          apply(cbind(raw.p.up[,4], raw.p.dn[,4]), 1, min),
                          apply(cbind(raw.p.up[,5], raw.p.dn[,5]), 1, min))
                                        # data frame will have horrible names

  raw.p.two <- raw.p.one
  raw.p.two[,2:6] <- 2 * (raw.p.one[,2:6])
  #  raw.p.two[,2:6] <- sapply(raw.p.two[,2:6],
  #                            function(x) x[x > 1] <- 1) # names ?

  raw.p.two[,2:6][raw.p.two[,2:6]>1] <- 1


  colnames(raw.p.one) <- colnames(raw.p.two) <-
    colnames(rawp.up ) <-  colnames(rawp.dn) <-
      c("Probe.ID","E2","Williams","Marcus","M","ModM")

  # put into list; bad use of return function (complaints)
  res <- list(raw.p.one = raw.p.one,
              raw.p.two = raw.p.two,
              rawp.up = rawp.up,
              rawp.dn = rawp.dn)
  rm(exp.E.up, exp.W.up, exp.WC.up, exp.M.up, exp.I.up,
     exp.E.dn, exp.W.dn, exp.WC.dn, exp.M.dn, exp.I.dn)
  gc()

  return(res)
}

