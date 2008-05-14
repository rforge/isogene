IsoGenemSAM <-
function(x, y, fudge.factor)
{
 y <- as.matrix(y)

  if(dim(y)[[1]]==1){
    print("Warning: There is only one gene in the data set.")
  } else {  
  ordx <- order(x) # compute once
  x <- x[ordx]
  y <- y[,ordx]  # reverse order (sorted two times)

  unx <- unique(x) # compute once
  
  ydf <- as.data.frame(t(y))
  y.m <- do.call("cbind", unclass(by(ydf, x, mean)))
  
  y.m.tot <- matrix(rep(rowMeans(y), length(x)), ncol = length(x))
  
  y.is.u <- t(apply(y.m, 1, function(x) isoreg(unx, x)$yf))
  y.is.d <- t(apply(y.m, 1, function(x) rev(isoreg(rev(unx), x)$yf)))

  n.p <- table(x)
  n.g <- length(n.p)

  ###################################################
   
  rep.iso.d <- y.is.d[, rep(1:length(n.p),n.p)]
  rep.iso.u <- y.is.u[, rep(1:length(n.p),n.p)]
  
  y.m.all <- y.m[, rep(1:length(n.p), n.p)]

  ########################################################


  SST0 <- rowSums((y - rowMeans(y))^2)

  SSIS.u1 <- rowSums((rep.iso.u - y)^2)
  SSIS.d1 <- rowSums((rep.iso.d - y)^2)

  SST <- rowSums((y - y.m.all)^2)

  direction=NULL

  direction <- ifelse(SSIS.u1 <= SSIS.d1, "u", "d")
           
  iso.u = y.m
  SSIS.dir= NULL

  iso.u[direction=="u",] <- y.is.u[direction=="u",]
  iso.u[direction=="d",] <- y.is.d[direction=="d",]

  SSIS.dir[direction=="u"] <- SSIS.u1[direction=="u"]
  SSIS.dir[direction=="d"] <- SSIS.d1[direction=="d"]

     
  Esquare <- (SST0-SSIS.dir)/(SST0+fudge.factor[[1]])
  w <- (iso.u[,n.g] - y.m[,1]) / (sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1)) +fudge.factor[[2]])
  w.c <- (iso.u[,n.g] - iso.u[,1]) / (sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1)) +fudge.factor[[3]])
  m <- (iso.u[,n.g] - iso.u[,1]) / (sqrt(SSIS.dir/(sum(n.p)-n.g)) +fudge.factor[[4]])
  i <- (iso.u[,n.g] - iso.u[,1]) / (sqrt(SSIS.dir/(sum(n.p) - apply(iso.u, 1, function(x) length(unique(x))))) + fudge.factor[[5]])



  res <-  list(E2 = Esquare,
                Williams = as.numeric(w),
                Marcus = as.numeric(w.c),
                M = as.numeric(m),
                ModM = as.numeric(i),
                direction = direction)
    return(res)
  }
}
