#pava <- function (x, wt=rep(1,length(x)))
#  Compute the isotonic regression of numeric vector 'x', with
#  weights 'wt', with respect to simple order.  The pool-adjacent-
#  violators algorithm is used.  Returns a vector of the same length
#  as 'x' containing the regression.

#  02 Sep 1994 / R.F. Raubertas
#{
#   n <- length(x)
#   if (n <= 1) return (x)
#   if (any(is.na(x)) || any(is.na(wt))) {
#      stop ("Missing values in 'x' or 'wt' not allowed")
#   }
#   lvlsets <- (1:n)
#   repeat {
#      viol <- (as.vector(diff(x)) < 0)  # Find adjacent violators
#      if (!(any(viol))) break
#
#      i <- min( (1:(n-1))[viol])        # Pool first pair of violators
#      lvl1 <- lvlsets[i]
#      lvl2 <- lvlsets[i+1]
#      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
#      x[ilvl] <- sum(x[ilvl]*wt[ilvl]) / sum(wt[ilvl])
#      lvlsets[ilvl] <- lvl1
#   }
#  return(x)
#}


IsoGene1 <- function(x, y){
  ordx <- order(x)
  x <- x[ordx]
  y <- y[ordx]

  unx <- unique(x)
  
  n.p <- table(x)
  n.g <- length(n.p)           
     
  y.m <- tapply(y, as.factor(x), mean)
  y.m.tot <- rep(mean(y), length(unx))
  
   y.is.u <- pava(y = y.m, w = n.p )
   y.is.d <- rev(pava(y = rev(y.m), w = rev(n.p)))
##  y.is.u <- isoreg(unx,  y.m)$yf          ##equal weights
##  y.is.d <- rev(isoreg(rev(unx), y.m)$yf)
  


  iso.u <- rep.iso.u <- rep.iso.d <- y.m.all <- NULL

  rep.iso.d <- rep(y.is.d, n.p)
  rep.iso.u <- rep(y.is.u, n.p)
  y.m.all <- rep(y.m, n.p)

  SST0 <- sum((y - mean(y))^2)

  SSIS.u1<-sum((rep.iso.u-y)^2)
  SSIS.d1<-sum((rep.iso.d-y)^2)

  SST<-sum((y-y.m.all)^2)

  direction <- if (SSIS.u1 <= SSIS.d1) "u" else "d" # no ifelse
         
  lambda1.up <- SSIS.u1 / SST0 # no brackets
  Esquare.up <- 1 - lambda1.up
  iso.u <- y.is.u
    
  w.up <- (y.is.u[n.g]-y.m[1])/ sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))
  w.c.up <- (y.is.u[n.g]-y.is.u[1])/ sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))
  m.up <- (y.is.u[n.g]-y.is.u[1])/sqrt(SSIS.u1/(sum(n.p)-n.g))
  i.up <- (y.is.u[n.g]-y.is.u[1])/sqrt(SSIS.u1/(sum(n.p)-length(unique(y.is.u))))

  lambda1.dn <- SSIS.d1 / SST0
  Esquare.dn <- 1 - lambda1.dn
  iso.u <- y.is.d
  
  n.pSum <- sum(n.p)
  w.dn <- (y.is.d[n.g]-y.m[1])/ sqrt(2*SST / (n.pSum-n.g)/(n.g-1))
  w.c.dn <- (y.is.d[n.g]-y.is.d[1]) / sqrt(2*SST/(n.pSum-n.g)/(n.g-1))
  m.dn <- (y.is.d[n.g]-y.is.d[1]) / sqrt(SSIS.d1/(n.pSum-n.g))
  i.dn <- (y.is.d[n.g]-y.is.d[1]) / sqrt(SSIS.d1/(n.pSum-length(unique(y.is.d))))
  
  res <-  list(E2.up = Esquare.up,
               Williams.up = as.numeric(w.up),
               Marcus.up = as.numeric(w.c.up),
               M.up = as.numeric(m.up),
               ModM.up = as.numeric(i.up),
               E2.dn = Esquare.dn,
               Williams.dn = as.numeric(w.dn),
               Marcus.dn = as.numeric(w.c.dn),
               M.dn = as.numeric(m.dn),
               ModM.dn = as.numeric(i.dn),
               direction = direction)
  return(res)
}
