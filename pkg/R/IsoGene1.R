IsoGene1 <- function(x, y){
  ordx <- order(x)
  x <- x[ordx]
  y <- y[ordx]  # reverse order (sorted two times)

  unx <- unique(x)
      
  y.m <- tapply(y, as.factor(x), mean)
  y.m.tot <- rep(mean(y), length(unx))
  y.is.u <- isoreg(unx,  y.m)$yf
  y.is.d <- rev(isoreg(rev(unx), y.m)$yf)
  n.p <- table(x)
  n.g <- length(n.p)

  iso.u <- rep.iso.u <- rep.iso.d <- y.m.all <- NULL

  ##################################################
   
  rep.iso.d <- rep(y.is.d, n.p)
  rep.iso.u <- rep(y.is.u, n.p)
  y.m.all <- rep(y.m, n.p)

  ########################################################


  SST0 <- sum((y - mean(y))^2)

  SSIS.u1<-sum((rep.iso.u-y)^2)
  SSIS.d1<-sum((rep.iso.d-y)^2)

  SST<-sum((y-y.m.all)^2)

  direction <- if (SSIS.u1 <= SSIS.d1) "u" else "d" # no ifelse
         
  lambda1.up <- SSIS.u1 / SST0 # no brackets
  Esquare.up <- 1 - lambda1.up
  iso.u <- y.is.u
         
  ##################################################
    
  w.up <- (y.is.u[n.g]-y.m[1])/ sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))
  w.c.up <- (y.is.u[n.g]-y.is.u[1])/ sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))
  m.up <- (y.is.u[n.g]-y.is.u[1])/sqrt(SSIS.u1/(sum(n.p)-n.g))
  i.up <- (y.is.u[n.g]-y.is.u[1])/sqrt(SSIS.u1/(sum(n.p)-length(unique(y.is.u))))
  #  dfi=length(unique(y.is.u)) 

  lambda1.dn <- SSIS.d1 / SST0 # no brackets needed
  Esquare.dn <- 1 - lambda1.dn
  iso.u <- y.is.d

  ##################################################
   
  w.dn <- (y.is.d[n.g]-y.m[1])/ sqrt(2*SST / (sum(n.p)-n.g)/(n.g-1))
  w.c.dn <- (y.is.d[n.g]-y.is.d[1]) / sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))
  m.dn <- (y.is.d[n.g]-y.is.d[1]) / sqrt(SSIS.d1/(sum(n.p)-n.g))
  i.dn <- (y.is.d[n.g]-y.is.d[1]) / sqrt(SSIS.d1/(sum(n.p)-length(unique(y.is.d))))
  # dfi=length(unique(y.is.d))   

            ## do not return a list directly like that
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
