IsopvaluePlot <-
function(x,y,niter,seed,stat=c("E2","Williams","Marcus","M","ModifM") ){

      Probe.ID <- row.names(y)
      obs <- IsoGene1(x,y)

      if(stat=="E2") obs.up=obs[[1]]
      if(stat=="Williams") obs.up=obs[[2]]
      if(stat=="Marcus") obs.up=obs[[3]]
      if(stat=="M") obs.up=obs[[4]]
      if(stat=="ModifM") obs.up=obs[[5]]
      
      if(stat=="E2") obs.dn=obs[[6]]
      if(stat=="Williams") obs.dn=obs[[7]]
      if(stat=="Marcus") obs.dn=obs[[8]]
      if(stat=="M") obs.dn=obs[[9]]
      if(stat=="ModifM") obs.dn=obs[[10]]
      
      exp.up <- exp.dn <- c(1:niter)
    
      set.seed(seed)

      ### replicate ?
      
      x.niter <- t(sapply(1:niter, function(i) sample(x)))
      # dim(x.niter) ## of no use
    
    
    for ( j in 1: niter) {
      exps <- IsoGene1(x.niter[j,], y)
#     exps <- lam.niter # why two assignments needed ?
       
      if(stat=="E2") exp.up[j] = exps[[1]]
      if(stat=="Williams") exp.up[j] = exps[[2]]
      if(stat=="Marcus") exp.up[j] = exps[[3]]
      if(stat=="M") exp.up[j] = exps[[4]]
      if(stat=="ModifM") exp.up[j] = exps[[5]]
      
      if(stat=="E2") exp.dn[j] = exps[[6]]
      if(stat=="Williams") exp.dn[j] = exps[[7]]
      if(stat=="Marcus") exp.dn[j] = exps[[8]]
      if(stat=="M") exp.dn[j] = exps[[9]]
      if(stat=="ModifM") exp.dn[j] = exps[[10]]
            
       cat(paste(j, ". "))
    }
    rawp.up <- sum(obs.up < exp.up) / niter
    rawp.dn <- sum(obs.dn > exp.dn) / niter
    if (stat == "E2") rawp.dn <- sum(obs.dn<exp.dn)/niter # no braces needed
      
    par(mfrow = c(2,1))
            
    hist(exp.up,main = "", nclass = 1000, col = 0, probability = TRUE, # always use TRUE (not T)
         xlim=c(min(exp.up, obs.up), max(exp.up, obs.up)), xlab = paste(stat))
    dx <- density(exp.up, from = min(exp.up), to = max(exp.up))
    abline(v = obs.up, col = 7, lwd = 3)
    #  lines(c(obs.up,obs.up),c(0,max(dx$y/2)),lwd=3,col=7)
    lines(dx$x,dx$y, lwd = 3, col = 5)
    title(paste("Gene",Probe.ID,":p-value^{up}=",rawp.up,sep='')) 
  
    hist(exp.dn, main = "", nclass = 1000, col = 0, probability = TRUE,
         xlim = c(min(exp.dn,obs.dn), max(exp.dn,obs.dn)), xlab = paste(stat))
    dx <- density(exp.dn, from = min(exp.dn), to = max(exp.dn))
    abline(v = obs.dn, col = 7, lwd = 3)
    # lines(c(obs.dn,obs.dn),c(0,max(dx$y/2)),lwd=3,col=7)
    lines(dx$x, dx$y, lwd = 3, col = 5)
    title(paste("Gene", Probe.ID, ":p-value^{down}=", rawp.dn, sep = '')) 
}
