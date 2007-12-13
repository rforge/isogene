IsoTestBH <-                        ## naming of the function (BH while all sorts of tests)
function (rp, FDR, type =  c("BH", "BY"),
          stat = c("E2", "Williams", "Marcus", "M", "ModifM")){

  Probe.ID <- rp[,1]

  type <- match.arg(type)
  stat <- match.arg(stat)
    
  if (stat == "E2") rpraw <- rp[,2]
  else if (stat == "Williams") rpraw <- rp[,3]
  else if (stat == "Marcus") rpraw <- rp[,4]
  else if (stat == "M") rpraw <- rp[,5]
  else rpraw <- rp[,6]
    # (stat == "ModifM") rpraw = rp[,6]

  # library(multtest) ## never use inside a function
  procs <- c("BH", "BY") # only "BH" and "BY" are needed
  res <- mt.rawp2adjp(rpraw, procs)
  adjp <- res$adjp[order(res$index), ]

  # use names Bonferroni   Holm Hochberg     SidakSS     SidakSD          BH
   
  if (type == "BH")  {
    place.keep33 <- which(adjp[,2] <= FDR)
  } else { # type == "BY"  
    place.keep33 <- which(adjp[,3] <= FDR)
  }


  sign.Probe.ID <- Probe.ID[place.keep33] ### keep3 will not be available if
                                   ### another type than BH or BY is chosen (so
                                   ### removed other types
  
  if (type == "BH")  {
  sign.genes <- data.frame(sign.Probe.ID,
                           place.keep33,
                           adjp[adjp[,2] <= FDR,1],
                           adjp[adjp[,2] <= FDR,2])
       } else {
  sign.genes <- data.frame(sign.Probe.ID,
                           place.keep33,
                           adjp[adjp[,3] <= FDR,1],
                           adjp[adjp[,3] <= FDR,3])     
    }    
  
    
  names(sign.genes) <- c("Probe.ID", "row.name", "raw p-values",
                         paste(type, "adjusted p values", sep = " "))
  return(sign.genes)
}
