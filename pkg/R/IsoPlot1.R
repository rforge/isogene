IsoPlot1 <- function (x, y, type = c("Continuous", "Ordinal")) 
{    options(warn=-1)
     if (is.na(match(type, c("Ordinal", "Continuous")))) 
        {options(warn=0)
	   print("Error: The dose can be only Continuous or Ordinal")}
     
     else{
     y1 <- as.numeric(y)[order(x)]
     m.y <- tapply(y1, as.factor(x[order(x)]), mean) 
          
     if ( missing(type) | type== "Continuous" ){
     plot(x[order(x)], y1, lwd=2, xlab = "Doses", ylab = "Gene Expression")
     points(unique(x), m.y[order(unique(x))], pch = "+", cex = 0.85) }
    
     
     if (type == "Ordinal") {
	     miny <- min(y)
       maxy <- max(y)
       unx <- sort(unique(x))
    	 catx <- factor(x , levels = sort(unique.default(x )), labels=unx ,ordered =TRUE)
    	 a <-c(1:length(unx))
    	 plot(a ,  ylim=c(miny,maxy), pch="", ylab = "Gene Expression", xlab ="Doses",axes = FALSE)
    	 axis(1, sort(unique(a)),  as.character(unx))
    	 axis(2)
    	 points(catx,y, lwd=2)
       points(unique(catx), m.y[order(unique(x))], pch = "+", cex = 0.85)

       }

     title(paste("Gene: ", row.names(y), sep = ""))
     
     }
}
