# GLOBAL VARIABLES: #

# PROCEDURE CODE: #
bandfn <- function(u,kernel){

 
  tu <- nrow(u)
  p <- ncol(u)
  e <- u
  te <- tu
  eb <- as.matrix(e[1:(te-1),])
  ef <- as.matrix(e[2:te,])
  ae <- as.matrix(colSums(eb*ef)/colSums(eb^2))
  ee <- ef - eb*(matrix(1,nrow(eb),1)%*%t(ae))
  se <- as.matrix(colMeans(ee^2))
  ad <- sum((se/((1-ae)^2))^2)
  a1 <- 4*sum((ae*se/(((1-ae)^3)*(1+ae)))^2)/ad
  a2 <- 4*sum((ae*se/((1-ae)^4))^2)/ad
  if (kernel == 2){                   # Quadratic Spectral #
     bandw <- 1.3221*((a2*te)^.2)
  }   
  if (kernel == 1){                   # Parzen #
     bandw <- 2.6614*((a2*te)^.2)               
     if (bandw > (te-2)) bandw <- te-2
  }
  if (kernel == 7){                   # Bartlett #   
     bandw <- 1.1447*((a1*te)^.333)
     if (bandw > (te-2)) bandw <- te-2
  }

  return(bandw)
}
