 


lr_var <- function(u, kernel = 2, aband = 1, bband = 1/3, urprint = 1, white = 0 ){
  tu <- nrow(u)
  p <- ncol(u)
  bandw = bband
  if (white == 1){
    te <- tu-1
    au <- qr.solve(as.matrix(u[1:te,]),as.matrix(u[2:tu,]))
    e <- as.matrix(u[2:tu,]) - as.matrix(u[1:te,])%*%au
  }else{
    e <- u
    te <- tu
  }
  
  # Estimate Covariances #
  if (bandw >0){
    if (kernel == 1){                             # Parzen kernel #
      tm <- floor(bandw)
      if (tm > 0){
        jb <- as.matrix(seq(1,tm,1)/bandw)
        kern <- (1 - (jb^2)*6 + (jb^3)*6)*(jb <= .5)
        kern <- kern + ((1-jb)^3)*(jb > .5)*2
      }
      intker = 0.75
    } else if (kernel == 2){                       ####Biweight
      tm <- floor(bandw)
      if (tm >0){
        jb = as.matrix(seq(1,tm,1)/bandw)
        #   kern <- .75*(1-jb^2)
        kern <- (1-jb^2)^2
      }
      intker = 16/15
    }else if (kernel == 3){                 ####Triweight
      tm <- floor(bandw)
      if (tm >0){
        jb = as.matrix(seq(1,tm,1)/bandw) 
        kern <- (1-jb^2)^3
      }
      intker = 32/35
    }else if (kernel == 4){                 ####Cosine
      tm <- floor(bandw)
      if (tm >0){
        jb = as.matrix(seq(1,tm,1)/bandw) 
        kern <-  cos(jb*pi/2)  
      }
      intker = 4/pi
    }else if (kernel == 5){                ####Tukey Hanning 
      tm <- floor(bandw)
      if (tm >0){
        jb = as.matrix(seq(1,tm,1)/bandw) 
        kern <-  (1+cos(jb*pi))/2  
      }
      intker = 1
    }else if (kernel == 6){               ###quadweight   
       tm <- floor(bandw)
      if (tm >0){
        jb = as.matrix(seq(1,tm,1)/bandw) 
        kern <- (1-jb^2)^4
      }
      intker = 256/315
    }else if (kernel == 7){                             # Bartlett kernel #
    tm <- floor(bandw)
    if (tm > 0) kern <- as.matrix(1 - seq(1,tm,1)/bandw)
	  intker = 1
    }
  
    
    
    lam <- matrix(0,p,p)
    for (j in 1:tm){
      kj <- kern[j]
      lam <- lam + (t(as.matrix(e[1:(te-j),]))%*%as.matrix(e[(1+j):te,]))*kj
    }
    omega <- (t(e)%*%e + lam + t(lam)) 
    
    if (white == 1){
      eau <- solve(diag(p) - au)
      omega <- t(eau)%*%omega%*%eau
    }
  } else {
    omega = crossprod(e)
    intker = 0
  }
  list(omega=omega,bandw=bandw,intker = intker)
}
