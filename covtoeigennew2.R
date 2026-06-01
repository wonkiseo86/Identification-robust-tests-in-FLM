covtoeigen = function (fdobj1, fdobj2, nharm1, nharm2 , harmfdPar1=fdPar(fdobj1), harmfdPar2 = fdPar(fdobj2), pfac = NULL) {
  coef1 = fdobj1$coefs  
  coef2 = fdobj2$coefs
  
  nlengthobs = dim(coef1)[2]
  
  eigvalc = list()
  ctemp = ((coef2)%*%t(coef1) / nlengthobs) 
  
  basisobj1 = fdobj1$basis
  basisobj2 = fdobj2$basis
  
  
  harmbasis1 <- harmfdPar1$fd$basis
  nbasis1 <- basisobj1$nbasis
  nhbasis1 <- harmbasis1$nbasis
  Lfdobj1 <- harmfdPar1$Lfd
  lambda1 <- harmfdPar1$lambda
  Lmat1 <- eval.penalty(harmbasis1, 0)              ##matrix consisting of inner product of basis functions used for approximating functional data and harmonic functions
  if (lambda1 > 0) {
    Rmat1 <- eval.penalty(harmbasis1, Lfdobj1)
    Lmat1 <- Lmat1 + lambda1 * Rmat1
  }
  
  harmbasis2 <- harmfdPar2$fd$basis
  nbasis2 <- basisobj2$nbasis
  nhbasis2 <- harmbasis2$nbasis
  Lfdobj2 <- harmfdPar2$Lfd
  lambda2 <- harmfdPar2$lambda
  Lmat2 <- eval.penalty(harmbasis2, 0)              ##matrix consisting of inner product of basis functions used for approximating functional data and harmonic functions
  if (lambda2 > 0) {
    Rmat2 <- eval.penalty(harmbasis2, Lfdobj2)
    Lmat2 <- Lmat2 + lambda2 * Rmat2
  }
  
  
  Lmat1 <- (Lmat1 + t(Lmat1))/2
  Mmat1 <- chol(Lmat1)
  Mmatinv1 <- solve(Mmat1)
  Jmat1 = inprod(harmbasis1, basisobj1)
  MIJW1 = crossprod(Mmatinv1, Jmat1)
  
  Lmat2 <- (Lmat2 + t(Lmat2))/2
  Mmat2 <- chol(Lmat2)
  Mmatinv2 <- solve(Mmat2)
  Jmat2 = inprod(harmbasis2, basisobj2)
  MIJW2 = crossprod(Mmatinv2, Jmat2)
  
  Wmat1 = t(ctemp)%*%crossprod(MIJW2)%*%ctemp
  Cmat1 = MIJW1 %*% Wmat1 %*% t(MIJW1)
  Cmat1 <- (Cmat1 + t(Cmat1))/2
  result1 <- eigen(Cmat1)
  
  Wmat2 = (ctemp)%*%crossprod(MIJW1)%*%t(ctemp)
  Cmat2 = MIJW2 %*% Wmat2 %*% t(MIJW2)
  Cmat2 <- (Cmat2 + t(Cmat2))/2
  result2 <- eigen(Cmat2)
  
  
  if (is.null(pfac)==TRUE){
    eigvecc1 <- as.matrix(result1$vectors[, 1:nharm1])
    sumvecc1 <- apply(eigvecc1, 2, sum)
    eigvecc1[, sumvecc1 < 0] <- -eigvecc1[, sumvecc1 < 0]
    result1$values[which(result1$values==0)] = .9*min(abs(result1$values[which(result1$values>0)]))
    result1$values[which(result1$values<0)] = -result1$values[which(result1$values<0)]
    eigvecc2 <- MIJW2%*%(ctemp)%*%t(MIJW1)%*%eigvecc1%*%diag(1/sqrt(result1$values[1:nharm1]))  
   }else {
     eigvecc2 <- as.matrix(result2$vectors[, 1:nharm2])
     sumvecc2 <- apply(eigvecc2, 2, sum)
     eigvecc2[, sumvecc2 < 0] <- -eigvecc2[, sumvecc2 < 0]
     result2$values[which(result1$values==0)] = .9*min(abs(result2$values[which(result2$values>0)]))
     result2$values[which(result2$values<0)] = -result2$values[which(result2$values<0)]
     eigvecc1 <- MIJW1%*%t(ctemp)%*%t(MIJW2)%*%eigvecc2%*%diag(1/sqrt(result2$values[1:nharm2])) 
  }
  
 
  harmcoef1 <- Mmatinv1 %*% eigvecc1  
  harmcoef2 <- Mmatinv2 %*% eigvecc2 
  harmfd1 <- fd(harmcoef1, harmbasis1)
  harmfd2 <- fd(harmcoef2, harmbasis2)
  
  eigvalc <-  (diag(t(eigvecc1)%*% MIJW1%*%t(ctemp)%*%t(MIJW2)%*%eigvecc2))
  
  covtoeigen <- list(harmfd1,harmfd2,eigvalc)
  
  names(covtoeigen) <- c("harmonics1", "harmonics2", "values")
  return(covtoeigen)
  
}



