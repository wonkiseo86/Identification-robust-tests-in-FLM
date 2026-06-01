## Description: Simulation codes for "Identification-Robust Testing in Endogenous Functional Linear Regression with Weak or Irrelevant Auxiliary Variables".
## The code can produce the simulation results for the model WITH intercept.

## Loading packages
library(fda);library(sde);library(truncnorm)

## Setwd and functions
setwd("Set folder containing all the files")
inner = dget("inprod.R")
operator = dget("operator.R")
lrvar = dget("lr_var_v2.R")
bandfn = dget("bandfn.R")

consts=function(p){a=1/(2*p+3);b=2/(p+2);1/sqrt(a*b)}

## Integration weights
make_int_weights <- function(grid) {
  dx <- grid[2] - grid[1]
  w <- rep(dx, length(grid))
  w[1] <- dx / 2
  w[length(grid)] <- dx / 2
  w
}

## Inner products between columns of B and vector x
mat_inner <- function(B, x, w) {as.numeric(crossprod(B, x * w))}

## L2 norms of columns
l2norm_cols <- function(M, w) {colSums(M^2 * w)}

## Coefficient matrix for all finite p values
make_W <- function(PSET, nn) {
  W <- sapply(PSET, function(pp) {
    const_fn(pp, nn) *
      vapply(seq_len(nn), function(j) coef_fn(pp, j, nn), numeric(1))})
  colnames(W) <- paste0("p", PSET)
  W
}

## Compute endpoint statistic and all finite-p statistics
calc_stats <- function(ldata, W, w, nobs) {
  ## p = infinity statistic
  lm <- rowMeans(ldata)
  stat_inf <- nobs * sum(w * lm^2)

  ## finite-p statistics
  LL <- ldata %*% W
  stat_p <- nobs * l2norm_cols(LL, w)

  c(inf = stat_inf, stat_p)
}

## Functions for simulating critical values
calc_cv2 <- function(eigen1, eigen2, band2, cviter, alpha = 0.95) {
  band2 <- min(band2, length(eigen1), length(eigen2))
  idx <- seq_len(band2)

  atem <- matrix(rnorm(cviter * band2)^2, nrow = cviter, ncol = band2)
  CV <- atem %*% cbind(eigen1[idx], eigen2[idx])

  c(cv1 = as.numeric(quantile(CV[, 1], alpha, names = FALSE)),
    cv2 = as.numeric(quantile(CV[, 2], alpha, names = FALSE)))
}

## 1 = informative IV, 2 = weakly informative IV, X = exogeneity-based benchmark Z_t=X_t
calc_cv3 <- function(eigen1, eigen2, eigenX, band2, cviter, alpha = 0.95) {
  band2 <- min(band2, length(eigen1), length(eigen2), length(eigenX))
  idx <- seq_len(band2)

  atem <- matrix(rnorm(cviter * band2)^2, nrow = cviter, ncol = band2)

  CV <- atem %*% cbind(eigen1[idx], eigen2[idx], eigenX[idx])

  c(cv1 = as.numeric(quantile(CV[, 1], alpha, names = FALSE)),
    cv2 = as.numeric(quantile(CV[, 2], alpha, names = FALSE)),
    cvX = as.numeric(quantile(CV[, 3], alpha, names = FALSE)))
}

## Precalcuated coefficients functions used for constructing test statistics
coef_fn=function(pp,j,nn)
{
	if(pp==0){aa=(nn-j+1)/(nn^2)}
	if(pp==1){aa=(-j + nn +1)*(j+ nn )/(2*nn^3)}
	if(pp==2){aa=(-j + nn +1)*(2*j^2 + 2*j*nn-j+2*nn^2 + nn)/(6*nn^4)}
	if(pp==3){aa=(-j + nn +1)*(j+nn)*(j^2-j + nn^2 + nn)/(4*nn^5)}
	if(pp==4){aa=(-6*j^5 +15*j^4 - 10*j^3 + j + 6*nn^5 +15*nn^4 + 10*nn^3 - nn)/(30*nn^6)}
	if(pp==5){aa=(-2*j^6 +6*j^5 - 5*j^4 +j^2 + (nn^2)*((nn+1)^2)*(2*(nn^2)+2*nn-1))/(12*nn^7)}
	if(pp==6){aa=(-6*j^7 +21*j^6 - 21*j^5 + 7*j^3 - j + 6*nn^7 + 21*nn^6 + 21*nn^5 - 7*nn^3 + nn)/(42*nn^8)}
	if(pp==7){aa=(-3*j^8 +12*j^7 - 14*j^6 + 7*j^4 - 2*j^2 + 3*nn^8 + 12*nn^7 + 14*nn^6 - 7*nn^4 + 2*nn^2)/(24*nn^9)}
	if(pp==8){aa=(-10*j^9 +45*j^8 - 60*j^7 + 42*j^5 - 20*j^3 + 3*j + 10*nn^9 + 45*nn^8 + 60*nn^7 - 42*nn^5 + 20*nn^3 - 3*nn)/(90*nn^10)}
    if(pp==9){aa=(-2*j^10 +10*j^9 - 15*j^8 + 14*j^6 - 10*j^4 + 3*j^2 + 2*nn^10 + 10*nn^9 + 15*nn^8 - 14*nn^6 + 10*nn^4 - 3*nn^2)/(20*nn^11)}
	if(pp==10){aa=(-6*j^11 + 33*j^10 - 55*j^9 + 66*j^7 -66*j^5 + 33*j^3 -5*j + nn*( 6*nn^10+33*nn^9+55*nn^8-66*nn^6 + 66*nn^4-33*nn^2+5 ))/(66*nn^12)}
return(aa)
}

const_fn=function(p,nn)
{aa=0
for (j in 1:nn){aa=aa+(coef_fn(p,j,nn)^2)*nn}
return(sqrt(1/aa))
}

## Brownian Briges as functional errors
BBridge2=function(x0=0,y=0,t00=0,T0=1,N0=nt-1)
{
vvtem=BBridge(x=x0,t0=t00,T=T0,N=N0)
vvtem
}


for (upermean in c(0))   # upermean corresponds to beta_u in the paper;  To reproduce all columns, run the code with upermean = 0, 0.1, and 0.25.
{
for (kernf in c(1,7))   #This specifies the kernel function k in the paper;  1:Parzen, 7: Barttlet
{

BANDCOL1=NULL
BANDCOL2=NULL
BANDCOLX=NULL

set.seed(99999)
ceil=1
maxiter=2000 ;  maxbasis=50
TSET=c(100,200,400) 
bw2=0.333; ranbasis=4; ranini=1 ;PSET=c(7,3,1,0); pcsforbw=5;  FARmean1=0.4; FARvar1=0.4 ;  ranfirst=1 ; scalran=0.2 ; limlim=0.8; limlim2 = -0.2;  upersd= -limlim2; cviter=1000
FAR = 1; scalevar=1;  fardec=0.95 ; facz=0.25  #Note: PSET :descending order. TSET : ascending order. 
ulower=limlim2; uupper=limlim
intrange=3 ; intrange2 = 3 ;  FARmean2=0.6; FARvar2=0.6; minlim=limlim2; maxlim =limlim ; #; pp=Inf #(exp5)

WLIST <- setNames(lapply(TSET, function(nn) make_W(PSET, nn)), as.character(TSET))

####################################################
## Other key simulation paramters and related settings ###
#################################################### 
lbnumber2=100
nt = 200 ; t = (0:(nt-1))/(nt-1); w_int <- make_int_weights(t)

LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}
LBF=cbind(rep(1,length(t)),LBF)
lb=LBF
for(i in 2:lbnumber2){  
  for(j in 1:i)  { 
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber2){
  LBF[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
} 




nnbasis=50
inidiscard=50
kpoint = t
basis_fn = create.bspline.basis(rangeval = c(0,1), nbasis = nnbasis )

dimt=length(TSET)

TPOWER=array(NA,dim=c(maxiter,4,dimt)); TPOWERA=array(NA,dim=c(maxiter,4,dimt)); TPOWERB=array(NA,dim=c(maxiter,4,dimt)); TPOWERC=array(NA,dim=c(maxiter,4,dimt));  TPOWERD=array(NA,dim=c(maxiter,4,dimt))
TSTAT=array(NA,dim=c(maxiter,4,dimt)); TSTATA=array(NA,dim=c(maxiter,4,dimt)); TSTATB=array(NA,dim=c(maxiter,4,dimt)); TSTATC=array(NA,dim=c(maxiter,4,dimt)); TSTATD=array(NA,dim=c(maxiter,4,dimt));
TSIZE=array(NA,dim=c(maxiter,4,dimt));TSIZEA=array(NA,dim=c(maxiter,4,dimt)); TSIZEB=array(NA,dim=c(maxiter,4,dimt)); TSIZEC=array(NA,dim=c(maxiter,4,dimt)); TSIZED=array(NA,dim=c(maxiter,4,dimt)); 

TPOWER2=array(NA,dim=c(maxiter,4,dimt)); TPOWERA2=array(NA,dim=c(maxiter,4,dimt)); TPOWERB2=array(NA,dim=c(maxiter,4,dimt)); TPOWERC2=array(NA,dim=c(maxiter,4,dimt));  TPOWERD2=array(NA,dim=c(maxiter,4,dimt))
TSTAT2=array(NA,dim=c(maxiter,4,dimt)); TSTATA2=array(NA,dim=c(maxiter,4,dimt)); TSTATB2=array(NA,dim=c(maxiter,4,dimt)); TSTATC2=array(NA,dim=c(maxiter,4,dimt)); TSTATD2=array(NA,dim=c(maxiter,4,dimt));
TSIZE2=array(NA,dim=c(maxiter,4,dimt));TSIZEA2=array(NA,dim=c(maxiter,4,dimt)); TSIZEB2=array(NA,dim=c(maxiter,4,dimt)); TSIZEC2=array(NA,dim=c(maxiter,4,dimt)); TSIZED2=array(NA,dim=c(maxiter,4,dimt)); 

## Exogeneity-based benchmark: set Z_t = X_t, with demeaning for the intercept model
TPOWERX=array(NA,dim=c(maxiter,4,dimt)); TPOWERAX=array(NA,dim=c(maxiter,4,dimt)); TPOWERBX=array(NA,dim=c(maxiter,4,dimt)); TPOWERCX=array(NA,dim=c(maxiter,4,dimt)); TPOWERDX=array(NA,dim=c(maxiter,4,dimt))
TSTATX=array(NA,dim=c(maxiter,4,dimt)); TSTATAX=array(NA,dim=c(maxiter,4,dimt)); TSTATBX=array(NA,dim=c(maxiter,4,dimt)); TSTATCX=array(NA,dim=c(maxiter,4,dimt)); TSTATDX=array(NA,dim=c(maxiter,4,dimt));
TSIZEX=array(NA,dim=c(maxiter,4,dimt)); TSIZEAX=array(NA,dim=c(maxiter,4,dimt)); TSIZEBX=array(NA,dim=c(maxiter,4,dimt)); TSIZECX=array(NA,dim=c(maxiter,4,dimt)); TSIZEDX=array(NA,dim=c(maxiter,4,dimt));

#start_time <- Sys.time()
for (iteration in 1:maxiter)
{
#######################################################
####### Data generation ###############################
#######################################################
T=max(TSET)
eigenfactors=(1:(lbnumber2+1))^(-1)
eta = matrix(rep(NA,((T+inidiscard)*nt)),ncol=T+inidiscard) 
eta[,1] = BBridge2(x=0,y=0,t0=0,T=1,N=nt-1)
eta2 = eta; eta22 = eta; eta3 = eta; x_mat = eta; z_mat=eta; z_mat1=eta; z_mat2=eta
eta3[1,1]=rnorm(1)

##persistence setting
if(FAR==1){dcoefficients1=runif(lbnumber2+1,minlim,maxlim)*fardec^(0:(lbnumber2))}
if(FAR==0){dcoefficients1=runif(lbnumber2+1,minlim,maxlim)*fardec^(0:(lbnumber2))}

seed1=rnorm(intrange,0,1)
seed1=seed1/sqrt(sum(seed1^2))


#####################
ffindex=1:maxbasis
seed2=sample(ranini:ranbasis,1)
fixindex=sample(1:5,seed2)

LBFtem=as.matrix(LBF[,ffindex],ncol=maxbasis)
	
upersist=0 # Not used and set to zero always; used eariler for code checking
if(ranfirst==1){coeffs_first=runif(maxbasis,1+limlim2, 1-limlim2)}else{coeffs_first=rep(1,maxbasis)}

seedx=rnorm(intrange,0,1); seedx=seedx/sqrt(sum(seedx^2))
seedz=rnorm(intrange,0,1); seedz=seedz/sqrt(sum(seedz^2))
meanx=t((t(seedx))%*%t(LBF[,1:length(seedx)]))
meanz=t((t(seedz))%*%t(LBF[,1:length(seedz)]))
meany=rnorm(1,0,1)
	
  for (jiter in 2:((T+inidiscard))){
	eta2[,jiter] = sqrt(scalevar)*BBridge2(x=0,y=0,t0=0,T=1,N=nt-1); eta22[,jiter] = sqrt(scalevar)*BBridge2(x=0,y=0,t0=0,T=1,N=nt-1) ; 
	randtem=rnorm(1)
	eta3[1,jiter] =  upersist*eta3[1,jiter-1] + randtem
	
    x_signal=eta2[,jiter]  + operator(LBF[,1:length(dcoefficients1)],dcoefficients1, x_mat[,jiter-1]-meanx ,t)
    x_mat[,jiter] = meanx + upermean*eta3[1,jiter]   + x_signal
   
    scores <- mat_inner(LBFtem, x_signal, w_int)
	mask <- as.numeric(seq_len(maxbasis) %in% fixindex)

	gtem1 <- as.numeric(LBFtem %*% (coeffs_first * scores))
	gtem2 <- as.numeric(LBFtem %*% (coeffs_first * scores * mask))

	z_mat1[, jiter] <- meanz + gtem1 + facz * eta22[, jiter]  # good aux var.
	z_mat2[, jiter] <- meanz + gtem2 + facz * eta22[, jiter]  # poor aux var.
	}
	devfn=t((t(seed1))%*%t(LBF[,1:length(seed1)]))
	
xdata=x_mat[,(inidiscard+1):(T+inidiscard)]
zdata1=z_mat1[,(inidiscard+1):(T+inidiscard)]	
zdata2=z_mat2[,(inidiscard+1):(T+inidiscard)]
eta3=eta3[,(inidiscard+1):(T+inidiscard)]

x_mat=xdata
z_mat1=zdata1
z_mat2=zdata2

## Exogeneity-based benchmark: Z_t = X_t
zdataX=x_mat

for (nobs in TSET)
{
kappas <- c(20, 10, 5, 0)

devfn_vec <- as.numeric(devfn)

x_centered_for_y <- x_mat[, seq_len(nobs), drop = FALSE] -
  matrix(meanx, nrow = length(meanx), ncol = nobs)

xscore <- as.numeric(crossprod(devfn_vec * w_int, x_centered_for_y))

YY <- meany +
  outer(xscore, kappas / sqrt(nobs)) +
  matrix(eta3[1, seq_len(nobs)], nrow = nobs, ncol = length(kappas))

for (jjk in 1:ncol(YY))
{
ydata <- YY[seq_len(nobs), jjk]
ydata_c <- ydata - mean(ydata)

Z1n <- zdata1[, seq_len(nobs), drop = FALSE]
Z2n <- zdata2[, seq_len(nobs), drop = FALSE]
ZXn <- zdataX[, seq_len(nobs), drop = FALSE]

Z1n_c <- Z1n - rowMeans(Z1n)
Z2n_c <- Z2n - rowMeans(Z2n)
ZXn_c <- ZXn - rowMeans(ZXn)

ldata  <- sweep(Z1n_c, 2, ydata_c, "*")
ldata2 <- sweep(Z2n_c, 2, ydata_c, "*")
ldataX <- sweep(ZXn_c, 2, ydata_c, "*")
	
W <- WLIST[[as.character(nobs)]]

stat1 <- calc_stats(ldata,  W, w_int, nobs)
stat2 <- calc_stats(ldata2, W, w_int, nobs)
statX <- calc_stats(ldataX, W, w_int, nobs)

teststat  <- stat1["inf"]
teststat2 <- stat2["inf"]
teststatX <- statX["inf"]

teststatA  <- stat1[paste0("p", PSET[1])]
teststatB  <- stat1[paste0("p", PSET[2])]
teststatC  <- stat1[paste0("p", PSET[3])]
teststatD  <- stat1[paste0("p", PSET[4])]

teststatA2 <- stat2[paste0("p", PSET[1])]
teststatB2 <- stat2[paste0("p", PSET[2])]
teststatC2 <- stat2[paste0("p", PSET[3])]
teststatD2 <- stat2[paste0("p", PSET[4])]

teststatAX <- statX[paste0("p", PSET[1])]
teststatBX <- statX[paste0("p", PSET[2])]
teststatCX <- statX[paste0("p", PSET[3])]
teststatDX <- statX[paste0("p", PSET[4])]
	
	xx_mat=ldata - rowMeans(ldata)
	hh2=t(LBF[2:(nt),])%*%xx_mat[2:(nt),]*(t[2]-t[1])
    xcoef=t(hh2)
	EVs=eigen(crossprod(xcoef)/nobs)$vectors[,1:pcsforbw]
	bandinput=xcoef%*%EVs
	if (ceil==1){band1 =  min(1+ceiling(bandfn(bandinput,kernf)),nobs-2)}else{band1 =min(1+round(bandfn(bandinput,kernf)),nobs-2)}	
	lrx1 = lrvar(xcoef, aband = 1, bband = band1 ,kernel = kernf)
    eigen1 = eigen(lrx1$omega / nobs)$values 	
	BANDCOL1=append(BANDCOL1,band1)
	
	xx_mat=ldata2 - rowMeans(ldata2)
    hh2=t(LBF[2:(nt),])%*%xx_mat[2:(nt),]*(t[2]-t[1])
    xcoef=t(hh2)
	EVs=eigen(crossprod(xcoef)/nobs)$vectors[,1:pcsforbw]
	bandinput=xcoef%*%EVs
	if (ceil==1){band1 =  min(1+ceiling(bandfn(bandinput,kernf)),nobs-2)}else{band1 =min(1+round(bandfn(bandinput,kernf)),nobs-2)}
	lrx1 = lrvar(xcoef, aband = 1, bband = band1 ,kernel = kernf)
    eigen2 = eigen(lrx1$omega / nobs)$values 	
	BANDCOL2=append(BANDCOL2,band1)
	
   ## Exogeneity-based benchmark: Z_t = X_t
	xx_mat=ldataX - rowMeans(ldataX)
	hh2=t(LBF[2:(nt),])%*%xx_mat[2:(nt),]*(t[2]-t[1])
	xcoef=t(hh2)
	EVs=eigen(crossprod(xcoef)/nobs)$vectors[,1:pcsforbw]
	bandinput=xcoef%*%EVs
	if (ceil==1){band1 =  min(1+ceiling(bandfn(bandinput,kernf)),nobs-2)}else{band1 =min(1+round(bandfn(bandinput,kernf)),nobs-2)}
	lrx1 = lrvar(xcoef, aband = 1, bband = band1 ,kernel = kernf)
	eigenX = eigen(lrx1$omega / nobs)$values 	
	BANDCOLX=append(BANDCOLX,band1)
	
  	if (ceil == 1) {
	band2 <- 5 + ceiling(nobs^bw2)
	} else {
	band2 <- 5 + round(nobs^bw2)
	}

	cvs <- calc_cv3(eigen1, eigen2, eigenX, band2, cviter, alpha = 0.95)

	cv1 <- cvs["cv1"]
	cv2 <- cvs["cv2"]
	cvX <- cvs["cvX"]
	  
	  if(nobs==TSET[1]){tindex=1};  if(nobs==TSET[2]){tindex=2}; if(nobs==TSET[3]){tindex=3}; 
	  TPOWER[iteration,jjk,tindex]=((teststat>cv1)) ;   
	TPOWER2[iteration,jjk,tindex]=((teststat2>cv2))
	TPOWERX[iteration,jjk,tindex]=((teststatX>cvX))

	TPOWERA[iteration,jjk,tindex]=((teststatA>cv1));  
	TPOWERA2[iteration,jjk,tindex]=((teststatA2>cv2))
	TPOWERAX[iteration,jjk,tindex]=((teststatAX>cvX))

	TPOWERB[iteration,jjk,tindex]=((teststatB>cv1));  
	TPOWERB2[iteration,jjk,tindex]=((teststatB2>cv2))	
	TPOWERBX[iteration,jjk,tindex]=((teststatBX>cvX))

	TPOWERC[iteration,jjk,tindex]=((teststatC>cv1));  
	TPOWERC2[iteration,jjk,tindex]=((teststatC2>cv2))
	TPOWERCX[iteration,jjk,tindex]=((teststatCX>cvX))

	TPOWERD[iteration,jjk,tindex]=((teststatD>cv1));  
	TPOWERD2[iteration,jjk,tindex]=((teststatD2>cv2))
	TPOWERDX[iteration,jjk,tindex]=((teststatDX>cvX))

	TSTAT[iteration,jjk,tindex]=teststat;   
	TSTAT2[iteration,jjk,tindex]=teststat2
	TSTATX[iteration,jjk,tindex]=teststatX

	TSTATA[iteration,jjk,tindex]=teststatA;   
	TSTATA2[iteration,jjk,tindex]=teststatA2  
	TSTATAX[iteration,jjk,tindex]=teststatAX

	TSTATB[iteration,jjk,tindex]=teststatB;   
	TSTATB2[iteration,jjk,tindex]=teststatB2
	TSTATBX[iteration,jjk,tindex]=teststatBX

	TSTATC[iteration,jjk,tindex]=teststatC;   
	TSTATC2[iteration,jjk,tindex]=teststatC2
	TSTATCX[iteration,jjk,tindex]=teststatCX

	TSTATD[iteration,jjk,tindex]=teststatD;   
	TSTATD2[iteration,jjk,tindex]=teststatD2
	TSTATDX[iteration,jjk,tindex]=teststatDX
	}
}
	   if (iteration%%60==0) {print(paste("kernf=",kernf," bw=",bw2," PSET=",PSET[1],PSET[2],PSET[3],PSET[4]," FAR=",FAR," ranbasis=",ranbasis," ranini=",ranini," FARmean=",FARmean1," FARsd=",FARvar1," upermean=",upermean," upersd=",upersd," ranfirst=",ranfirst," limlim=",limlim," limlim2=",limlim2," cviter=",cviter," intercept",sep=""))}
	   if (iteration%%15==0) {print(round(c(colMeans(TPOWER[1:iteration,,1]),colMeans(TPOWERA[1:iteration,,1]),colMeans(TPOWERB[1:iteration,,1]),colMeans(TPOWERC[1:iteration,,1]),colMeans(TPOWERD[1:iteration,,1])),digits=2))}
	   if (iteration%%15==0) {print(round(c(colMeans(TPOWER2[1:iteration,,1]),colMeans(TPOWERA2[1:iteration,,1]),colMeans(TPOWERB2[1:iteration,,1]),colMeans(TPOWERC2[1:iteration,,1]),colMeans(TPOWERD2[1:iteration,,1])),digits=2))}
       if (iteration%%15==0) {print(round(c(colMeans(TPOWERX[1:iteration,,1]),colMeans(TPOWERAX[1:iteration,,1]),colMeans(TPOWERBX[1:iteration,,1]),colMeans(TPOWERCX[1:iteration,,1]),colMeans(TPOWERDX[1:iteration,,1])), digits=2))}   
	  if (iteration%%15==0) {print(round(c(colMeans(TPOWER[1:iteration,,1]),colMeans(TPOWER[1:iteration,,2]),colMeans(TPOWER[1:iteration,,3])),digits=2))}
	  if (iteration%%15==0) {print(round(c(colMeans(TPOWER2[1:iteration,,1]),colMeans(TPOWER2[1:iteration,,2]),colMeans(TPOWER2[1:iteration,,3])),digits=2))}

}

# code used to save the results. 
#save(file=paste("results/ARtest","kernf",kernf,"bw",bw2,"PSET",PSET[1],PSET[2],PSET[3],PSET[4],"upermean",upermean,"upersd",upersd,"ranbasis",ranbasis,"ranini",ranini,"ceil",ceil,"pcsforbw",pcsforbw,"ranfirst",ranfirst,"limlim",limlim,
 #               "limlim2",limlim2,"cviter",cviter,"facz",facz, "_sizepower_interceptver1.RData",sep=""),TPOWER,TPOWERA,TPOWERB,TPOWERC,TPOWERD, TSTAT,TSTATA,TSTATB,TSTATC,TSTATD, TPOWER2,TPOWERA2,TPOWERB2,TPOWERC2,TPOWERD2,
#     TSTAT2,TSTATA2,TSTATB2,TSTATC2,TSTATD2,TPOWERX,TPOWERAX,TPOWERBX,TPOWERCX,TPOWERDX, TSTATX,TSTATAX,TSTATBX,TSTATCX,TSTATDX, BANDCOL1,BANDCOL2,BANDCOLX)
}
}


############################################################
## Table 4: proposed-test rejection rates.
##
## TPOWER, TPOWERA, TPOWERB, TPOWERC, and TPOWERD correspond to the informative design with p = infinity, 7, 3, 1, and 0, respectively. TPOWER2, TPOWERA2, TPOWERB2, TPOWERC2, and TPOWERD2 are the corresponding objects for the weakly informative design.
##
## Since kappas = c(20,10,5,0), columns are reordered as c(4,3,2,1) to match the table order kappa = 0, 5, 10, 20.
##
## To reproduce the beta_u=0.1 and beta_u=0.25 blocks in Table 2, run the code separately with upermean = 0.1 and 0.25. To reproduce the Bartlett and Parzen blocks, run it with kernf = 7 and 1, respectively.
############################################################

lapply(seq_along(TSET), function(i) {
  list(
    T = TSET[i],
    informative = round(100 * rbind(
      "p=inf" = colMeans(TPOWER[,  c(4,3,2,1), i], na.rm = TRUE),
      "p=7"   = colMeans(TPOWERA[, c(4,3,2,1), i], na.rm = TRUE),
      "p=3"   = colMeans(TPOWERB[, c(4,3,2,1), i], na.rm = TRUE),
      "p=1"   = colMeans(TPOWERC[, c(4,3,2,1), i], na.rm = TRUE),
      "p=0"   = colMeans(TPOWERD[, c(4,3,2,1), i], na.rm = TRUE)
    ), 1),
    weakly_informative = round(100 * rbind(
      "p=inf" = colMeans(TPOWER2[,  c(4,3,2,1), i], na.rm = TRUE),
      "p=7"   = colMeans(TPOWERA2[, c(4,3,2,1), i], na.rm = TRUE),
      "p=3"   = colMeans(TPOWERB2[, c(4,3,2,1), i], na.rm = TRUE),
      "p=1"   = colMeans(TPOWERC2[, c(4,3,2,1), i], na.rm = TRUE),
      "p=0"   = colMeans(TPOWERD2[, c(4,3,2,1), i], na.rm = TRUE)
    ), 1)
  )
})

