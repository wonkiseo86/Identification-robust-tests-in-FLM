#This code is made from moderate modification of the publicly available code of Seong and Seo (2025) for the simulation experiments mentioned in the paper.

library(fda);library(AER);library(Matrix);library(xtable);library(sde);library(mvtnorm)
set.seed(123456)

library(sde)

#set WD
setwd("Set folder containing all the files")


for (goodIV in c(0))
{
useiid=1
setaval=0.5

covtoeigen = dget('covtoeigennew2.R')
innersingle = dget('inprod.R')
bandfn = dget("bandfn.R")
lrvar = dget("lr_var_v2.R")

##Evaluation points##
kgrid = 50
kpoint = seq(from = 0, to = 1, length.out = kgrid)


## Fourier basis 
lbnumber2=50
nt = 50 ; t = kpoint; 

LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(innersingle(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(innersingle(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}
LBF=cbind(rep(1,length(t)),LBF)
lb=LBF
for(i in 2:lbnumber2){  
  for(j in 1:i)  { 
    if (j != i) {lb[,i] = lb[,i]-(innersingle(lb[,i],lb[,j],t)/innersingle(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber2){
  LBF[,i] = lb[,i]/(sqrt(innersingle(lb[,i],lb[,i],t)))
} 

## functions
inner = function (f, g, grid) {
  h = c(f)*g
  return(  colSums((0.5 * h[1:(length(grid) - 1),] + 0.5 * h[2:(length(grid)),]) * 
                     (grid[2] - grid[1]))     )
}

norm = function (f, grid) {
  return(innersingle(f,f,grid))
}

operator2 = function(f,grid){
  temp = NULL
  for (ss in grid){
    temp = append(temp,innersingle(1-abs(ss-grid)^2,f,grid))
  }
  return(temp)
}

load("trcz.RData")
trc_z = mean(trace_z)


nbb = 51
regfac = 2
rsq = .5
titer = 2000
nobslist = c(100,200,400)
mobs = max(nobslist) 


dcutlist = c(.8,.85,.9,.95)

ndlist = length(dcutlist)
 
##Functional Observation Name##
fzname = rep("", mobs)
fxname = rep("", mobs)
fyname = rep("", mobs)
for (i in 1:mobs) {
  fzname[i] <- paste("fz", i, sep = "")
  fxname[i] <- paste("fx", i, sep = "")
  fyname[i] <- paste("fy", i, sep = "")
}

## Partitioned support to draw the confidence band ##
nu_mat  = NULL
fig_cri = 5
for(i in 0:(kgrid/fig_cri -1)){
  nu_mat = cbind(nu_mat,1*( .1*i  <kpoint & kpoint<= (.1*(i+1))) )
} 

## result matrices ##
lnobs = length(nobslist)
npcurves= 25
nsetaq = 3

for (seta in c(setaval)) {
  set.seed(12345)
  npcurve2 = npcurves * 2 + 1
  test_cfive = array(NA, dim = c(titer, npcurve2, lnobs))
  test_cftsls = array(NA, dim = c(titer, npcurve2, lnobs))
  test_cftik = array(NA, dim = c(titer, npcurve2, lnobs))

  ngrid = 50
  gpoint = seq(from = 0, to = 1, length.out = ngrid)
  nnbasis1 = 51
  nnbasis2 = 51
  basis_fn1 = create.fourier.basis(rangeval = c(0, 1), nbasis = nnbasis1)
  basis_fn2 = create.fourier.basis(rangeval = c(0, 1), nbasis = nnbasis2)

  seta1 = 1e-1
  sig_etaexp = seta * c(0.9^(0:(nbb - 1)))
  sig_etasps = seta * c(.99^(0:(nsetaq - 1)), seta1^(1:(nbb - nsetaq)))

  sig_eta = sqrt(sum(sig_etaexp^2) / sum(sig_etasps^2)) * sig_etasps
  trc_eta = sum(sig_eta^2)

  bval = sqrt((rsq / (1 - rsq)) * (1 / (6 * (trc_z + trc_eta))))

  cval = .8

  for (iter in 1:titer) {
    rn_coef = rnorm(11, 0, 0.5^(0:10))
    test_tmp = eval.basis(kpoint, create.fourier.basis(rangeval = c(0, 1), nbasis = 11)) %*% rn_coef
    true_curve = operator2(test_tmp, kpoint)
    gam_seq = (rnorm(11, 0, 0.5^(0:10)))^2
    c_gam_val = seq(0.01, 0.5, length.out = npcurves) / sum(gam_seq)

    gam_seq = apply(matrix(c_gam_val, ncol = 1), 1, function(x) { return(sqrt(x * gam_seq)) })
    gam_seq = cbind(-gam_seq[, c((dim(gam_seq)[2]):1)], 0, gam_seq)
    pert_true_curve = eval.basis(kpoint, create.fourier.basis(rangeval = c(0, 1), nbasis = 11)) %*% gam_seq

    power_curve = true_curve + pert_true_curve

    dz = matrix(NA, nrow = ngrid, ncol = mobs)
    dzz = matrix(NA, nrow = ngrid, ncol = mobs)
    dx = matrix(NA, nrow = ngrid, ncol = mobs)
    dy = matrix(NA, nrow = ngrid, ncol = mobs)
    for (j in 1:mobs) {

      alpha = 3 * runif(1) + 2
      beta = 3 * runif(1) + 2
 ############################################################################################################################################   
      eta1 = BBridge(0, 0, 0, 1, ngrid - 1)
      eta2 = BBridge(0, 0, 0, 1, ngrid - 1)
      eta3 = BBridge(0, 0, 0, 1, ngrid - 1)
      v = cval * eta2 + sqrt(1 - cval^2) * eta1
      dz[, j] = dbeta(gpoint, alpha, beta) + eval.basis(gpoint, create.fourier.basis(rangeval = c(0, 1), nbasis = nbb)) %*% rnorm(nbb, 0, sig_eta)
      # dx[,j] = innersingle(dz[,j], LBF[,4], gpoint) + innersingle(dz[,j], LBF[,5], gpoint) + innersingle(dz[,j], LBF[,6], gpoint) + innersingle(dz[,j], LBF[,7], gpoint) + eta2
      dx[, j] = bval[1] * dz[, j] + eta2
      dzz[, j] = innersingle(dz[, j], LBF[,2], gpoint) * LBF[,2] + eta3 #+ innersingle(dz[, j], LBF[,3], gpoint) * LBF[,3] + innersingle(dz[, j], LBF[,1], gpoint) * LBF[,1] 
      dy[, j] = operator2(dx[, j], gpoint) + v
    }
   ################################################################################################################################################ 
    if(goodIV==1){dz=dz}
	if(goodIV==0){dz = dzz}
    if(goodIV==2){dz = dx}
   ################################################################################################################################################ 
   
    colnames(dx[, ]) = fxname
    colnames(dy[, ]) = fyname
    colnames(dz) = fzname

    fz = Data2fd(gpoint, dz, basis_fn2)
    fz_m = fz$coefs

    for (model in 1:length(nobslist)) {

      nobs = nobslist[model]
      fx = Data2fd(gpoint, dx[, 1:nobs], basis_fn1)
      fy = Data2fd(gpoint, dy[, 1:nobs], basis_fn1)
      fz = fd(fz_m[, 1:nobs], basis_fn2)

      fx = center.fd(fx)
      fz = center.fd(fz)
      fy = center.fd(fy)

      z_fd = eval.fd(kpoint, fz)
      y_fd = eval.fd(kpoint, fy)
      x_fd = eval.fd(kpoint, fx)

      covzx = covtoeigen(fdobj1 = fx, fdobj2 = fz, nharm1 = nnbasis1, nharm2 = nnbasis2)
      covzy = covtoeigen(fdobj1 = fy, fdobj2 = fz, nharm1 = nnbasis1, nharm2 = nnbasis2)
      covz = pca.fd(fz, nharm = nnbasis2)

      covx = pca.fd(fx, nharm = nnbasis1)

      reg_alpha = regfac * sum(covzx$values^2) / sqrt(nobs)
      reg_alpha_1 = regfac * sum(covz$values^2) / sqrt(nobs)
      reg_alpha_x = regfac * sum(covx$values^2) / sqrt(nobs)

      kzz = max(c(1, which((covz$values^2) > reg_alpha_1)))
      kxx = max(c(1, which((covx$values^2) > reg_alpha_x)))

      psivec = c((covz$values^(-1/2))[1:kzz], rep(0, nnbasis2 - kzz))
      wmatz = apply(eval.fd(gpoint, covz$harmonics), 2, inner, eval.fd(gpoint, fz), gpoint) %*% diag(psivec)

      psicoef = covz$harmonics$coefs
      ztildecoef = wmatz %*% t(psicoef)
      fztilde = fd(t(ztildecoef), covz$harmonics$basis)

      covzxtilde = covtoeigen(fdobj1 = fx, fdobj2 = fztilde, nharm1 = nnbasis1, nharm2 = nnbasis2)
      covzytilde = covtoeigen(fdobj1 = fy, fdobj2 = fztilde, nharm1 = nnbasis1, nharm2 = nnbasis2)

      reg_alpha_2 = regfac * sum(covzxtilde$values^2) / sqrt(nobs)

      kzx = max(c(1, which((covzx$values^2) > reg_alpha)))
      kzxtilde = min(max(c(1, which((covzxtilde$values^2) > reg_alpha_2)), kzz))

      x_lammat = diag((1 / covzx$values))[1:kzx, 1:kzx]

      phi_fd = eval.fd(kpoint, covzx$harmonics2[1:min(nnbasis1, nnbasis2)])
      xi_fd = eval.fd(kpoint, covzx$harmonics1[1:nnbasis1])
      z_phi_inmat = apply(as.matrix(phi_fd[, 1:kzx]), 2, inner, z_fd, kpoint)
      wmat = y_fd %*% z_phi_inmat %*% x_lammat / nobs

      x_tmp = wmat %*% t(apply(as.matrix(xi_fd[, 1:kzx]), 2, inner, x_fd, kpoint))
      u_hat = y_fd - x_tmp

      test_function = function(test_f, true_test1, res, grid) {
        nobs_tmp = dim(res)[2]
        chatuu = crossprod(inner(test_f, res, grid)) / nobs_tmp
        test_num = (z_fd) %*% inner(test_f, y_fd, grid) / nobs_tmp ### C_{yz}psi

        test_num1 = (z_fd) %*% inner(true_test1, x_fd, grid) / nobs_tmp

        Tstat = (nobs_tmp / (chatuu)) * norm(test_num - test_num1, kpoint)
        return(Tstat)
      }

      test_function2 = function(test_f, true_test1, res, grid) {
        nobs_tmp = dim(res)[2]
        zdata = z_fd - rowMeans(z_fd)
        xdata = x_fd - rowMeans(x_fd)
        ydata = inner(test_f, y_fd, grid)
        ydata = ydata - mean(ydata)
        ldata = zdata
        for (i in 1:nobs_tmp) {
          ldata[, i] = (ydata[i]) * (zdata[, i]) - innersingle(xdata[, i], true_test1, kpoint) * (zdata[, i])
        }
        lm1 = rowMeans(ldata)
        teststat = innersingle(lm1, lm1, kpoint) * (nobs_tmp)

        xx_mat = ldata - rowMeans(ldata)
        hh2 = t(LBF[2:(nt), ]) %*% xx_mat[2:(nt), ] * (t[2] - t[1])
        xcoef = t(hh2)

        if (useiid == 0) {
          lrx1 = crossprod(xcoef)
          eigen1 = eigen(lrx1 / nobs_tmp)$values
          EVs = eigen(crossprod(xcoef) / nobs_tmp)$vectors[, 1:5]
          bandinput = xcoef %*% EVs
          band1 = min(1 + ceiling(bandfn(bandinput, 1)), nobs_tmp - 2)
          lrx1 = lrvar(xcoef, aband = 1, bband = band1, kernel = 1)
          eigen1 = eigen( (0.5*lrx1$omega+0.5*t(lrx1$omega))/ nobs_tmp)$values
        }

        if (useiid == 1) {
          lrx1 = crossprod(xcoef)
          eigen1 = eigen(lrx1 / nobs_tmp)$values
        }

        band2 = 5 + ceiling(nobs_tmp^(0.333))
        CV1 = NULL
        for (jj in 1:1000) {
          atem = rnorm(band2, 0, 1)^2
          cvlong1 = eigen1[1:band2] * atem
          cv1 = sum(cvlong1[1:band2])
          CV1 = append(CV1, cv1)
        }
        cv1 = quantile(CV1, 0.95)
		(teststat > cv1)
        return(teststat > cv1)
      }

      litermax = 500
      dseq = min(ceiling(nobs^(1/3)), nnbasis2)
      mattmp = diag(sqrt(sqrt((covz$values[1:dseq])^2))) %*% matrix(rnorm(dseq * litermax), nrow = dseq, ncol = litermax)
      emp_critical = c(sort(diag(crossprod(mattmp)))[(litermax * 0.95)])

      test_cfive[iter, , model] = apply(power_curve, 2, function(x) { return(test_function(test_tmp, x, u_hat, kpoint) > emp_critical) })
      test_cftsls[iter, , model] = apply(power_curve, 2, function(x) { return(test_function2(test_tmp, x, u_hat, kpoint)) })
    }
    if (iter %% 500 == 0) {
      print(round(c(iter, apply(test_cfive[1:iter, (npcurves - 1):(npcurves + 3), (lnobs - 1):lnobs], 3, colMeans)), digits = 3))
      print(round(c(iter, apply(test_cftsls[1:iter, (npcurves - 1):(npcurves + 3), (lnobs - 1):lnobs], 3, colMeans)), digits = 3))
    }
  }
 
# code used to save the results at the local (Not used)   
 # save(file=paste("results/comparison/ARtest","TSET",nobslist[1],nobslist[2],nobslist[3],"seta",seta,"useiid",useiid,"goodIV",goodIV,"_comparison_symmetrized.RData",sep=""), test_cfive, test_cftsls)

}


 RESULT=(round(c(iter, apply(test_cfive[1:iter, c((npcurves + 1), (npcurves + 6), (npcurves + 11), (npcurves + 16)), 1:3], 3, colMeans)), digits = 3))
 RESULTK=(round(c(iter, apply(test_cftsls[1:iter, c((npcurves + 1), (npcurves + 6), (npcurves + 11), (npcurves + 16)), 1:3], 3, colMeans)), digits = 3))
 rbind(RESULT*100)
rbind(RESULTK*100)		

		options(digits=2)
		 TABLEA=xtable(rbind(RESULTK*100), row.names = FALSE, col.names = FALSE, quote = FALSE,digits=1)
		 TABLEB=xtable(rbind(RESULT*100), row.names = FALSE, col.names = FALSE, quote = FALSE,digits=1) 
		 TABLEA
		 TABLEB
		 
		 }
		 
		 
############################################################
## Report results by sample size
############################################################

sel_curve <- c(
  npcurves + 1,
  npcurves + 6,
  npcurves + 11,
  npcurves + 16
)

make_RESULT <- function(A, iter, sel_curve, nobslist, digits = 1) {
  out <- apply(
    A[seq_len(iter), sel_curve, seq_along(nobslist), drop = FALSE],
    c(2, 3),
    mean,
    na.rm = TRUE
  )

  ## selected curves x sample-size -> sample-size x selected curves
  out <- t(out)

  out <- 100 * out

  rownames(out) <- paste0("T=", nobslist)
  colnames(out) <- paste0("curve", c(1, 6, 11, 16))

  ## numeric version, useful for further calculation
  out_num <- round(out, digits = digits)

  ## display version, keeps trailing zeros
  out_fmt <- matrix(
    formatC(out_num, format = "f", digits = digits),
    nrow = nrow(out_num),
    ncol = ncol(out_num),
    dimnames = dimnames(out_num)
  )

  return(list(num = out_num, fmt = out_fmt))
}


############################################################
## How to check entries in Table~\ref{tab3}
##
## The variable goodIV determines the design:
##     goodIV = 1: informative design with Z_t,
##     goodIV = 0: weakly informative design with Z_t^circ.
##
## Run this code separately with goodIV = 1 and goodIV = 0 to obtain the two design blocks reported in Table~\ref{tab3}.
##
## The selected columns
##     sel_curve = c(npcurves+1, npcurves+6, npcurves+11, npcurves+16)
## correspond to
##     magni^2 = 0, 0.05, 0.10, 0.15.
##
## RESULTK_fmt reports the proposed test based on test_cftsls.
## RESULT_fmt  reports the functional IV test based on test_cfive.
##
## Rows correspond to sample sizes T = 100, 200, 400.
############################################################

RES1 <- make_RESULT(test_cfive,  iter, sel_curve, nobslist, digits = 1)
RES2 <- make_RESULT(test_cftsls, iter, sel_curve, nobslist, digits = 1)

RESULT  <- RES1$num
RESULTK <- RES2$num

RESULT_fmt  <- RES1$fmt
RESULTK_fmt <- RES2$fmt

RESULTK_fmt  # Proposed test
RESULT_fmt   # Functional IV test

TABLEA <- xtable(
  data.frame(T = nobslist, RESULTK_fmt, check.names = FALSE),
  caption = "Rejection rates for the proposed test",
  label = "tab:resultK",
  align = c("l", "c", rep("c", ncol(RESULTK_fmt)))
)

TABLEB <- xtable(
  data.frame(T = nobslist, RESULT_fmt, check.names = FALSE),
  caption = "Rejection rates for the comparison test",
  label = "tab:result",
  align = c("l", "c", rep("c", ncol(RESULT_fmt)))
)

print(TABLEA, include.rownames = FALSE, sanitize.text.function = identity)
print(TABLEB, include.rownames = FALSE, sanitize.text.function = identity)


