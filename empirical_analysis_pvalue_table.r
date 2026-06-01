
############################################################
## Empirical p-value table generator
## Keeping the original variable names and code structure as much as possible
##
## This version loops over:
##   kernf    = 7 (Bartlett), 1 (Parzen)
##   variden  = 2, 3, 4, 1, 5, 0 (CLR, LHR, LRHR, LCDF, PDF, QF)
##   maxmoment= 1, 2, 3, 4
##
## It stores the final p-values in TABLE_PVAL and prints the LaTeX table
## in the format used in the paper.
############################################################


############################################################
## Loading packages
############################################################

library(fda)
library(tseries)
library(sandwich)
library(sde)
library(variables)
library(basefun)
library(polynom)
library(fracdiff)
library(LongMemoryTS)
library(arfima)
library(truncnorm)


############################################################
## Loading functions
############################################################

setwd("Set folder containing all the files")

inner  = dget("inprod.R")
operator = dget("operator.R")
lrvar = dget("lr_var_v2.R")
bandfn = dget("bandfn.R")

############################################################
## Read MATLAB data directly
############################################################

############################################################
## Helper functions for reading MATLAB data
############################################################

library(R.matlab)

get_mat_object <- function(mat, target_name) {
  nm <- names(mat)

  clean <- function(x) {
    tolower(gsub("[^A-Za-z0-9]", "", x))
  }

  idx <- which(clean(nm) == clean(target_name))

  if (length(idx) == 0) {
    cat("Available objects:\n")
    print(nm)
    print(sapply(mat, dim))
    stop("Cannot find object corresponding to: ", target_name)
  }

  mat[[idx[1]]]
}

movmean_matlab <- function(x, k = c(6, 6)) {
  ## MATLAB movmean(x, [6 6]) with endpoint shrinkage.
  x <- as.numeric(x)
  n <- length(x)
  out <- numeric(n)

  for (i in seq_len(n)) {
    idx1 <- max(1, i - k[1])
    idx2 <- min(n, i + k[2])
    out[i] <- mean(x[idx1:idx2])
  }

  out
}

make_month_dummies <- function(n) {
  SEAS <- matrix(0, nrow = n, ncol = 12)
  month_index <- ((seq_len(n) - 1) %% 12) + 1
  SEAS[cbind(seq_len(n), month_index)] <- 1
  SEAS
}


############################################################
## p-value approximation
############################################################

weighted_chisq_pvalues_mc <- function(teststat, eigenvalues, band2, B = 500000,
                                      chunk = 100000, seed = 1,
                                      progress_every = 100000) {
  set.seed(seed)

  lambda <- Re(eigenvalues[seq_len(band2)])
  lambda <- pmax(lambda, 0)

  teststat <- as.numeric(teststat)
  exceed <- numeric(length(teststat))

  done <- 0
  start_time <- Sys.time()
  next_print <- progress_every

  while (done < B) {
    m <- min(chunk, B - done)

    z2 <- matrix(rnorm(m * band2)^2, nrow = m, ncol = band2)
    cvdraw <- as.numeric(z2 %*% lambda)

    exceed <- exceed + vapply(teststat, function(q) {
      sum(cvdraw >= q)
    }, numeric(1))

    done <- done + m

    if (done >= next_print || done == B) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      rate <- done / elapsed
      remaining <- (B - done) / rate

      cat(sprintf(
        "MC progress: %d / %d draws completed (%.1f%%). Elapsed: %.1f sec. ETA: %.1f sec.\n",
        done, B, 100 * done / B, elapsed, remaining
      ))

      next_print <- next_print + progress_every
    }
  }

  100 * exceed / B
}


############################################################
## Settings
############################################################

## 0: quantile function, 1: logit CDF, 2: CLR,
## 3: log hazard, 4: log reverse hazard, 5: original pdf, 6: cdf
variden_SET <- c(2, 3, 4, 1, 5, 0)
variden_LABEL <- c("CLR", "LHR", "LRHR", "LCDF", "PDF", "QF")
names(variden_LABEL) <- as.character(variden_SET)

maxmoment_SET <- c(1, 2, 3, 4)

## sumindi = 1: Z_t = X_{t-1}
## sumindi = 2: Z_t = X_{t-1} + 0.5 X_{t-2}
## sumindi = 3: Z_t = X_{t-1} + 0.5 X_{t-2} + 0.25 X_{t-3}
## sumindi = 0: Z_t = X_t
sumindi_SET <- c(1, 2, 3, 0)

sumindi_LABEL <- c(
  "1" = "Case $Z_t = X_{t-1}$",
  "2" = "Case $Z_t = \\sum_{j=1}^{\\ell} (0.5)^{j-1} X_{t-j}$ with $\\ell=2$",
  "3" = "Case $Z_t = \\sum_{j=1}^{\\ell} (0.5)^{j-1} X_{t-j}$ with $\\ell=3$",
  "0" = "Case $Z_t = X_t$ (no correction for endogeneity)"
)

## kernf = 7: Bartlett, kernf = 1: Parzen
kernf_SET <- c(7, 1)
kernf_LABEL <- c("7" = "Bartlett", "1" = "Parzen")

## Monte Carlo settings
cvinter <- 500000 #5000000
chunk_size <- 100000

## Tuning parameters
ceil <- 1
maxbasis <- 50
bw2 <- 0.333
ranbasis <- 4
ranini <- 1
PSET <- c(7, 3, 1, 0)
pcsforbw <- 5
FARmean1 <- 0.4
FARvar1 <- 0.4
ranfirst <- 1
scalran <- 0.2
limlim <- 0.8
limlim2 <- -0.2
upersd <- -limlim2
cviter <- 500000
FAR <- 1
scalevar <- 1
ulower <- limlim2
uupper <- limlim
intrange <- 3
intrange2 <- 3
FARmean2 <- 0.6
FARvar2 <- 0.6
flooring <- 10^(-7)
minlim <- limlim2
maxlim <- limlim


############################################################
## Final result arrays
############################################################

TABLE_PVAL <- array(
  NA_real_,
  dim = c(length(sumindi_SET), length(maxmoment_SET), length(variden_SET), length(kernf_SET)),
  dimnames = list(
    sumindi = as.character(sumindi_SET),
    maxmoment = as.character(maxmoment_SET),
    variden = variden_LABEL,
    kernf = kernf_LABEL[as.character(kernf_SET)]
  )
)

## This keeps p = 7 results too, if needed later.
TABLE_PVAL_A <- TABLE_PVAL


############################################################
## Functions retained from your code
############################################################

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
{
  aa=0
  for (j in 1:nn)
  {
    aa=aa+(coef_fn(p,j,nn)^2)*nn
  }
  return(sqrt(1/aa))
}


############################################################
## Main loop
############################################################

spec_counter <- 0
total_specs <- length(sumindi_SET) * length(kernf_SET) * length(variden_SET) * length(maxmoment_SET)

for (sumindi in sumindi_SET) {

  for (kernf in kernf_SET) {

    reportresult2 <- NULL
    reportresult2a <- NULL

    cat("\n============================================================\n")
    cat("Starting sumindi =", sumindi, ", kernf =", kernf, "(", kernf_LABEL[as.character(kernf)], ")\n")
    cat("============================================================\n")

    for (variden in variden_SET) {

      STATVAL1=NULL
      TESTVAL1=NULL
      PVAL1=NULL
      CVVAL1=NULL

      STATVAL2=NULL
      TESTVAL2=NULL
      PVAL2=NULL
      CVVAL2=NULL

      cat("\nTransform =", variden_LABEL[as.character(variden)], "\n")

      for (maxmoment in maxmoment_SET) {

        spec_counter <- spec_counter + 1

        cat(sprintf(
          "\n[%d/%d] sumindi=%s, kernf=%s, variden=%s, maxmoment=%s\n",
          spec_counter, total_specs, sumindi, kernf_LABEL[as.character(kernf)],
          variden_LABEL[as.character(variden)], maxmoment
        ))

## Temperature density: 300 x 601, rows are monthly PDFs
Temp_Density_mat <- readMat("KTemp_Density.mat")
Temp_Density <- get_mat_object(Temp_Density_mat, "Temp_Density")
X <- t(as.matrix(Temp_Density))   ## 601 x 300

## Temperature quantile function: 300 x 601, rows are monthly QFs
Temp_Quant_mat <- readMat("Temp_Quant_v202501.mat")
Temp_Quant <- get_mat_object(Temp_Quant_mat, "Temp_Quant")
X2 <- t(as.matrix(Temp_Quant))    ## 601 x 300

## Electricity demand data
Elec_mat <- readMat("Elec_demand_data.mat")

Pre_data <- get_mat_object(Elec_mat, "Pre_data")
R_weight <- get_mat_object(Elec_mat, "R_weight")

Pre_data <- as.matrix(Pre_data)

## MATLAB:
## D_Elec_demand = log((Pre_data(:,1)/10^6)./Pre_data(:,2));
D_Elec_demand <- log((Pre_data[, 1] / 10^6) / Pre_data[, 2])

## MATLAB:
## D_Elec_demand_MA = D_Elec_demand - movmean(D_Elec_demand, [6 6]);
D_Elec_demand_MA <- D_Elec_demand - movmean_matlab(D_Elec_demand, c(6, 6))

y <- as.vector(D_Elec_demand_MA)

## No need for sumstats_dens.csv
SEAS <- make_month_dummies(length(y))

        tt=seq(-20,40,length.out=601)

        meantem=NULL
        for (i in 1:ncol(X))
        {
          aatem= sum(tt*X[,i])*(tt[2]-tt[1])
          meantem=append(meantem,aatem)
        }

        vartem=NULL
        for (i in 1:ncol(X))
        {
          aatem= sum(((tt-meantem[i])^2)*X[,i])*(tt[2]-tt[1])
          vartem=append(vartem,aatem)
        }

        sdtem=sqrt(vartem)

        stdmoments = matrix(0,nrow=ncol(X),ncol=16)
        stdmoments[,1]=meantem
        stdmoments[,2]=sdtem

        for (jjj in 3:16)
        {
          for(i in 1:ncol(X))
          {
            aatem=sum(((tt-meantem[i])^jjj)*X[,i])*(tt[2]-tt[1])
            aatem=aatem/((vartem[i])^(jjj/2))
            stdmoments[i,jjj]= aatem
          }
        }

        if(variden==1)
        {
          for (jj in 1:ncol(X))
          {
            X[,jj] = cumsum(X[,jj])/sum(X[,jj])
            X[,jj] = log((X[,jj]+flooring)/(1-X[,jj]+flooring))
          }
        }

        if(variden==0)
        {
          X=X2
        }

        if(variden==2)
        {
          for (jj in 1:ncol(X))
          {
            X[,jj] = log((X[,jj]+flooring))-mean(log((X[,jj]+flooring)))
          }
        }

        if(variden==3)
        {
          for (jj in 1:ncol(X))
          {
            disttem= cumsum(X[,jj])/sum(X[,jj])
            X[,jj] = log((X[,jj]+flooring)/(1-disttem+flooring))
          }
        }

        if(variden==4)
        {
          for (jj in 1:ncol(X))
          {
            disttem= cumsum(X[,jj])/sum(X[,jj])
            X[,jj] = log((X[,jj]+flooring)/(disttem+flooring))
          }
        }

        if(variden==5)
        {
          X=X
        }

        if(variden==6)
        {
          for (jj in 1:ncol(X))
          {
            disttem= cumsum(X[,jj])/sum(X[,jj])
            X[,jj] = disttem
          }
        }

        TSET=length(y)

        CONST=matrix(0,nrow=length(TSET),ncol=length(PSET))
        for (ibb in 1:length(TSET))
        {
          for (jbb in 1:length(PSET))
          {
            CONST[ibb,jbb]=const_fn(PSET[jbb],TSET[ibb])
          }
        }

        lbnumber2=100
        nt = 601
        t = (0:(nt-1))/(nt-1)

        LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
        for (i in 1:(lbnumber2/2)){
          LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
          LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
        }

        LBF=cbind(rep(1,length(t)),LBF)
        lb=LBF

        for(i in 2:lbnumber2){
          for(j in 1:i)  {
            if (j != i) {
              lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]
            }
          }
        }

        for(i in 1:lbnumber2){
          LBF[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
        }

        #################################################
        #### Main empirical computation
        #################################################

        T=max(TSET)
        ydata = y[(sumindi+1):T]
        xdata=X[,(sumindi+1):ncol(X)]

        if(sumindi==0){
          zdata2=xdata
        }else{
          aatemtem=0
          for (jkj in 1:sumindi)
          {
            aatemtem=aatemtem + X[,jkj:(ncol(X)-(sumindi-jkj+1))]*2^{-((sumindi-jkj+1))} * 2
          }
          zdata2=aatemtem
        }

        Wfull=stdmoments[(sumindi+1):T,]

        if (maxmoment==0)
        {
          W=cbind(rep(1,nrow(Wfull)))
        }else{
          W=Wfull[,1:(maxmoment)]
        }

        ydata = lm(ydata~cbind(W,SEAS[(sumindi+1):T,]))$residuals

        xdatafourier=t(LBF)%*%xdata*(t[2]-t[1])
        zdata2fourier=t(LBF)%*%zdata2*(t[2]-t[1])

        xdatanew=xdatafourier
        zdata2new=zdata2fourier

        for (i in 1:nrow(xdatanew))
        {
          xdatanew[i,]=lm(xdatafourier[i,]~cbind(W,SEAS[(sumindi+1):T,]))$residuals
          zdata2new[i,]=lm(zdata2fourier[i,]~cbind(W,SEAS[(sumindi+1):T,]))$residuals
        }

        xdata_resid=LBF%*%xdatanew
        zdata2_resid=LBF%*%zdata2new

        xdata=xdata_resid
        zdata2=zdata2_resid

        nobs=ncol(xdata)

        ldata=zdata2
        ldata2=zdata2

        for (i in 1:length(ydata))
        {
          ldata2[,i]=(ydata[i]-mean(ydata))*(zdata2[,i]-rowMeans(zdata2))
        }

        lm2=rowMeans(ldata2)

        test_p1=NULL
        test_p2=NULL

        for (pp in PSET)
        {
          llm1=0
          llm2=0

          if (pp==PSET[1]){constfac=CONST[1,1]}
          if (pp==PSET[2]){constfac=CONST[1,2]}
          if (pp==PSET[3]){constfac=CONST[1,3]}
          if (pp==PSET[4]){constfac=CONST[1,4]}

          for(j in 1:ncol(ldata))
          {
            nn=ncol(ldata)

            if(pp==0){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==1){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==2){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==3){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==4){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==5){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==6){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==7){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==8){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==9){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
            if(pp==10){llm2 = llm2+ constfac*coef_fn(pp,j,nn)*ldata2[,j]}
          }

          if (pp==PSET[1]){teststatA2= inner(llm2,llm2,t)*(nobs)}
          if (pp==PSET[2]){teststatB2= inner(llm2,llm2,t)*(nobs)}
          if (pp==PSET[3]){teststatC2= inner(llm2,llm2,t)*(nobs)}
          if (pp==PSET[4]){teststatD2= inner(llm2,llm2,t)*(nobs)}
        }

        teststat2= inner(lm2,lm2,t)*(nobs)

        xx_mat=ldata2 - rowMeans(ldata2)
        hh2=t(LBF[2:(nt),])%*%xx_mat[2:(nt),]*(t[2]-t[1])
        xcoef=t(hh2)

        EVs=eigen(crossprod(xcoef)/nobs)$vectors[,1:pcsforbw]
        bandinput=xcoef%*%EVs

        if (ceil==1){
          band1 =  min(1+ceiling(bandfn(bandinput,kernf)),nobs-2)
        }else{
          band1 =min(1+round(bandfn(bandinput,kernf)),nobs-2)
        }

        lrx1 = lrvar(xcoef, aband = 1, bband = band1 ,kernel = kernf)
        eigen2 = eigen(lrx1$omega / nobs)$values

        if (ceil == 1) {
          band2 <- 5 + ceiling(nobs^(bw2))
        } else {
          band2 <- 5 + round(nobs^(bw2))
        }

        TESTSTAT2 <- c(teststat2,teststatA2,teststatB2,teststatC2,teststatD2)

        pval2 <- weighted_chisq_pvalues_mc(
          teststat = TESTSTAT2,
          eigenvalues = eigen2,
          band2 = band2,
          B = cvinter,
          chunk = chunk_size,
          seed = 1,
          progress_every = 100000
        )

        STATVAL2 <- rbind(STATVAL2, TESTSTAT2)
        TESTVAL2 <- rbind(TESTVAL2, pval2 < 5)
        PVAL2 <- rbind(PVAL2, pval2)
      }

      report2 <- PVAL2

      ## p = infinity
      reportresult2 <- rbind(reportresult2, report2[,1])

      ## p = 7, kept only for checking if needed
      reportresult2a <- rbind(reportresult2a, report2[,2])
    }

    ## reportresult2 has rows = transformations and columns = maxmoment.
    TABLE_PVAL[as.character(sumindi), , variden_LABEL[as.character(variden_SET)], kernf_LABEL[as.character(kernf)]] <- t(reportresult2)
    TABLE_PVAL_A[as.character(sumindi), , variden_LABEL[as.character(variden_SET)], kernf_LABEL[as.character(kernf)]] <- t(reportresult2a)

    cat("\nFinished sumindi =", sumindi, ", kernf =", kernf_LABEL[as.character(kernf)], "\n")
    print(round(t(reportresult2), 2))
  }
}


############################################################
## Print LaTeX table
############################################################

make_empirical_table <- function(TABLE_PVAL) {

  fmt2 <- function(x) {
    formatC(x, format = "f", digits = 2)
  }

  lines <- character(0)

  lines <- c(lines, "\\begin{table}[h!]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, "\\small")
  lines <- c(lines, "\\caption{Testing functional association of electricity demand and temperature distribution, p-values (\\%)}")
  lines <- c(lines, "\\label{tabemp1}")
  lines <- c(lines, "\\renewcommand*{\\arraystretch}{0.75}")
  lines <- c(lines, "\\setlength{\\aboverulesep}{0.1ex}")
  lines <- c(lines, "\\setlength{\\belowrulesep}{0.1ex}")
  lines <- c(lines, "\\setlength{\\cmidrulesep}{0.1ex}")
  lines <- c(lines, "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}c rrrrrr c rrrrrr@{}}")
  lines <- c(lines, "\\toprule")
  lines <- c(lines, " & \\multicolumn{6}{c}{Bartlett} & & \\multicolumn{6}{c}{Parzen} \\\\")
  lines <- c(lines, "\\cmidrule(lr){2-7} \\cmidrule(lr){9-14}")
  lines <- c(lines, "$m \\;\\backslash\\; X_t$ & CLR & LHR & LRHR & LCDF & PDF & QF & & CLR & LHR & LRHR & LCDF & PDF & QF \\\\")
  lines <- c(lines, "\\midrule")

  for (sumindi in sumindi_SET) {

    lines <- c(lines, "\\addlinespace[0.6ex]")
    lines <- c(lines, paste0("\\multicolumn{14}{c}{\\textit{", sumindi_LABEL[as.character(sumindi)], "}} \\\\"))
    lines <- c(lines, "\\addlinespace[0.3ex]")

    for (mm in maxmoment_SET) {

      vals_bart <- sapply(variden_LABEL[as.character(variden_SET)], function(vv) {
        fmt2(TABLE_PVAL[as.character(sumindi), as.character(mm), vv, "Bartlett"])
      })

      vals_parz <- sapply(variden_LABEL[as.character(variden_SET)], function(vv) {
        fmt2(TABLE_PVAL[as.character(sumindi), as.character(mm), vv, "Parzen"])
      })

      line <- paste(c(mm, vals_bart, "", vals_parz), collapse = " & ")
      line <- paste0(line, " \\\\")
      lines <- c(lines, line)
    }
  }

  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular*}")
  lines <- c(lines, "\\vspace{-0.3em}")
  lines <- c(lines, "\\begin{flushleft}")
  lines <- c(lines, "\\footnotesize")
  lines <- c(lines, paste0(
    "Notes: The table reports approximate $p$-values (\\%) for $H_0: \\theta = 0$ in model~\\eqref{eq_empiric}, ",
    "computed via ", format(cvinter, big.mark = ","), " Monte Carlo draws of ",
    "$\\sum_{j=1}^{d_T} \\hat\\lambda_{\\cont,j}\\nu_j^2$. ",
    "$m$ denotes the number of standardized moments included as scalar controls."
  ))
  lines <- c(lines, "\\end{flushleft}")
  lines <- c(lines, "\\end{table}")

  paste(lines, collapse = "\n")
}

LATEX_TABLE <- make_empirical_table(TABLE_PVAL)

cat(LATEX_TABLE)
#writeLines(LATEX_TABLE, "table_empirical_pvalues.tex")  # NOT USED

#save(file = "empirical_pvalues_all_sumindi.RData",TABLE_PVAL,TABLE_PVAL_A,sumindi_SET,kernf_SET,variden_SET,maxmoment_SET,cvinter) # for checking the results, not used. 


############################################################
## Plot temperature density functions and quantile functions
############################################################

Temp_Density_mat <- readMat("KTemp_Density.mat")
Temp_Density_plot <- get_mat_object(Temp_Density_mat, "Temp_Density")

Temp_Quant_mat <- readMat("Temp_Quant_v202501.mat")
Temp_Quant_plot <- get_mat_object(Temp_Quant_mat, "Temp_Quant")

X_pdf_plot <- t(as.matrix(Temp_Density_plot))
X_qf_plot  <- t(as.matrix(Temp_Quant_plot))

tt_plot <- seq(-20, 40, length.out = nrow(X_pdf_plot))
qq_plot <- seq(0, 1, length.out = nrow(X_qf_plot))

pdf_mean_plot <- rowMeans(X_pdf_plot, na.rm = TRUE)
qf_mean_plot  <- rowMeans(X_qf_plot, na.rm = TRUE)



############################################################
## Plot temperature density functions and quantile functions
## This produces the figure directly in the R plotting window.
############################################################

Temp_Density_mat <- readMat("KTemp_Density.mat")
Temp_Density_plot <- get_mat_object(Temp_Density_mat, "Temp_Density")

Temp_Quant_mat <- readMat("Temp_Quant_v202501.mat")
Temp_Quant_plot <- get_mat_object(Temp_Quant_mat, "Temp_Quant")

X_pdf_plot <- t(as.matrix(Temp_Density_plot))
X_qf_plot  <- t(as.matrix(Temp_Quant_plot))

tt_plot <- seq(-20, 40, length.out = nrow(X_pdf_plot))
qq_plot <- seq(0, 1, length.out = nrow(X_qf_plot))

pdf_mean_plot <- rowMeans(X_pdf_plot, na.rm = TRUE)
qf_mean_plot  <- rowMeans(X_qf_plot, na.rm = TRUE)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1.0))

matplot(
  tt_plot, X_pdf_plot,
  type = "l",
  lty = 1,
  col = "grey80",
  xlab = "Temperature",
  ylab = "Density",
  main = "Temperature density functions"
)
lines(tt_plot, pdf_mean_plot, lwd = 2)
box()

matplot(
  qq_plot, X_qf_plot,
  type = "l",
  lty = 1,
  col = "grey80",
  xlab = "Probability level",
  ylab = "Temperature",
  main = "Temperature quantile functions"
)
lines(qq_plot, qf_mean_plot, lwd = 2)
box()

par(old_par)




# code used to generate the pdf file of the reported figure. 
#pdf("figure_temperature_pdf_qf.pdf", width = 10, height = 4.5)
#par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1.0))
#matplot(tt_plot, X_pdf_plot, type = "l", lty = 1, col = "grey80", xlab = "Temperature", ylab = "Density", main = "Temperature density functions"); lines(tt_plot, pdf_mean_plot, lwd = 2);box()
#matplot(qq_plot, X_qf_plot,type = "l",lty = 1,col = "grey80",xlab = "Probability level", ylab = "Temperature", main = "Temperature quantile functions");lines(qq_plot, qf_mean_plot, lwd = 2);box()
#dev.off()


