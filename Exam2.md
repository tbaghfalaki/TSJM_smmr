
### The Moreno-Betancur et al. (2018) approach


#### Generating the data as described in https://github.com/tbaghfalaki/TSJM_smmr/blob/main/Exam1.md.

```
rm(list = ls())
library(mvtnorm)
library(R2jags)
library(doParallel)
library(foreach)
library(timeROC)
library(survival)
library(statmod)
library(PermAlgo)
library(survtd)

nsujet <- 1000
AUC_All=BS_All=matrix(0,4,5)
inde1=(nsujet/2+1):nsujet


# real values of parameters
N=100
resultsss1=list()
for(k1 in 1:N) {
  set.seed(k1)

  BetaLreal= Beta1 <- c(-0.5, 0.5, 0.5, 0.5)
  Beta2 <- c(-0.5, 0.5, 0.5, 0.5)
  Beta3 <- c(-0.5, 0.5, 0.5, 0.5)
  Beta4 <- c(-0.5, 0.5, 0.5, 0.5)
  alpha <- 1
  sigma <- sqrt(0.5) # residual error
  # Survival
  gammareal=gamma <- c(-1, -1, 1, 1)*0.2
  gapLongi <- 0.2 # gap between longi measurements
  gap <- 0.02 # used to generate a lot of time points because the permutation
  # algorithm chooses among those time points to define survival times
  followup <- 2 # follow-up time
  mestime <- seq(0, followup, gap) # measurement times
  t=timesLongi <- mestime[which(round(mestime - round(mestime / gapLongi, 0) * gapLongi, 6) == 0)] # visit times
  time <- rep(mestime, nsujet) # time column
  nmesindiv <- followup / gap + 1 # max. number of individual measurements
  nmesy <- nmesindiv * nsujet # max. total number of longi measurements
  # the number is reduced because of censoring by the terminal event
  idY <- rep(1:nsujet, each = nmesindiv) # individual id
  
  # random effects variance and covariance matrix
  rho <- 0.5
  rho1 <- 0.4
  Sigmar=Sigma <- Matrix::bdiag(
    matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2),
    matrix(c(1, rho, rho, 1), 2, 2)
  )
  Sigma[3:4,1:2] <- 0.1
  Sigma[5:6,1:2] <- 0.5
  Sigma[7:8,1:2] <- 0.1
  
  Sigma[5:6,3:4] <- 0.6
  Sigma[7:8,3:4] <- 0.4
  
  Sigma[7:8,5:6] <- 0.5
  Sigma=as.matrix(Sigma)
  
  Sigma=as.matrix(t(Sigma))
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  
  solve(Sigma)
  
  nmark <- 4
  
  Sigmar=Sigma_r <- Sigma <- matrix(Sigma, 2 * nmark, 2 * nmark)
  Realpara <- c(
    Beta1, Beta1, Beta1, Beta1, gamma,
    rep((sigma)^2, nmark), Sigma_r[1, ], Sigma_r[2, 2:8], Sigma_r[3, 3:8],
    Sigma_r[4, 4:8], Sigma_r[5, 5:8], Sigma_r[6, 6:8], Sigma_r[7, 7:8], Sigma_r[8, 8]
  )
  ############### real
  u <- rmvnorm(nsujet, rep(0, 2 * nmark), Sigma)
  b1_int <- rep(u[, 1], each = nmesindiv) # random intercept Y1
  b1_slo <- rep(u[, 2], each = nmesindiv) # random slope Y1
  b2_int <- rep(u[, 3], each = nmesindiv) # random intercept Y2
  b2_slo <- rep(u[, 4], each = nmesindiv) # random slope Y2
  b3_int <- rep(u[, 5], each = nmesindiv) # random intercept Y3
  b3_slo <- rep(u[, 6], each = nmesindiv) # random intercept Y3
  b4_int <- rep(u[, 7], each = nmesindiv) # random intercept Y3
  b4_slo <- rep(u[, 8], each = nmesindiv) # random intercept Y3
  
  
  x1 <- rnorm(nsujet, 1, 0.5)
  X1 <- rep(x1, each = nmesindiv) # continuous covariate
  x2 <- rbinom(nsujet, 1, 0.5)
  X2 <- rep(x2, each = nmesindiv) # binary covariate
  
  # linear predictors
  linPredY1 <- (Beta1[1] + b1_int) + (Beta1[2] + b1_slo) * time + Beta1[3] * X1 + Beta1[4] * X2
  linPredY2 <- (Beta2[1] + b2_int) + (Beta2[2] + b2_slo) * time + Beta2[3] * X1 + Beta2[4] * X2
  linPredY3 <- (Beta3[1] + b3_int) + (Beta3[2] + b3_slo) * time + Beta3[3] * X1 + Beta3[4] * X2
  linPredY4 <- (Beta4[1] + b4_int) + (Beta4[2] + b4_slo) * time + Beta4[3] * X1 + Beta4[4] * X2
  
  # continuous outcomes
  Y1 <- rnorm(nmesy, linPredY1, sigma)
  Y2 <- rnorm(nmesy, linPredY2, sigma)
  Y3 <- rnorm(nmesy, linPredY3, sigma)
  Y4 <- rnorm(nmesy, linPredY4, sigma)
  # Permutation algorithm to generate survival times that depends on the linear predictors
  DatTmp <- permalgorithm(nsujet, nmesindiv,
                          Xmat = matrix(c(linPredY1, linPredY2, linPredY3, linPredY4), nrow = nsujet * nmesindiv),
                          eventRandom = round(rexp(nsujet, 0.02) + 1, 0), # ~40% death
                          censorRandom = runif(nsujet, 1, nmesindiv), # uniform random censoring
                          XmatNames = c("linPredY1", "linPredY2", "linPredY3", "linPredY4"), # association
                          betas = c(gamma)
  ) # association parameters
  
  # extract last line for each Id (= death/censoring time)
  DatTmp2 <- DatTmp[c(which(diff(DatTmp[, "Id"]) == 1), dim(DatTmp)[1]), c("Id", "Event", "Stop")]
  DatTmp2$deathTimes <- mestime[DatTmp2$Stop + 1] # deathtimes
  survDat <- DatTmp2[, c("Id", "deathTimes", "Event")]
  DatTmp$time <- mestime[DatTmp$Start + 1] # measurements times of the biomarker
  DatTmp$Uid <- paste(DatTmp$Id, DatTmp$time) # unique identifier to match covariates and observed biomarker values
  longDat3 <- merge(DatTmp[, c("Uid", "Id", "time")], cbind("Uid" = paste(idY, time), X1, X2, Y1, Y2, Y3, Y4), by = c("Uid"))
  longDat <- sapply(longDat3[longDat3$time %in% timesLongi, -1], as.numeric)
  longDat <- as.data.frame(longDat[order(longDat[, "Id"], longDat[, "time"]), ])
  summary(survDat) # survival dataset
  summary(longDat) # longitudinal dataset
  head(longDat)
  names(survDat)
  names(longDat)
  median(table(longDat$Id))
  min(table(longDat$Id))
  max(table(longDat$Id))
  sum(survDat$Event) / length(survDat$Event)
  death=survDat$Event
  
  
  ### ######################
  DataAll=merge(longDat,survDat,by=c("Id"))
  names(DataAll)
  head(DataAll)
```



#### Implementing the Moreno-Betancur et al. (2018) approach

```
  start2 <- Sys.time()
  
  DataAll1=subset(DataAll,
                  DataAll$Id %in% c(1:(nsujet/2)))
  
  AAA=c()
  AAA=try({survtd(Surv(time=deathTimes,event=Event)~td(Y1)+td(Y2)+td(Y3)+td(Y4),
                  data=DataAll1,  id="Id", visit.time="time", model="Cox",
                  method="MIJM", M=5, G=20,time.trend=as.formula("~x"))
  }, 
  silent = TRUE)
  
  
  Est_MB=try({c(rep(0,16),AAA$Results[,1], rep(0,40)) }, 
             silent = TRUE)
  
  
  Est_MB_l=try({c(rep(0,16),
                  AAA$Results[,3],
                  rep(0,40)) }, 
               silent = TRUE)
  
  
  
  Est_MB_u=try({c(rep(0,16),
                  AAA$Results[,4],
                  rep(0,40)) }, 
               silent = TRUE)
  
  
  CRMB=try({rep(0,length(Realpara)) }, 
           silent = TRUE)
  
  
  sd_MB=try({c(rep(0,16),
               AAA$Results[,2],
               rep(0,40))   }, 
            silent = TRUE)
  
  
  try({for(kk in 1:length(Realpara)){
    if(Realpara[kk]>Est_MB_l[kk] & Realpara[kk]< Est_MB_u[kk])(CRMB[kk]=1)
  } }, 
  silent = TRUE)
  
  end2 <- Sys.time()
  TIME_MB=difftime(end2,start2,units ="mins")
  
  CRMB
  TIME_MB
  #######################################
  resultsss1[[k1]]=list(
    Est_MB=Est_MB,Est_MB_l=Est_MB_l,
    Est_MB_u=Est_MB_u,CRMB=CRMB,sd_MB=sd_MB,TIME_MB=TIME_MB)
  
  print(k1)
  
}
