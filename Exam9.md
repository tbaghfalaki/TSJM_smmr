### R code for simulation study to compare the proposed MMJM approach with block-diagonal covariance matrix for the random effects

```
rm(list = ls())
#setwd("/Users/taban/Desktop/apr/Feb8")
# setwd("/Users/taban/Desktop/D2025")
library(Matrix)
library(MASS);library(R2jags);library(doParallel);library(foreach)
library(timeROC);library(survival)
library(statmod);library(PermAlgo) # Permutation algorithm to generate survival times dependent on time-varying covariates
library(mvtnorm);library(pec);library(dplyr);library(DPCri)
source("DP4b.R")
cl <- makeCluster(10)
registerDoParallel(cl)
nsujet <- 1000
AUC_All=BS_All=matrix(0,4,5)
inde1=(nsujet/2+1):nsujet
Limp=100
S=c(0,0.25,0.5,0.75)
# real values of parameters
NN=100
resultsss <-foreach(ij=1:NN,.packages=c("MASS","survival","DPCri","dplyr","R2jags","statmod","mvtnorm","PermAlgo","Matrix")) %dopar% {
  # Y1 (continuous)
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
  
  Sigma=as.matrix(Sigma)
  
  Sigma=as.matrix(t(Sigma))

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
                          eventRandom = round(rexp(nsujet, 0.006) + 1, 0), # ~40% death
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
  surv.data=survDat
  long.data=longDat
  # Number of patients and number of longitudinal observations per patient
  n <- length(surv.data$Id)
  M <- table(long.data$Id)
  # Survival and censoring times
  st <- surv.data$deathTimes
  
  Rate=table(surv.data$Event)/n
  #boxplot(M~surv.data$Event)
  MM11=mean(M[surv.data$Event==1]);MM21=median(M[surv.data$Event==1])
  MM10=mean(M[surv.data$Event==0]);MM20=median(M[surv.data$Event==0])
  MM1=c(mean(M),MM11,MM10)
  MM2=c(median(M),MM21,MM20)
  names(MM1)=names(MM2)=c("all","cause1","censored")
  # Longitudinal information in matrix format
  time = Y1=Y2=Y3=Y4 <- X1 <- X2 <- matrix(NA, n, max(M))
  for (i in 1:n) {
    time[i, 1:M[i]] <- long.data$time[long.data$Id == i]
    Y1[i, 1:M[i]] <- long.data$Y1[long.data$Id == i]
    Y2[i, 1:M[i]] <- long.data$Y2[long.data$Id == i]
    Y3[i, 1:M[i]] <- long.data$Y3[long.data$Id == i]
    Y4[i, 1:M[i]] <- long.data$Y4[long.data$Id == i]
    
    X1[i, 1:M[i]] <- long.data$X1[long.data$Id == i]
    X2[i, 1:M[i]] <- long.data$X2[long.data$Id == i]
  }
  
  X1 <- X1[, 1]
  X2 <- X2[, 1]
  
  W <- cbind(1,surv.data$x1 ,surv.data$w1) # Fixed effects
  
  X <- array(1, dim = c(n, max(M), 4)) # Fixed effects
  X[, , 2] <- time
  X[, , 3] <- x1
  X[, , 4] <- x2
  
  Z <- array(1, dim = c(n, max(M), 2)) # Random effects
  Z[, , 2] <- time
  
  ########  BUGS code  ########
  peice=quantile(st,seq(.2,0.8,length=4))
  delta=rep(0,n)
  for(i in 1:n){
    if(st[i]<=peice[1])(delta[i]=1)
    if(st[i]>peice[1] & st[i]<= peice[2])(delta[i]=2)
    if(st[i]>peice[2] & st[i]<= peice[3])(delta[i]=3)
    if(st[i]>peice[3] & st[i]< peice[4])(delta[i]=4)
    if(st[i]>=peice[4])(delta[i]=5)
  }
  Delta=nnet::class.ind(delta)
  table(delta)
  
  INDTRAIN <- 1:(nsujet/2)
  INDVALID <- surv.data$Id[-INDTRAIN]
  
  long.data_t <- subset(
    long.data,
    long.data$Id %in% INDTRAIN
  )
  surv.data_t <- subset(
    surv.data,
    surv.data$Id %in% INDTRAIN
  )
  
  long.data_v <- subset(
    long.data,
    long.data$Id %in% INDVALID
  )
  surv.data_v <- subset(
    surv.data,
    surv.data$Id %in% INDVALID
  )
  
  
  
  ############################# MMJM   ############################# 
  start2 <- Sys.time()
  i.jags <- function() {
    list(gamma1 = rnorm(1), gamma2 = rnorm(1),gamma3 = rnorm(1),gamma4 = rnorm(1), 
         betaL1 = rnorm(dim(X)[3]),betaL2 = rnorm(dim(X)[3]),betaL3 = rnorm(dim(X)[3]),
         tau1 = 1,tau2 = 1,tau3 = 1, Omega1 = diag(runif(2)), Omega2 = diag(runif(2)),
         Omega3 = diag(runif(2)), Omega4 = diag(runif(2)))
  }
  d.jags <- list(n = nsujet/2, M = M, Time = st, Y1 = Y1, Y2=Y2, Y3=Y3, Y4=Y4,  
                 XL1 =X, ZL1 = Z,XL2 =X, ZL2 = Z,XL3 =X, ZL3 = Z,XL4 =X, ZL4 = Z,
                 death = death, mub = rep(0, 2), V = diag(1, 2),  zeros = rep(0,n),
                 NbetasL = dim(X)[3],x1=x1,x2=x2,delta=Delta,s=peice,J=length(peice)+1)
  
  parameters <- c("gamma1", "gamma2", "gamma3","gamma4", "betaL1","betaL2","betaL3", "betaL4", 
                  "sigma1", "sigma2","sigma3", "sigma4", "Sigma1","Sigma2","Sigma3","Sigma4", "h")
  
  model.file="jm4bd.R"
  
  #file.show(model.file)
  sim123 <- R2jags::jags(data=d.jags, inits=i.jags, parameters,n.chains = 1,
                         n.iter=10000, model.file=model.file)
  sim123
  
  end2 <- Sys.time()
  TimeMulti=difftime(end2,start2,units ="mins")
  
  Sigma_m=bdiag(sim123$BUGSoutput$mean$Sigma1,sim123$BUGSoutput$mean$Sigma2,
                sim123$BUGSoutput$mean$Sigma3,sim123$BUGSoutput$mean$Sigma4)
  betaL1_m=sim123$BUGSoutput$mean$betaL1
  betaL2_m=sim123$BUGSoutput$mean$betaL2
  betaL3_m=sim123$BUGSoutput$mean$betaL3
  betaL4_m=sim123$BUGSoutput$mean$betaL4
  
  gamma_m=c(sim123$BUGSoutput$mean$gamma1,sim123$BUGSoutput$mean$gamma2,sim123$BUGSoutput$mean$gamma3,
            sim123$BUGSoutput$mean$gamma4)
  sigma_m=c(sim123$BUGSoutput$mean$sigma1,sim123$BUGSoutput$mean$sigma2,sim123$BUGSoutput$mean$sigma3,
            sim123$BUGSoutput$mean$sigma4)
  h_m=sim123$BUGSoutput$mean$h
  
  Est_mul=gamma_m
  
  Est_mul_l=sim123$BUGSoutput$summary[c(34:37),3]
  Est_mul_u=sim123$BUGSoutput$summary[c(34:37),7]
  
  CRM=rep(0,length(gamma_m))
  for(kk in 1:length(gamma_m)){
    if(gammareal[kk]>Est_mul_l[kk] & gammareal[kk]< Est_mul_u[kk])(CRM[kk]=1)
  }
  
  ###### DP ####
  CRI_MMJM=matrix(0,length(S),2)
  for(kks in 1:length(S)){ 
    s=S[kks];Dt=S[2]
    Sigma_m=as.matrix(bdiag(sim123$BUGSoutput$mean$Sigma1,sim123$BUGSoutput$mean$Sigma2,
                            sim123$BUGSoutput$mean$Sigma3,sim123$BUGSoutput$mean$Sigma4),8,8)
    
    betaL1_m=sim123$BUGSoutput$mean$betaL1
    betaL2_m=sim123$BUGSoutput$mean$betaL2
    betaL3_m=sim123$BUGSoutput$mean$betaL3
    betaL4_m=sim123$BUGSoutput$mean$betaL4
    
    gamma_m=c(sim123$BUGSoutput$mean$gamma1,sim123$BUGSoutput$mean$gamma2,sim123$BUGSoutput$mean$gamma3,sim123$BUGSoutput$mean$gamma4)
    sigma_m=c(sim123$BUGSoutput$mean$sigma1,sim123$BUGSoutput$mean$sigma2,sim123$BUGSoutput$mean$sigma3,sim123$BUGSoutput$mean$sigma4)
    h_m=sim123$BUGSoutput$mean$h
    
    ###########################   ###########################    ###########################   ###########################
    
    i.jags <- function() {
      list(b=matrix(0,n/2,2*nmark))
    }
    
    M1=M
    for(ii in 1:n){
      M1[ii]=length(Y1[ii,][t<=s][!is.na(Y1[ii,][t<=s])])
    } 
    M1[M1==0]=1 
    
    d.jags <- list(n = nsujet/2, M = M1[(nsujet/2+1):nsujet], Time = rep(s,nsujet/2), 
                   Y1 = Y1[(nsujet/2+1):nsujet,], Y2=Y2[(nsujet/2+1):nsujet,], 
                   Y3=Y3[(nsujet/2+1):nsujet,], Y4=Y4[(nsujet/2+1):nsujet,], 
                   XL1 =X[(nsujet/2+1):nsujet,,], ZL1 = Z[(nsujet/2+1):nsujet,,],
                   XL2 =X[(nsujet/2+1):nsujet,,], ZL2 = Z[(nsujet/2+1):nsujet,,],XL3 =X[(nsujet/2+1):nsujet,,],
                   ZL3 = Z[(nsujet/2+1):nsujet,,],
                   XL4 =X[(nsujet/2+1):nsujet,,],
                   ZL4 = Z[(nsujet/2+1):nsujet,,],
                   death = death[(nsujet/2+1):nsujet], mub = rep(0, 2*nmark),  Nb = 2*nmark, zeros = rep(0,n/2),
                   x1=x1[(nsujet/2+1):nsujet],x2=x2[(nsujet/2+1):nsujet],delta=Delta[(nsujet/2+1):nsujet,],s=peice,
                   betaL1=betaL1_m,
                   betaL2=betaL2_m,
                   betaL3=betaL3_m,
                   betaL4=betaL4_m,
                   gamma1=gamma_m[1],gamma2=gamma_m[2],gamma3=gamma_m[3],gamma4=gamma_m[4],
                   sigma1=sigma_m[1],sigma2=sigma_m[2],sigma3=sigma_m[3],sigma4=sigma_m[4],
                   Sigma=Sigma_m,h=h_m)
    parameters <- c("b")
    model.file="jm4b.R"
    
    #file.show(model.file)
    sim123b <- R2jags::jags(data=d.jags, inits=i.jags, parameters,n.chains = 1,
                            n.iter=5000, model.file=model.file)
    bhatm=sim123b$BUGSoutput$mean$b
    
    DPtotal=c()
    for(i in inde1){ 
      DPtotal[i]=DP4b(s=s,Dt=Dt,betaL1=betaL1_m,betaL2=betaL2_m,betaL3=betaL3_m,betaL4=betaL4_m,
                      gamma=gamma_m,sigma=sigma_m,h=h_m,Sigma=Sigma_m,peice=peice,
                      y1=Y1[i,],y2=Y2[i,],y3=Y3[i,],y4=Y4[i,],x1=x1[i],x2=x2[i],b=bhatm[i-(nsujet/2),])
    }
    
    DPtotal=DPtotal[inde1]
    #########################
    
    CRI_MMJM[kks,]=Criteria(s = s, t = Dt, Survt = surv.data_v$deathTimes, 
                            CR =surv.data_v$Event, P = 1-DPtotal, cause = 1)$Cri[,1]
    
  }
  
  MMJM=list(Est_mul=Est_mul,CRM=CRM,CRI_MMJM=CRI_MMJM)
  
  
  list(MMJM=MMJM)
  
  
}

stopCluster(cl)

save(resultsss, file = "BD3n1.RData")

```


##### JAGS code
```
model{
  for(i in 1:n){
    #Longitudinalobservations
    for(j in 1:M[i]){
      Y1[i,j]~dnorm(mu1[i,j],tau1)
      mu1[i,j]<-inprod(betaL1[],XL1[i,j,])+inprod(b1[i,1:2],ZL1[i,j,])
      Y2[i,j]~dnorm(mu2[i,j],tau2)
      mu2[i,j]<-inprod(betaL2[],XL2[i,j,])+inprod(b2[i,1:2],ZL2[i,j,])
      Y3[i,j]~dnorm(mu3[i,j],tau3)
      mu3[i,j]<-inprod(betaL3[],XL3[i,j,])+inprod(b3[i,1:2],ZL3[i,j,])
      Y4[i,j]~dnorm(mu4[i,j],tau4)
      mu4[i,j]<-inprod(betaL4[],XL4[i,j,])+inprod(b4[i,1:2],ZL4[i,j,])
      
    }
    #Survival and censoring times
    #Hazard function
    
    
    Alpha0[i]<- gamma1*(betaL1[1]+betaL1[3]*x1[i]+betaL1[4]*x2[i]+b1[i,1])+
      gamma2*(betaL2[1]+betaL2[3]*x1[i]+betaL2[4]*x2[i]+b2[i,1])+
      gamma3*(betaL3[1]+betaL3[3]*x1[i]+betaL3[4]*x2[i]+b3[i,1])+
      gamma4*(betaL4[1]+betaL4[3]*x1[i]+betaL4[4]*x2[i]+b4[i,1])
    
    Alpha1[i]<- gamma1*(betaL1[2]+b1[i,2])+gamma2*(betaL2[2]+b2[i,2])+gamma3*(betaL3[2]+b3[i,2])+gamma4*(betaL4[2]+b4[i,2])
    
    
    
    haz[i]<- (h[1]*equals(delta[i,1],1)+h[2]*equals(delta[i,2],1)+h[3]*equals(delta[i,3],1)+
                h[4]*equals(delta[i,4],1)+h[5]*equals(delta[i,5],1))*exp(Alpha0[i]+Alpha1[i]*Time[i])
    #Cumulative hazard function 
    chaz1[i]<-h[1]*Time[i]*equals(delta[i,1],1)+
      (h[1]*s[1]+h[2]*(Time[i]-s[1]))*equals(delta[i,2],1)+
      (h[1]*s[1]+h[2]*(s[2]-s[1])+h[3]*(Time[i]-s[2]))*equals(delta[i,3],1)+
      (h[1]*s[1]+h[2]*(s[2]-s[1])+h[3]*(s[3]-s[2])+h[4]*(Time[i]-s[3]))*equals(delta[i,4],1)+
      (h[1]*s[1]+h[2]*(s[2]-s[1])+h[3]*(s[3]-s[2])+h[4]*(s[4]-s[3])+h[5]*(Time[i]-s[4]))*equals(delta[i,5],1)
    chaz2[i]<- h[1]*(exp(Alpha1[i]*Time[i])-1)*equals(delta[i,1],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*Time[i])-exp(Alpha1[i]*s[1])))*equals(delta[i,2],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*Time[i])-exp(Alpha1[i]*s[2])))*equals(delta[i,3],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*s[3])-exp(Alpha1[i]*s[2]))+h[4]*(exp(Alpha1[i]*Time[i])-exp(Alpha1[i]*s[3])))*equals(delta[i,4],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*s[3])-exp(Alpha1[i]*s[2]))+h[4]*(exp(Alpha1[i]*s[4])-exp(Alpha1[i]*s[3]))+h[5]*(exp(Alpha1[i]*Time[i])-exp(Alpha1[i]*s[4])))*equals(delta[i,5],1)
    
    
    
    
    
    chaz[i]<-exp(Alpha0[i])*chaz2[i]/Alpha1[i]
    
    #Log-survival function log(S)=-H(t) 
    logSurv[i]<-  -chaz[i]
    #Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-death[i]*log(haz[i])-logSurv[i]
    zeros[i]~dpois(phi[i])
    #Random effects
    b1[i,1:2]~dmnorm(mub[],Omega1[,])
    b2[i,1:2]~dmnorm(mub[],Omega2[,])
    b3[i,1:2]~dmnorm(mub[],Omega3[,])
    b4[i,1:2]~dmnorm(mub[],Omega4[,])
    
  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL1[l]~dnorm(0,0.001)
    betaL2[l]~dnorm(0,0.001)
    betaL3[l]~dnorm(0,0.001)
    betaL4[l]~dnorm(0,0.001)
    
  }
  
  for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }
  
  
  gamma1~dnorm(0,0.001)
  gamma2~dnorm(0,0.001)
  gamma3~dnorm(0,0.001)
  gamma4~dnorm(0,0.001)
  
  alpha~dgamma(0.01,0.01)
  tau1~dgamma(0.01,0.01)
  tau2~dgamma(0.01,0.01)
  tau3~dgamma(0.01,0.01)
  tau4~dgamma(0.01,0.01)
  
  #Derive dquantity
  sigma1<-1/tau1
  sigma2<-1/tau2
  sigma3<-1/tau3
  sigma4<-1/tau4
  
  Omega1[1:2,1:2]~dwish(V[,],2)
  Sigma1[1:2,1:2]<-inverse(Omega1[,])
  
  Omega2[1:2,1:2]~dwish(V[,],2)
  Sigma2[1:2,1:2]<-inverse(Omega2[,])
  
  Omega3[1:2,1:2]~dwish(V[,],2)
  Sigma3[1:2,1:2]<-inverse(Omega3[,])
  
  Omega4[1:2,1:2]~dwish(V[,],2)
  Sigma4[1:2,1:2]<-inverse(Omega4[,])
}


```



