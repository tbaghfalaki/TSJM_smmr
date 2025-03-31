### R code for simulation study to compare the proposed TSJM approach with the multi-marker joint modeling (MMJM) approach, multiple two-stage (MTS), and the two-stage approach

```
rm(list = ls())
# setwd("/Users/taban/Desktop/apr/Feb8") # A path for saving all R code in there

# Loading necessary R packages  
library(MASS)
library(R2jags)
library(doParallel)
library(foreach)
library(timeROC)
library(survival)
library(statmod)
library(PermAlgo)  
library(mvtnorm)
library(pec)
library(dplyr)
library(DPCri)
source("DP4b.R")
# Considering parallel computation for improved performance in simulation study
cl <- makeCluster(10)
registerDoParallel(cl)
nsujet <- 1000
AUC_All <- BS_All <- matrix(0, 4, 5)
inde1 <- (nsujet / 2 + 1):nsujet
Limp <- 100
S <- c(0, 0.25, 0.5, 0.75)
# real values of parameters
NN <- 30
resultsss <- foreach(ij = 1:NN, .packages = c("MASS", "survival", "DPCri", "dplyr", "R2jags", "statmod", "mvtnorm", "PermAlgo")) %dopar% {
```
#### Generating data for joint modeling of multivariate longitudinal markers and time-to-event data. For more information, refer to: https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/multivariate


Data are generated using multi-marker joint models with $K=4$ longitudinal markers and a survival model. To be able to compare the performance to that of MMJM and limit  computation time, the simulation scenario considered only four markers.
More specifically, the data generation models were as follows:


$Y_{ik}(t)=\eta_{ik}(t|\bm{\beta}_k,\bm{b}_{ik})+\epsilon_{ikt}$

$=\beta_{0k}+\beta_{1k}t+\beta_{2k}x_{1i}+\beta_{3k}x_{2i}+b_{0ki}+b_{1ki} t+\varepsilon_{ikt}$


where 
$\varepsilon_{ikt} \sim N(0, \sigma_k^2),$
$x_1$ 
is generated from a binary distribution with a success probability of 0.5, and 
$x_2$ 
is generated from a normal distribution with a mean of zero and a standard deviation of 0.5. Also, 
$\bm{b}_i=(b_{01i},b_{11i},...,b_{04i},b_{14i})^{\top} \sim N(\bm{0},\bm{\Sigma}),$
$i=1,...,n$ and $t=0,0.2,0.4,...,2$.

```
  BetaLreal <- Beta1 <- Beta2 <- Beta3 <- Beta4 <- c(-0.5, 0.5, 0.5, 0.5)
  alpha <- 1
  sigma <- sqrt(0.5) # residual error
  # Survival
  gammareal <- gamma <- c(-1, -1, 1, 1) * 0.2
  gapLongi <- 0.2 # gap between longi measurements
  gap <- 0.02 # used to generate a lot of time points because the permutation
  # algorithm chooses among those time points to define survival times
  followup <- 2 # follow-up time
  mestime <- seq(0, followup, gap) # measurement times
  t <- timesLongi <- mestime[which(round(mestime - round(mestime / gapLongi, 0) * gapLongi, 6) == 0)] # visit times
  time <- rep(mestime, nsujet) # time column
  nmesindiv <- followup / gap + 1 # max. number of individual measurements
  nmesy <- nmesindiv * nsujet # max. total number of longi measurements
  # the number is reduced because of censoring by the terminal event
  idY <- rep(1:nsujet, each = nmesindiv) # individual id

  # random effects variance and covariance matrix
  rho <- 0.5
  rho1 <- 0.4
  Sigmar <- Sigma <- Matrix::bdiag(
    matrix(c(1, rho, rho, 1), 2, 2), 
    matrix(c(1, rho, rho, 1), 2, 2), 
    matrix(c(1, rho, rho, 1), 2, 2),
    matrix(c(1, rho, rho, 1), 2, 2)
  )
  Sigma[3:4, 1:2] <- 0.1
  Sigma[5:6, 1:2] <- 0.5
  Sigma[7:8, 1:2] <- 0.1

  Sigma[5:6, 3:4] <- 0.6
  Sigma[7:8, 3:4] <- 0.4

  Sigma[7:8, 5:6] <- 0.5
  Sigma <- as.matrix(Sigma)

  Sigma <- as.matrix(t(Sigma))
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

  solve(Sigma)

  nmark <- 4

  Sigmar <- Sigma_r <- Sigma <- matrix(Sigma, 2 * nmark, 2 * nmark)
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
  death <- survDat$Event
  surv.data <- survDat
  long.data <- longDat
  # Number of patients and number of longitudinal observations per patient
  n <- length(surv.data$Id)
  M <- table(long.data$Id)
  # Survival and censoring times
  st <- surv.data$deathTimes

  Rate <- table(surv.data$Event) / n
  # boxplot(M~surv.data$Event)
  MM11 <- mean(M[surv.data$Event == 1])
  MM21 <- median(M[surv.data$Event == 1])
  MM10 <- mean(M[surv.data$Event == 0])
  MM20 <- median(M[surv.data$Event == 0])
  MM1 <- c(mean(M), MM11, MM10)
  MM2 <- c(median(M), MM21, MM20)
  names(MM1) <- names(MM2) <- c("all", "cause1", "censored")
  # Longitudinal information in matrix format
  time <- Y1 <- Y2 <- Y3 <- Y4 <- X1 <- X2 <- matrix(NA, n, max(M))
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

  W <- cbind(1, surv.data$x1, surv.data$w1) # Fixed effects

  X <- array(1, dim = c(n, max(M), 4)) # Fixed effects
  X[, , 2] <- time
  X[, , 3] <- x1
  X[, , 4] <- x2

  Z <- array(1, dim = c(n, max(M), 2)) # Random effects
  Z[, , 2] <- time

  ########  BUGS code  ########
  peice <- quantile(st, seq(.2, 0.8, length = 4))
  delta <- rep(0, n)
  for (i in 1:n) {
    if (st[i] <= peice[1]) (delta[i] <- 1)
    if (st[i] > peice[1] & st[i] <= peice[2]) (delta[i] <- 2)
    if (st[i] > peice[2] & st[i] <= peice[3]) (delta[i] <- 3)
    if (st[i] > peice[3] & st[i] < peice[4]) (delta[i] <- 4)
    if (st[i] >= peice[4]) (delta[i] <- 5)
  }
  Delta <- nnet::class.ind(delta)
  table(delta)
```
##### Splitting data into training and validation sets for model evaluation

```
  INDTRAIN <- 1:(nsujet / 2)
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
```
#### The "TRUE" model (gold-standard) estimates association parameters using both the true values of the parameters and the random effects from the mixed model.

```
  mureal1 <- mureal2 <- mureal3 <- mureal4 <- rep(0, length(long.data_t$Id))

  for (i in 1:length(long.data_t$Id)) {
    mureal1[i] <- (BetaLreal[1] + u[long.data$Id[i], 1]) + (BetaLreal[2] + u[long.data$Id[i], 2]) * long.data$time[i] +
      BetaLreal[3] * long.data$X1[i] + BetaLreal[4] * long.data$X2[i]

    mureal2[i] <- (BetaLreal[1] + u[long.data$Id[i], 3]) + (BetaLreal[2] + u[long.data$Id[i], 4]) * long.data$time[i] +
      BetaLreal[3] * long.data$X1[i] + BetaLreal[4] * long.data$X2[i]

    mureal3[i] <- (BetaLreal[1] + u[long.data$Id[i], 5]) + (BetaLreal[2] + u[long.data$Id[i], 6]) * long.data$time[i] +
      BetaLreal[3] * long.data$X1[i] + BetaLreal[4] * long.data$X2[i]

    mureal4[i] <- (BetaLreal[1] + u[long.data$Id[i], 7]) + (BetaLreal[2] + u[long.data$Id[i], 8]) * long.data$time[i] +
      BetaLreal[3] * long.data$X1[i] + BetaLreal[4] * long.data$X2[i]
  }


  temp <- subset(surv.data_t, select = c(Id, deathTimes, Event)) # baseline data
  pbc21 <- tmerge(temp, temp, id = Id, endpt = event(deathTimes, Event))
  long.data1 <- tmerge(pbc21, long.data_t,
    id = Id, lY1 = tdc(time, mureal1),
    lY2 = tdc(time, mureal2), lY3 = tdc(time, mureal3), lY4 = tdc(time, mureal4)
  )
  head(long.data)

  step2_gold <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
    data = long.data1, id = Id, cluster = Id, ties = "breslow"
  )

  Est_gold_l <- step2_gold$coefficients - qnorm(0.975) * sqrt(diag(step2_gold$var))
  Est_gold_u <- step2_gold$coefficients + qnorm(0.975) * sqrt(diag(step2_gold$var))

  CRgold <- rep(0, length(gammareal))
  for (kk in 1:length(gammareal)) {
    if (gammareal[kk] > Est_gold_l[kk] & gammareal[kk] < Est_gold_u[kk]) (CRgold[kk] <- 1)
  }
```
#### Computing the "TRUE" dynamic prediction based on the gold-standard model
```
  # Calculate linear predictors for the validation set
  long.data_v <- long.data_v %>%
    mutate(
      lY1 = (BetaLreal[1] + u[Id, 1]) + (BetaLreal[2] + u[Id, 2]) * time + BetaLreal[3] * X1 + BetaLreal[4] * X2,
      lY2 = (BetaLreal[1] + u[Id, 3]) + (BetaLreal[2] + u[Id, 4]) * time + BetaLreal[3] * X1 + BetaLreal[4] * X2,
      lY3 = (BetaLreal[1] + u[Id, 5]) + (BetaLreal[2] + u[Id, 6]) * time + BetaLreal[3] * X1 + BetaLreal[4] * X2,
      lY4 = (BetaLreal[1] + u[Id, 7]) + (BetaLreal[2] + u[Id, 8]) * time + BetaLreal[3] * X1 + BetaLreal[4] * X2
    )



  temp_v <- subset(surv.data_v, select = c(Id, deathTimes, Event))
  pbc21_v <- tmerge(temp_v, temp_v, id = Id, endpt = event(deathTimes, Event))
  long.data1_v <- tmerge(pbc21_v, long.data_v,
    id = Id, lY1 = tdc(time, lY1),
    lY2 = tdc(time, lY2), lY3 = tdc(time, lY3), lY4 = tdc(time, lY4)
  )



  CRI_gold <- matrix(0, length(S), 2)

  for (kks in 1:length(S)) {
    s <- S[kks]
    Dt <- S[2]

    # Extract the baseline hazard function
    base_hazard <- basehaz(step2_gold, centered = FALSE)

    # Estimate cumulative hazard at time s=s for each individual

    cum_hazard_s <- rep(0, nrow(surv.data_v))
    for (i in 1:nrow(surv.data_v)) {
      individual_data <- long.data1_v[long.data1_v$Id == surv.data_v$Id[i] & long.data1_v$tstop <= s, ]
      if (nrow(individual_data) > 0) {
        linear_predictor <- sum(coef(step2_gold) * tail(individual_data[, c("lY1", "lY2", "lY3", "lY4")], 1))
        baseline_hazard_s <- base_hazard$hazard[base_hazard$time <= s]
        cum_hazard_s[i] <- tail(baseline_hazard_s, 1) * exp(linear_predictor)
      }
    }

    # Estimate survival probability at time s=s
    surv_prob_s <- exp(-cum_hazard_s)




    #####

    # Estimate cumulative hazard at time s=0.3 for each individual
    cum_hazard_sDt <- rep(0, nrow(surv.data_v))
    for (i in 1:nrow(surv.data_v)) {
      individual_data <- long.data1_v[long.data1_v$Id == surv.data_v$Id[i] & long.data1_v$tstop <= (s + Dt), ]
      if (nrow(individual_data) > 0) {
        linear_predictor <- sum(coef(step2_gold) * tail(individual_data[, c("lY1", "lY2", "lY3", "lY4")], 1))
        baseline_hazard_sDt <- base_hazard$hazard[base_hazard$time <= (s + Dt)]
        cum_hazard_sDt[i] <- tail(baseline_hazard_sDt, 1) * exp(linear_predictor)
      }
    }

    # Estimate survival probability at time s=0.3
    surv_prob_sDt <- exp(-cum_hazard_sDt)

    DP <- 1 - surv_prob_sDt / surv_prob_s


    length(DP)

    CRI_gold[kks, ] <- Criteria(
      s = s, t = Dt, Survt = surv.data_v$deathTimes,
      CR = surv.data_v$Event, P = DP, cause = 1
    )$Cri[, 1]
  }


  Gold <- list(Est_gold = step2_gold$coefficients, CRgold = CRgold, CRI_gold = CRI_gold)
```
#### R code for implementing the multi-marker joint model (MMJM)

```
  start2 <- Sys.time()
  i.jags <- function() {
    list(
      gamma1 = rnorm(1), gamma2 = rnorm(1), gamma3 = rnorm(1), gamma4 = rnorm(1),
      betaL1 = rnorm(dim(X)[3]), betaL2 = rnorm(dim(X)[3]), betaL3 = rnorm(dim(X)[3]),
      tau1 = 1, tau2 = 1, tau3 = 1, Omega = diag(runif(8))
    )
  }
  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y1, Y2 = Y2, Y3 = Y3, Y4 = Y4,
    XL1 = X, ZL1 = Z, XL2 = X, ZL2 = Z, XL3 = X, ZL3 = Z, XL4 = X, ZL4 = Z,
    death = death, mub = rep(0, 8), V = diag(1, 8), Nb = 8, zeros = rep(0, n),
    NbetasL = dim(X)[3], x1 = x1, x2 = x2, delta = Delta, s = peice, J = length(peice) + 1
  )

  parameters <- c(
    "gamma1", "gamma2", "gamma3", "gamma4", "betaL1", "betaL2", "betaL3", "betaL4",
    "sigma1", "sigma2", "sigma3", "sigma4", "Sigma", "h"
  )

  model.file <- "jm4.R"

  # file.show(model.file)
  sim123 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 10000, model.file = model.file
  )
  sim123

  end2 <- Sys.time()
  TimeMulti <- difftime(end2, start2, units = "mins")

  Sigma_m <- sim123$BUGSoutput$mean$Sigma
  betaL1_m <- sim123$BUGSoutput$mean$betaL1
  betaL2_m <- sim123$BUGSoutput$mean$betaL2
  betaL3_m <- sim123$BUGSoutput$mean$betaL3
  betaL4_m <- sim123$BUGSoutput$mean$betaL4

  gamma_m <- c(
    sim123$BUGSoutput$mean$gamma1, sim123$BUGSoutput$mean$gamma2, sim123$BUGSoutput$mean$gamma3,
    sim123$BUGSoutput$mean$gamma4
  )
  sigma_m <- c(
    sim123$BUGSoutput$mean$sigma1, sim123$BUGSoutput$mean$sigma2, sim123$BUGSoutput$mean$sigma3,
    sim123$BUGSoutput$mean$sigma4
  )
  h_m <- sim123$BUGSoutput$mean$h

  Est_mul <- gamma_m

  Est_mul_l <- sim123$BUGSoutput$summary[c(82:85), 3]
  Est_mul_u <- sim123$BUGSoutput$summary[c(82:85), 7]

  CRM <- rep(0, length(gamma_m))
  for (kk in 1:length(gamma_m)) {
    if (gammareal[kk] > Est_mul_l[kk] & gammareal[kk] < Est_mul_u[kk]) (CRM[kk] <- 1)
  }
```
#### Computing the dynamic prediction based on the multi-marker joint model (MMJM)

```
  CRI_MMJM <- matrix(0, length(S), 2)
  for (kks in 1:length(S)) {
    s <- S[kks]
    Dt <- S[2]
    Sigma_m <- sim123$BUGSoutput$mean$Sigma
    betaL1_m <- sim123$BUGSoutput$mean$betaL1
    betaL2_m <- sim123$BUGSoutput$mean$betaL2
    betaL3_m <- sim123$BUGSoutput$mean$betaL3
    betaL4_m <- sim123$BUGSoutput$mean$betaL4

    gamma_m <- c(sim123$BUGSoutput$mean$gamma1, sim123$BUGSoutput$mean$gamma2, sim123$BUGSoutput$mean$gamma3, sim123$BUGSoutput$mean$gamma4)
    sigma_m <- c(sim123$BUGSoutput$mean$sigma1, sim123$BUGSoutput$mean$sigma2, sim123$BUGSoutput$mean$sigma3, sim123$BUGSoutput$mean$sigma4)
    h_m <- sim123$BUGSoutput$mean$h

    ###########################   ###########################    ###########################   ###########################

    i.jags <- function() {
      list(b = matrix(0, n / 2, 2 * nmark))
    }

    M1 <- M
    for (ii in 1:n) {
      M1[ii] <- length(Y1[ii, ][t <= s][!is.na(Y1[ii, ][t <= s])])
    }
    M1[M1 == 0] <- 1

    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Time = rep(s, nsujet / 2),
      Y1 = Y1[(nsujet / 2 + 1):nsujet, ], Y2 = Y2[(nsujet / 2 + 1):nsujet, ],
      Y3 = Y3[(nsujet / 2 + 1):nsujet, ], Y4 = Y4[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      XL2 = X[(nsujet / 2 + 1):nsujet, , ], ZL2 = Z[(nsujet / 2 + 1):nsujet, , ], XL3 = X[(nsujet / 2 + 1):nsujet, , ],
      ZL3 = Z[(nsujet / 2 + 1):nsujet, , ],
      XL4 = X[(nsujet / 2 + 1):nsujet, , ],
      ZL4 = Z[(nsujet / 2 + 1):nsujet, , ],
      death = death[(nsujet / 2 + 1):nsujet], mub = rep(0, 2 * nmark), Nb = 2 * nmark, zeros = rep(0, n / 2),
      x1 = x1[(nsujet / 2 + 1):nsujet], x2 = x2[(nsujet / 2 + 1):nsujet], delta = Delta[(nsujet / 2 + 1):nsujet, ], s = peice,
      betaL1 = betaL1_m,
      betaL2 = betaL2_m,
      betaL3 = betaL3_m,
      betaL4 = betaL4_m,
      gamma1 = gamma_m[1], gamma2 = gamma_m[2], gamma3 = gamma_m[3], gamma4 = gamma_m[4],
      sigma1 = sigma_m[1], sigma2 = sigma_m[2], sigma3 = sigma_m[3], sigma4 = sigma_m[4],
      Sigma = Sigma_m, h = h_m
    )
    parameters <- c("b")
    model.file <- "jm4b.R"

    # file.show(model.file)
    sim123b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhatm <- sim123b$BUGSoutput$mean$b

    DPtotal <- c()
    for (i in inde1) {
      DPtotal[i] <- DP4b(
        s = s, Dt = Dt, betaL1 = betaL1_m, betaL2 = betaL2_m, betaL3 = betaL3_m, betaL4 = betaL4_m,
        gamma = gamma_m, sigma = sigma_m, h = h_m, Sigma = Sigma_m, peice = peice,
        y1 = Y1[i, ], y2 = Y2[i, ], y3 = Y3[i, ], y4 = Y4[i, ], x1 = x1[i], x2 = x2[i], b = bhatm[i - (nsujet / 2), ]
      )
    }

    DPtotal <- DPtotal[inde1]
    #########################

    CRI_MMJM[kks, ] <- Criteria(
      s = s, t = Dt, Survt = surv.data_v$deathTimes,
      CR = surv.data_v$Event, P = 1 - DPtotal, cause = 1
    )$Cri[, 1]
  }

  MMJM <- list(Est_mul = Est_mul, CRM = CRM, CRI_MMJM = CRI_MMJM)
```

#### R code for implementing the multiple two-stage (MTS)

```
  ###### MTS ######
  ### @@@@@@@@@ usual two-step
  start2 <- Sys.time()

  i.jags <- function() {
    list(
      betaL1 = rnorm(dim(X)[3]), tau1 = 1, Omega = diag(runif(2))
    )
  }

  parameters <- c("betaL1", "sigma1", "Sigma", "linearpred1", "b")

  model.file <- "Mix.txt"
  # file.show(model.file)
  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y1, XL1 = X, ZL1 = Z,
    mub = rep(0, 2), V = diag(1, 2), Nb = 2,
    NbetasL = dim(X)[3], x1 = x1, x2 = x2
  )
  mix1 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 5000, model.file = model.file
  )

  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y2, XL1 = X, ZL1 = Z,
    mub = rep(0, 2), V = diag(1, 2), Nb = 2,
    NbetasL = dim(X)[3], x1 = x1, x2 = x2
  )

  mix2 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 5000, model.file = model.file
  )

  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y3, XL1 = X, ZL1 = Z,
    mub = rep(0, 2), V = diag(1, 2), Nb = 2,
    NbetasL = dim(X)[3], x1 = x1, x2 = x2
  )
  mix3 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 5000, model.file = model.file
  )

  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y4, XL1 = X, ZL1 = Z,
    mub = rep(0, 2), V = diag(1, 2), Nb = 2,
    NbetasL = dim(X)[3], x1 = x1, x2 = x2
  )
  mix4 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 5000, model.file = model.file
  )

  ########### second stage
  Beta1 <- mix1$BUGSoutput$mean$betaL1
  Beta2 <- mix2$BUGSoutput$mean$betaL1
  Beta3 <- mix3$BUGSoutput$mean$betaL1
  Beta4 <- mix4$BUGSoutput$mean$betaL1

  b_1 <- mix1$BUGSoutput$sims.list$b
  b_2 <- mix2$BUGSoutput$sims.list$b
  b_3 <- mix3$BUGSoutput$sims.list$b
  b_4 <- mix4$BUGSoutput$sims.list$b

  u_mix <- apply(abind::abind(b_1, b_2, b_3, b_4), c(2, 3), mean)


  lPredY1 <- lPredY2 <- lPredY3 <- lPredY4 <- rep(0, length(long.data_t$Id))

  for (i in 1:length(long.data_t$Id)) {
    lPredY1[i] <- (Beta1[1] + u_mix[long.data_t$Id[i], 1]) + (Beta1[2] + u_mix[long.data_t$Id[i], 2]) * long.data_t$time[i] +
      Beta1[3] * long.data_t$X1[i] + Beta1[4] * long.data_t$X2[i]

    lPredY2[i] <- (Beta2[1] + u_mix[long.data_t$Id[i], 3]) + (Beta2[2] + u_mix[long.data_t$Id[i], 4]) * long.data_t$time[i] +
      Beta2[3] * long.data_t$X1[i] + Beta2[4] * long.data_t$X2[i]

    lPredY3[i] <- (Beta3[1] + u_mix[long.data_t$Id[i], 5]) + (Beta3[2] + u_mix[long.data_t$Id[i], 6]) * long.data_t$time[i] +
      Beta3[3] * long.data_t$X1[i] + Beta3[4] * long.data_t$X2[i]

    lPredY4[i] <- (Beta4[1] + u_mix[long.data_t$Id[i], 7]) + (Beta4[2] + u_mix[long.data_t$Id[i], 8]) * long.data_t$time[i] +
      Beta4[3] * long.data_t$X1[i] + Beta4[4] * long.data_t$X2[i]
  }


  temp <- subset(surv.data_t, select = c(Id, deathTimes, Event)) # baseline data
  pbc21 <- tmerge(temp, temp, id = Id, endpt = event(deathTimes, Event))
  long.data1 <- tmerge(pbc21, long.data_t,
    id = Id, lY1 = tdc(time, lPredY1),
    lY2 = tdc(time, lPredY2), lY3 = tdc(time, lPredY3),
    lY4 = tdc(time, lPredY4)
  )
  head(long.data)
  head(long.data1)

  cox_mts <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
    data = long.data1, id = Id, cluster = Id
  )
  #####################


  start3 <- Sys.time()
  indexchain <- 1:length(mix1$BUGSoutput$sims.list$betaL1[, 1])
  samplec <- sample(indexchain, Limp)
  betaL1 <- mix1$BUGSoutput$sims.list$betaL1[samplec, ]
  betaL2 <- mix2$BUGSoutput$sims.list$betaL1[samplec, ]
  betaL3 <- mix3$BUGSoutput$sims.list$betaL1[samplec, ]
  betaL4 <- mix4$BUGSoutput$sims.list$betaL1[samplec, ]

  b_1 <- mix1$BUGSoutput$sims.list$b
  b_2 <- mix2$BUGSoutput$sims.list$b
  b_3 <- mix3$BUGSoutput$sims.list$b
  b_4 <- mix4$BUGSoutput$sims.list$b

  u_mix <- abind::abind(b_1, b_2, b_3, b_4)[samplec, , ]

  SD2 <- matrix(0, length(samplec), nmark)
  TEM <- matrix(0, length(samplec), nmark)



  for (kkk in 1:length(samplec)) {
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u <- u_mix[kkk, , ]

    Beta1 <- betaL1[kkk, ]
    Beta2 <- betaL2[kkk, ]
    Beta3 <- betaL3[kkk, ]
    Beta4 <- betaL4[kkk, ]



    lPredY1 <- lPredY2 <- lPredY3 <- lPredY4 <- rep(0, length(long.data_t$Id))
    for (i in 1:length(long.data_t$Id)) {
      lPredY1[i] <- (Beta1[1] + u[long.data$Id[i], 1]) + (Beta1[2] + u[long.data$Id[i], 2]) * long.data$time[i] +
        Beta1[3] * long.data$X1[i] + Beta1[4] * long.data$X2[i]

      lPredY2[i] <- (Beta2[1] + u[long.data$Id[i], 3]) + (Beta2[2] + u[long.data$Id[i], 4]) * long.data$time[i] +
        Beta2[3] * long.data$X1[i] + Beta2[4] * long.data$X2[i]

      lPredY3[i] <- (Beta3[1] + u[long.data$Id[i], 5]) + (Beta3[2] + u[long.data$Id[i], 6]) * long.data$time[i] +
        Beta3[3] * long.data$X1[i] + Beta3[4] * long.data$X2[i]

      lPredY4[i] <- (Beta4[1] + u[long.data$Id[i], 7]) + (Beta4[2] + u[long.data$Id[i], 8]) * long.data$time[i] +
        Beta4[3] * long.data$X1[i] + Beta4[4] * long.data$X2[i]
    }






    temp <- subset(surv.data_t, select = c(Id, deathTimes, Event)) # baseline data
    pbc21 <- tmerge(temp, temp, id = Id, endpt = event(deathTimes, Event))
    long.data1 <- tmerge(pbc21, long.data_t,
      id = Id, lY1 = tdc(time, lPredY1),
      lY2 = tdc(time, lPredY2), lY3 = tdc(time, lPredY3), lY4 = tdc(time, lPredY4)
    )

    RRR <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
      data = long.data1, cluster = Id
    )




    TEM[kkk, ] <- summary(RRR)$coefficients[, 1]
    SD2[kkk, ] <- summary(RRR)$coefficients[, 3]
  }
  #########
  Wv <- apply(SD2^2, 2, mean)
  Bv <- apply(TEM, 2, var)
  sdnew <- sqrt(Wv + (Limp + 1) / Limp * Bv)

  ############# Results2 ########
  Est_mts_l <- cox_mts$coefficients - qnorm(.975) * sdnew
  Est_mts_u <- cox_mts$coefficients + qnorm(.975) * sdnew
  CRMTS <- rep(0, length(gammareal))

  for (kk in 1:length(gammareal)) {
    if (gammareal[kk] > Est_mts_l[kk] & gammareal[kk] < Est_mts_u[kk]) (CRMTS[kk] <- 1)
  }
```
#### Computing the dynamic prediction based on  the multiple two-stage (MTS)

```
  ### generation of random effects
  CRI_MTS <- matrix(0, length(S), 2)
  for (kks in 1:length(S)) {
    s <- S[kks]
    Dt <- S[2]
    i.jags <- function() {
      list(b = matrix(0, n / 2, 2))
    }

    parameters <- c("b")

    model.file <- "b_sts.R"
    # file.show(model.file)
    betaL1_mix <- mix1$BUGSoutput$mean$betaL1
    betaL2_mix <- mix2$BUGSoutput$mean$betaL1
    betaL3_mix <- mix3$BUGSoutput$mean$betaL1
    betaL4_mix <- mix4$BUGSoutput$mean$betaL1

    Sigma_mix <- as.matrix(Matrix::bdiag(
      mix1$BUGSoutput$mean$Sigma,
      mix2$BUGSoutput$mean$Sigma,
      mix3$BUGSoutput$mean$Sigma,
      mix4$BUGSoutput$mean$Sigma
    ))
    sigma_mix <- c(mix1$BUGSoutput$mean$sigma1, mix2$BUGSoutput$mean$sigma1, mix3$BUGSoutput$mean$sigma1, mix4$BUGSoutput$mean$sigma1)

    ########### Y1
    bhat_sts <- matrix(0, n / 2, 2 * nmark)
    M1 <- M
    for (ii in 1:n) {
      M1[ii] <- length(Y1[ii, ][t <= s][!is.na(Y1[ii, ][t <= s])])
    }
    M1[M1 == 0] <- 1
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Y1 = Y1[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      mub = rep(0, 2), Nb = 2,
      betaL1 = betaL1_mix,
      sigma1 = sigma_mix[1],
      Sigma = Sigma_mix[1:2, 1:2]
    )

    mix1b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_sts[, 1:2] <- mix1b$BUGSoutput$mean$b
    ########### Y2
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Y1 = Y2[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      mub = rep(0, 2), Nb = 2,
      betaL1 = betaL2_mix,
      sigma1 = sigma_mix[2],
      Sigma = Sigma_mix[3:4, 3:4]
    )

    mix2b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_sts[, 3:4] <- mix2b$BUGSoutput$mean$b
    ########### Y3
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Y1 = Y3[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      mub = rep(0, 2), Nb = 2,
      betaL1 = betaL3_mix,
      sigma1 = sigma_mix[3],
      Sigma = Sigma_mix[5:6, 5:6]
    )

    mix3b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_sts[, 5:6] <- mix3b$BUGSoutput$mean$b
    ########### Y4
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Y1 = Y4[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      mub = rep(0, 2), Nb = 2,
      betaL1 = betaL3_mix,
      sigma1 = sigma_mix[3],
      Sigma = Sigma_mix[5:6, 5:6]
    )

    mix4b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_sts[, 7:8] <- mix4b$BUGSoutput$mean$b

    ###########
    # Calculate linear predictors for the validation set



    bhat_sts <- rbind(bhat_sts, bhat_sts)

    long.data_v <- long.data_v %>%
      mutate(
        lY1 = (betaL1[1] + bhat_sts[Id, 1]) + (betaL1[2] + bhat_sts[Id, 2]) * time + betaL1[3] * X1 + betaL1[4] * X2,
        lY2 = (betaL2[1] + bhat_sts[Id, 3]) + (betaL2[2] + bhat_sts[Id, 4]) * time + betaL2[3] * X1 + betaL2[4] * X2,
        lY3 = (betaL3[1] + bhat_sts[Id, 5]) + (betaL3[2] + bhat_sts[Id, 6]) * time + betaL3[3] * X1 + betaL3[4] * X2,
        lY4 = (betaL4[1] + bhat_sts[Id, 7]) + (betaL4[2] + bhat_sts[Id, 8]) * time + betaL4[3] * X1 + betaL4[4] * X2
      )



    temp_v <- subset(surv.data_v, select = c(Id, deathTimes, Event))
    pbc21_v <- tmerge(temp_v, temp_v, id = Id, endpt = event(deathTimes, Event))
    long.data1_v <- tmerge(pbc21_v, long.data_v,
      id = Id, lY1 = tdc(time, lY1),
      lY2 = tdc(time, lY2), lY3 = tdc(time, lY3), lY4 = tdc(time, lY4)
    )



    # Extract the baseline hazard function
    base_hazard <- basehaz(cox_mts, centered = FALSE)

    # Estimate cumulative hazard at time s=s for each individual

    cum_hazard_s <- rep(0, nrow(surv.data_v))
    for (i in 1:nrow(surv.data_v)) {
      individual_data <- long.data1_v[long.data1_v$Id == surv.data_v$Id[i] & long.data1_v$tstop <= s, ]
      if (nrow(individual_data) > 0) {
        linear_predictor <- sum(coef(cox_mts) * tail(individual_data[, c("lY1", "lY2", "lY3", "lY4")], 1))
        baseline_hazard_s <- base_hazard$hazard[base_hazard$time <= s]
        cum_hazard_s[i] <- tail(baseline_hazard_s, 1) * exp(linear_predictor)
      }
    }

    # Estimate survival probability at time s=s
    surv_prob_s <- exp(-cum_hazard_s)




    #####

    # Estimate cumulative hazard at time s=0.3 for each individual
    cum_hazard_sDt <- rep(0, nrow(surv.data_v))
    for (i in 1:nrow(surv.data_v)) {
      individual_data <- long.data1_v[long.data1_v$Id == surv.data_v$Id[i] & long.data1_v$tstop <= (s + Dt), ]
      if (nrow(individual_data) > 0) {
        linear_predictor <- sum(coef(cox_mts) * tail(individual_data[, c("lY1", "lY2", "lY3", "lY4")], 1))
        baseline_hazard_sDt <- base_hazard$hazard[base_hazard$time <= (s + Dt)]
        cum_hazard_sDt[i] <- tail(baseline_hazard_sDt, 1) * exp(linear_predictor)
      }
    }

    # Estimate survival probability at time s=0.3
    surv_prob_sDt <- exp(-cum_hazard_sDt)

    DP <- 1 - surv_prob_sDt / surv_prob_s


    length(DP)

    CRI_MTS[kks, ] <- Criteria(
      s = s, t = Dt, Survt = surv.data_v$deathTimes,
      CR = surv.data_v$Event, P = DP, cause = 1
    )$Cri[, 1]
  }

  MTS <- list(Est = cox_mts$coefficients, CRMTS = CRMTS, CRI_MTS = CRI_MTS)
```
#### R code for implementing TSJM (without TSJM R package)

```
  start1 <- Sys.time()

  i.jags <- function() {
    list(
      gamma1 = rnorm(1),
      betaL1 = rnorm(4), tau1 = 1, Omega = diag(runif(2))
    )
  }

  parameters <- c(
    "linearpred1", "betaL1", "b",
    "gamma1", "sigma1", "Sigma", "h"
  )

  model.file <- "JM1p.txt"
  # file.show(model.file)
  nsujet <- n
  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y1,
    XL1 = X, ZL1 = Z,
    death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
    NbetasL = 4, x1 = x1, x2 = x2, delta = Delta, s = peice, J = length(peice) + 1
  )


  sim1 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 10000, model.file = model.file
  )
  ###########
  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y2,
    XL1 = X, ZL1 = Z,
    death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
    NbetasL = 4, x1 = x1, x2 = x2, delta = Delta, s = peice, J = length(peice) + 1
  )


  sim2 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 10000, model.file = model.file
  )
  ###########
  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y3,
    XL1 = X, ZL1 = Z,
    death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
    NbetasL = 4, x1 = x1, x2 = x2, delta = Delta, s = peice, J = length(peice) + 1
  )

  sim3 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 10000, model.file = model.file
  )
  ###########
  d.jags <- list(
    n = nsujet / 2, M = M, Time = st, Y1 = Y4,
    XL1 = X, ZL1 = Z,
    death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
    NbetasL = 4, x1 = x1, x2 = x2, delta = Delta, s = peice, J = length(peice) + 1
  )

  sim4 <- R2jags::jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 10000, model.file = model.file
  )
  sim1$BUGSoutput$mean$gamma1
  sim2$BUGSoutput$mean$gamma1
  sim3$BUGSoutput$mean$gamma1
  sim4$BUGSoutput$mean$gamma1
  ##########################################################
  Beta1 <- sim1$BUGSoutput$mean$betaL1
  Beta2 <- sim2$BUGSoutput$mean$betaL1
  Beta3 <- sim3$BUGSoutput$mean$betaL1
  Beta4 <- sim4$BUGSoutput$mean$betaL1

  b_1 <- sim1$BUGSoutput$sims.list$b
  b_2 <- sim2$BUGSoutput$sims.list$b
  b_3 <- sim3$BUGSoutput$sims.list$b
  b_4 <- sim4$BUGSoutput$sims.list$b

  u_sim <- apply(abind::abind(b_1, b_2, b_3, b_4), c(2, 3), mean)


  lPredY1 <- lPredY2 <- lPredY3 <- lPredY4 <- rep(0, length(long.data_t$Id))

  for (i in 1:length(long.data_t$Id)) {
    lPredY1[i] <- (Beta1[1] + u_sim[long.data_t$Id[i], 1]) + (Beta1[2] + u_sim[long.data_t$Id[i], 2]) * long.data_t$time[i] +
      Beta1[3] * long.data_t$X1[i] + Beta1[4] * long.data_t$X2[i]

    lPredY2[i] <- (Beta2[1] + u_sim[long.data_t$Id[i], 3]) + (Beta2[2] + u_sim[long.data_t$Id[i], 4]) * long.data_t$time[i] +
      Beta2[3] * long.data_t$X1[i] + Beta2[4] * long.data_t$X2[i]

    lPredY3[i] <- (Beta3[1] + u_sim[long.data_t$Id[i], 5]) + (Beta3[2] + u_sim[long.data_t$Id[i], 6]) * long.data_t$time[i] +
      Beta3[3] * long.data_t$X1[i] + Beta3[4] * long.data_t$X2[i]

    lPredY4[i] <- (Beta4[1] + u_sim[long.data_t$Id[i], 7]) + (Beta4[2] + u_sim[long.data_t$Id[i], 8]) * long.data_t$time[i] +
      Beta4[3] * long.data_t$X1[i] + Beta4[4] * long.data_t$X2[i]
  }


  temp <- subset(surv.data_t, select = c(Id, deathTimes, Event)) # baseline data
  pbc21 <- tmerge(temp, temp, id = Id, endpt = event(deathTimes, Event))
  long.data1 <- tmerge(pbc21, long.data_t,
    id = Id, lY1 = tdc(time, lPredY1),
    lY2 = tdc(time, lPredY2), lY3 = tdc(time, lPredY3),
    lY4 = tdc(time, lPredY4)
  )
  head(long.data)
  head(long.data1)

  cox_model_tsjm <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
    data = long.data1, id = Id, cluster = Id
  )
  #####################
  ### ###@@@@@@@@@ Step 2-Rubin formula ####
  start3 <- Sys.time()
  indexchain <- 1:length(sim1$BUGSoutput$sims.list$betaL1[, 1])
  samplec <- sample(indexchain, Limp)
  betaL1 <- sim1$BUGSoutput$sims.list$betaL1[samplec, ]
  betaL2 <- sim2$BUGSoutput$sims.list$betaL1[samplec, ]
  betaL3 <- sim3$BUGSoutput$sims.list$betaL1[samplec, ]
  betaL4 <- sim4$BUGSoutput$sims.list$betaL1[samplec, ]

  b_1 <- sim1$BUGSoutput$sims.list$b
  b_2 <- sim2$BUGSoutput$sims.list$b
  b_3 <- sim3$BUGSoutput$sims.list$b
  b_4 <- sim4$BUGSoutput$sims.list$b

  u_ujm <- abind::abind(b_1, b_2, b_3, b_4)[samplec, , ]

  SD2 <- matrix(0, length(samplec), nmark)
  TEM <- matrix(0, length(samplec), nmark)



  for (kkk in 1:length(samplec)) {
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u <- u_ujm[kkk, , ]

    Beta1 <- betaL1[kkk, ]
    Beta2 <- betaL2[kkk, ]
    Beta3 <- betaL3[kkk, ]
    Beta4 <- betaL4[kkk, ]



    lPredY1 <- lPredY2 <- lPredY3 <- lPredY4 <- rep(0, length(long.data_t$Id))
    for (i in 1:length(long.data_t$Id)) {
      lPredY1[i] <- (Beta1[1] + u[long.data$Id[i], 1]) + (Beta1[2] + u[long.data$Id[i], 2]) * long.data$time[i] +
        Beta1[3] * long.data$X1[i] + Beta1[4] * long.data$X2[i]

      lPredY2[i] <- (Beta2[1] + u[long.data$Id[i], 3]) + (Beta2[2] + u[long.data$Id[i], 4]) * long.data$time[i] +
        Beta2[3] * long.data$X1[i] + Beta2[4] * long.data$X2[i]

      lPredY3[i] <- (Beta3[1] + u[long.data$Id[i], 5]) + (Beta3[2] + u[long.data$Id[i], 6]) * long.data$time[i] +
        Beta3[3] * long.data$X1[i] + Beta3[4] * long.data$X2[i]

      lPredY4[i] <- (Beta4[1] + u[long.data$Id[i], 7]) + (Beta4[2] + u[long.data$Id[i], 8]) * long.data$time[i] +
        Beta4[3] * long.data$X1[i] + Beta4[4] * long.data$X2[i]
    }






    temp <- subset(surv.data_t, select = c(Id, deathTimes, Event)) # baseline data
    pbc21 <- tmerge(temp, temp, id = Id, endpt = event(deathTimes, Event))
    long.data1 <- tmerge(pbc21, long.data_t,
      id = Id, lY1 = tdc(time, lPredY1),
      lY2 = tdc(time, lPredY2), lY3 = tdc(time, lPredY3), lY4 = tdc(time, lPredY4)
    )
    head(long.data)
    head(long.data1)

    RRR <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
      data = long.data1, cluster = Id
    )




    TEM[kkk, ] <- summary(RRR)$coefficients[, 1]
    SD2[kkk, ] <- summary(RRR)$coefficients[, 3]
  }

  Wv <- apply(SD2^2, 2, mean)
  Bv <- apply(TEM, 2, var)
  sdnew <- sqrt(Wv + (Limp + 1) / Limp * Bv)

  ############# Results2 ########
  Est_unijm_l <- cox_model_tsjm$coefficients - qnorm(.975) * sdnew
  Est_unijm_u <- cox_model_tsjm$coefficients + qnorm(.975) * sdnew
  CRUJM <- rep(0, length(gammareal))

  for (kk in 1:length(gammareal)) {
    if (gammareal[kk] > Est_unijm_l[kk] & gammareal[kk] < Est_unijm_u[kk]) (CRUJM[kk] <- 1)
  }


  Est_unijm_l <- apply(TEM, 2, mean) - qnorm(.975) * sdnew
  Est_unijm_u <- apply(TEM, 2, mean) + qnorm(.975) * sdnew

  CRUJM1 <- rep(0, length(gammareal))

  for (kk in 1:length(gammareal)) {
    if (gammareal[kk] > Est_unijm_l[kk] & gammareal[kk] < Est_unijm_u[kk]) (CRUJM1[kk] <- 1)
  }
```
#### Computing the dynamic prediction based on the TSJM

```
  ### generation of random effects
  CRI_TSJM <- matrix(0, length(S), 2)
  for (kks in 1:length(S)) {
    s <- S[kks]
    Dt <- S[2]
    i.jags <- function() {
      list(b = matrix(0, n / 2, 2))
    }

    betaL1 <- sim1$BUGSoutput$sims.list$betaL1[samplec, ]
    betaL2 <- sim2$BUGSoutput$sims.list$betaL1[samplec, ]
    betaL3 <- sim3$BUGSoutput$sims.list$betaL1[samplec, ]
    betaL4 <- sim4$BUGSoutput$sims.list$betaL1[samplec, ]


    parameters <- c("b")

    ########### Y1
    i.jags <- function() {
      list(b = matrix(0, n / 2, 2))
    }

    parameters <- c("b")

    model.file <- "b_unijm.R"
    # file.show(model.file)
    ########### Y1
    ###########################
    #################################
    ############   ############   ############  tsjm_censored  ############   ############   ############
    M1 <- M
    for (ii in 1:n) {
      M1[ii] <- length(Y1[ii, ][t <= s][!is.na(Y1[ii, ][t <= s])])
    }
    M1[M1 == 0] <- 1
    ########### Y1
    bhat_tsjm <- matrix(0, n / 2, 2 * nmark)
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Time = rep(s, (nsujet / 2)), Y1 = Y1[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      death = death[(nsujet / 2 + 1):nsujet], mub = rep(0, 2), zeros = rep(0, n / 2),
      x1 = x1[(nsujet / 2 + 1):nsujet], x2 = x2[(nsujet / 2 + 1):nsujet], delta = Delta[(nsujet / 2 + 1):nsujet, ], s = peice,
      betaL1 = sim1$BUGSoutput$mean$betaL1,
      gamma1 = sim1$BUGSoutput$mean$gamma1,
      sigma1 = sim1$BUGSoutput$mean$sigma1,
      Sigma = sim1$BUGSoutput$mean$Sigma, h = sim1$BUGSoutput$mean$h
    )

    sim1b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_tsjm[, 1:2] <- sim1b$BUGSoutput$mean$b
    ########### Y2
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Time = rep(s, (nsujet / 2)), Y1 = Y2[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      death = death[(nsujet / 2 + 1):nsujet], mub = rep(0, 2), zeros = rep(0, n / 2),
      x1 = x1[(nsujet / 2 + 1):nsujet], x2 = x2[(nsujet / 2 + 1):nsujet], delta = Delta[(nsujet / 2 + 1):nsujet, ], s = peice,
      betaL1 = sim2$BUGSoutput$mean$betaL1,
      gamma1 = sim2$BUGSoutput$mean$gamma1,
      sigma1 = sim2$BUGSoutput$mean$sigma1,
      Sigma = sim2$BUGSoutput$mean$Sigma, h = sim2$BUGSoutput$mean$h
    )

    sim2b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_tsjm[, 3:4] <- sim2b$BUGSoutput$mean$b
    ########### Y3
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Time = rep(s, (nsujet / 2)), Y1 = Y3[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      death = death[(nsujet / 2 + 1):nsujet], mub = rep(0, 2), zeros = rep(0, n / 2),
      x1 = x1[(nsujet / 2 + 1):nsujet], x2 = x2[(nsujet / 2 + 1):nsujet], delta = Delta[(nsujet / 2 + 1):nsujet, ], s = peice,
      betaL1 = sim3$BUGSoutput$mean$betaL1,
      gamma1 = sim3$BUGSoutput$mean$gamma1,
      sigma1 = sim3$BUGSoutput$mean$sigma1,
      Sigma = sim3$BUGSoutput$mean$Sigma, h = sim3$BUGSoutput$mean$h
    )

    sim3b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_tsjm[, 5:6] <- sim3b$BUGSoutput$mean$b
    ########### Y4
    d.jags <- list(
      n = nsujet / 2, M = M1[(nsujet / 2 + 1):nsujet], Time = rep(s, (nsujet / 2)), Y1 = Y4[(nsujet / 2 + 1):nsujet, ],
      XL1 = X[(nsujet / 2 + 1):nsujet, , ], ZL1 = Z[(nsujet / 2 + 1):nsujet, , ],
      death = death[(nsujet / 2 + 1):nsujet], mub = rep(0, 2), zeros = rep(0, n / 2),
      x1 = x1[(nsujet / 2 + 1):nsujet], x2 = x2[(nsujet / 2 + 1):nsujet], delta = Delta[(nsujet / 2 + 1):nsujet, ], s = peice,
      betaL1 = sim4$BUGSoutput$mean$betaL1,
      gamma1 = sim4$BUGSoutput$mean$gamma1,
      sigma1 = sim4$BUGSoutput$mean$sigma1,
      Sigma = sim4$BUGSoutput$mean$Sigma, h = sim4$BUGSoutput$mean$h
    )

    sim4b <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = 1,
      n.iter = 5000, model.file = model.file
    )
    bhat_tsjm[, 7:8] <- sim4b$BUGSoutput$mean$b

    ###########
    # Calculate linear predictors for the validation set



    bhat_tsjm <- rbind(bhat_tsjm, bhat_tsjm)

    long.data_v <- long.data_v %>%
      mutate(
        lY1 = (betaL1[1] + bhat_tsjm[Id, 1]) + (betaL1[2] + bhat_tsjm[Id, 2]) * time + betaL1[3] * X1 + betaL1[4] * X2,
        lY2 = (betaL2[1] + bhat_tsjm[Id, 3]) + (betaL2[2] + bhat_tsjm[Id, 4]) * time + betaL2[3] * X1 + betaL2[4] * X2,
        lY3 = (betaL3[1] + bhat_tsjm[Id, 5]) + (betaL3[2] + bhat_tsjm[Id, 6]) * time + betaL3[3] * X1 + betaL3[4] * X2,
        lY4 = (betaL4[1] + bhat_tsjm[Id, 7]) + (betaL4[2] + bhat_tsjm[Id, 8]) * time + betaL4[3] * X1 + betaL4[4] * X2
      )



    temp_v <- subset(surv.data_v, select = c(Id, deathTimes, Event))
    pbc21_v <- tmerge(temp_v, temp_v, id = Id, endpt = event(deathTimes, Event))
    long.data1_v <- tmerge(pbc21_v, long.data_v,
      id = Id, lY1 = tdc(time, lY1),
      lY2 = tdc(time, lY2), lY3 = tdc(time, lY3), lY4 = tdc(time, lY4)
    )



    # Extract the baseline hazard function
    base_hazard <- basehaz(cox_model_tsjm, centered = FALSE)

    # Estimate cumulative hazard at time s=s for each individual

    cum_hazard_s <- rep(0, nrow(surv.data_v))
    for (i in 1:nrow(surv.data_v)) {
      individual_data <- long.data1_v[long.data1_v$Id == surv.data_v$Id[i] & long.data1_v$tstop <= s, ]
      if (nrow(individual_data) > 0) {
        linear_predictor <- sum(coef(cox_model_tsjm) * tail(individual_data[, c("lY1", "lY2", "lY3", "lY4")], 1))
        baseline_hazard_s <- base_hazard$hazard[base_hazard$time <= s]
        cum_hazard_s[i] <- tail(baseline_hazard_s, 1) * exp(linear_predictor)
      }
    }

    # Estimate survival probability at time s=s
    surv_prob_s <- exp(-cum_hazard_s)

    #####

    # Estimate cumulative hazard at time s=0.3 for each individual
    cum_hazard_sDt <- rep(0, nrow(surv.data_v))
    for (i in 1:nrow(surv.data_v)) {
      individual_data <- long.data1_v[long.data1_v$Id == surv.data_v$Id[i] & long.data1_v$tstop <= (s + Dt), ]
      if (nrow(individual_data) > 0) {
        linear_predictor <- sum(coef(cox_model_tsjm) * tail(individual_data[, c("lY1", "lY2", "lY3", "lY4")], 1))
        baseline_hazard_sDt <- base_hazard$hazard[base_hazard$time <= (s + Dt)]
        cum_hazard_sDt[i] <- tail(baseline_hazard_sDt, 1) * exp(linear_predictor)
      }
    }

    # Estimate survival probability at time s=0.3
    surv_prob_sDt <- exp(-cum_hazard_sDt)

    DP <- 1 - surv_prob_sDt / surv_prob_s


    CRI_TSJM[kks, ] <- Criteria(
      s = s, t = Dt, Survt = surv.data_v$deathTimes,
      CR = surv.data_v$Event, P = DP, cause = 1
    )$Cri[, 1]
  }


  TSJM <- list(Est = cox_model_tsjm$coefficients, CRUJM = CRUJM, CRI_TSJM = CRI_TSJM, Est1 = apply(TEM, 2, mean), CRUJM1 = CRUJM1)



  list(Gold = Gold, MMJM = MMJM, MTS = MTS, TSJM = TSJM)
}

stopCluster(cl)

save(resultsss, file = "result_d1hh2.RData")


```
