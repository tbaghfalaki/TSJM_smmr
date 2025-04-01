The R code for computing the computational time based on TSJM package. For this purpose we consider TSC function of TSJM package. Also, we consider 100 imputations for computationg statnaded deviation of the assoctian parameters. For more detail of this function refer to [here](https://github.com/tbaghfalaki/TSJM/blob/main/Exam3.md) and [here](https://github.com/tbaghfalaki/TSJM/blob/main/Exam4.md).


```
rm(list = ls())
# Load necessary libraries
library(survival) # Survival analysis
library(PermAlgo) # Permutation algorithm for survival times with time-varying covariates
library(mvtnorm) # Multivariate normal distribution
library(TSJM) # Joint modeling for longitudinal and survival data


TS0 <- TS1 <- TS2 <- TS3 <- list()
Timetsjm0 <- Timetsjm1 <- Timetsjm2 <- Timetsjm3 <- c()
cox_model_tsjm0 <- cox_model_tsjm1 <- cox_model_tsjm2 <- cox_model_tsjm3 <- list()

nsujet <- 1000
NN <- 100
for (ij in 1:NN) {
  set.seed(ij)
  # Y1 (continuous)
  Beta1 <- c(-0.5, 0.5, 0.5, 0.5)
  Beta2 <- c(-0.5, 0.5, 0.5, 0.5)
  Beta3 <- c(-0.5, 0.5, 0.5, 0.5)
  Beta4 <- c(-0.5, 0.5, 0.5, 0.5)
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
    matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2),
    matrix(c(1, rho, rho, 1), 2, 2)
  )
  Sigma[Sigma == 0] <- rho1
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
  
  
  w1 <- rnorm(nsujet, 1, 0.5)
  W1 <- rep(x1, each = nmesindiv) # continuous covariate
  w2 <- rbinom(nsujet, 1, 0.5)
  W2 <- rep(x2, each = nmesindiv) # binary covariate
  
  
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
  survDat$w1 <- w1
  survDat$w2 <- w2
  
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
  MM2000 <- median(M[surv.data$Event == 0])
  MM1 <- c(mean(M), MM11, MM10)
  MM2 <- c(median(M), MM21, MM2000)
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
  
  INDTRAIN <- 1:(nsujet / 2)
  INDVALID <- surv.data$Id[-INDTRAIN]
  
  dataLong_t <- subset(
    long.data,
    long.data$Id %in% INDTRAIN
  )
  dataSurv_t <- subset(
    surv.data,
    surv.data$Id %in% INDTRAIN
  )
  
  dataLong_v <- subset(
    long.data,
    long.data$Id %in% INDVALID
  )
  dataSurv_v <- subset(
    surv.data,
    surv.data$Id %in% INDVALID
  )
  
  
  head(dataLong_t)
  head(dataSurv_t)
  
  
  formFixed <- list(Y1 ~ time, Y2 ~ time, Y3 ~ time, Y4 ~ time)
  formRandom <- list(~time, ~time, ~time, ~time)
  formGroup <- list(~Id, ~Id, ~Id, ~Id)
  formSurv <- survival::Surv(deathTimes, Event) ~ w1 + w2
  
  
  
  model <- list("linear", "linear", "linear", "linear")
  
  start0 <- Sys.time()
  # TSC function
  TS0[[ij]] <- TSC(formFixed, formRandom, formGroup, formSurv,
                   nmark = 4, K1 = 15,
                   model = model, n.chains1 = 1, n.iter1 = 5000, n.burnin1 = 2000, n.thin1 = 1,
                   Obstime = "time", ncl = 4, Limp = 100,
                   DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
  )
  
  # Merge survival data for modeling
  surv_data <- survival::tmerge(dataSurv_t, dataSurv_t, id = Id, endpt = event(deathTimes, Event))
  
  head(dataSurv_t)
  # Merge longitudinal data with survival data
  long.data1 <- survival::tmerge(surv_data, dataLong_t,
                                 id = Id, lY1 = tdc(time, TS0[[ij]]$lPredY[, 1]),
                                 lY2 = tdc(time, TS0[[ij]]$lPredY[, 2]), lY3 = tdc(time, TS0[[ij]]$lPredY[, 3]),
                                 lY4 = tdc(time, TS0[[ij]]$lPredY[, 4])
  )
  
  # Fit a Cox proportional hazards model using the joint model's longitudinal predictions
  cox_model_tsjm0[[ij]] <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
                                 data = long.data1, id = Id
  )
  
  
  # Prepare data for multiple imputation
  surv_data_new <- survival::tmerge(dataSurv_t, dataSurv_t, id = Id, endpt = event(deathTimes, Event))
  
  
  Limp <- 100
  SD2 <- matrix(0, Limp, length(cox_model_tsjm0[[ij]]$coefficients))
  TEM <- matrix(0, Limp, length(cox_model_tsjm0[[ij]]$coefficients))
  
  # Perform multiple imputation to estimate variability of coefficients
  for (k in 1:Limp) {
    lPredY <- TS0[[ij]]$LPredY[[k]]
    long.data1 <- survival::tmerge(surv_data_new, dataLong_t,
                                   id = Id,
                                   lY1 = tdc(time, TS0[[ij]]$lPredY[, 1]),
                                   lY2 = tdc(time, TS0[[ij]]$lPredY[, 2]),
                                   lY3 = tdc(time, TS0[[ij]]$lPredY[, 3]),
                                   lY4 = tdc(time, TS0[[ij]]$lPredY[, 4])
    )
    
    cox_model <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
                       data = long.data1, id = Id
    )
    
    TEM[k, ] <- summary(cox_model)$coefficients[, 1]
    SD2[k, ] <- summary(cox_model)$coefficients[, 3]
  }
  
  # Compute the within-imputation variance (Wv) and between-imputation variance (Bv)
  Wv <- apply(SD2^2, 2, mean)
  Bv <- apply(TEM, 2, var)
  sdnew <- sqrt(Wv + (Limp + 1) / Limp * Bv)
  
  # Combine results into a single summary table
  Res <- cbind(
    cox_model_tsjm0[[ij]]$coefficients, sdnew,
    cox_model_tsjm0[[ij]]$coefficients - qnorm(.95) * sdnew,
    cox_model_tsjm0[[ij]]$coefficients + qnorm(.95) * sdnew
  )
  
  colnames(Res) <- c("coefficients", "sd", "L_CI", "U_CI")
  print(Res)
  
  end0 <- Sys.time()
  Timetsjm0[ij] <- difftime(end0, start0, units = "mins")
  TS0
  
  
  ##################################################
  
  start1 <- Sys.time()
  
  # TSC function
  TS1[[ij]] <- TSC(formFixed, formRandom, formGroup, formSurv,
                   nmark = 4, K1 = 15,
                   model = model, n.chains1 = 1, n.iter1 = 5000, n.burnin1 = 2000, n.thin1 = 1,
                   Obstime = "time", ncl = 2, Limp = 100,
                   DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
  )
  
  # Merge survival data for modeling
  surv_data <- survival::tmerge(dataSurv_t, dataSurv_t, id = Id, endpt = event(deathTimes, Event))
  
  head(dataSurv_t)
  # Merge longitudinal data with survival data
  long.data1 <- survival::tmerge(surv_data, dataLong_t,
                                 id = Id, lY1 = tdc(time, TS1[[ij]]$lPredY[, 1]),
                                 lY2 = tdc(time, TS1[[ij]]$lPredY[, 2]), lY3 = tdc(time, TS1[[ij]]$lPredY[, 3]),
                                 lY4 = tdc(time, TS1[[ij]]$lPredY[, 4])
  )
  
  # Fit a Cox proportional hazards model using the joint model's longitudinal predictions
  cox_model_tsjm1[[ij]] <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
                                 data = long.data1, id = Id
  )
  
  
  
  
  
  
  # Prepare data for multiple imputation
  surv_data_new <- survival::tmerge(dataSurv_t, dataSurv_t, id = Id, endpt = event(deathTimes, Event))
  
  
  Limp <- 100
  SD2 <- matrix(0, Limp, length(cox_model_tsjm1[[ij]]$coefficients))
  TEM <- matrix(0, Limp, length(cox_model_tsjm1[[ij]]$coefficients))
  
  # Perform multiple imputation to estimate variability of coefficients
  for (k in 1:Limp) {
    lPredY <- TS1[[ij]]$LPredY[[k]]
    long.data1 <- survival::tmerge(surv_data_new, dataLong_t,
                                   id = Id,
                                   lY1 = tdc(time, TS1[[ij]]$lPredY[, 1]),
                                   lY2 = tdc(time, TS1[[ij]]$lPredY[, 2]),
                                   lY3 = tdc(time, TS1[[ij]]$lPredY[, 3]),
                                   lY4 = tdc(time, TS1[[ij]]$lPredY[, 4])
    )
    
    cox_model <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
                       data = long.data1, id = Id
    )
    
    TEM[k, ] <- summary(cox_model)$coefficients[, 1]
    SD2[k, ] <- summary(cox_model)$coefficients[, 3]
  }
  
  # Compute the within-imputation variance (Wv) and between-imputation variance (Bv)
  Wv <- apply(SD2^2, 2, mean)
  Bv <- apply(TEM, 2, var)
  sdnew <- sqrt(Wv + (Limp + 1) / Limp * Bv)
  
  # Combine results into a single summary table
  Res <- cbind(
    cox_model_tsjm1[[ij]]$coefficients, sdnew,
    cox_model_tsjm1[[ij]]$coefficients - qnorm(.95) * sdnew,
    cox_model_tsjm1[[ij]]$coefficients + qnorm(.95) * sdnew
  )
  
  colnames(Res) <- c("coefficients", "sd", "L_CI", "U_CI")
  print(Res)
  
  end1 <- Sys.time()
  Timetsjm1[ij] <- difftime(end1, start1, units = "mins")
  TS1
  
  ##################################################
  
  
  start2 <- Sys.time()
  # TSC function
  TS2[[ij]] <- TSC(formFixed, formRandom, formGroup, formSurv,
                   nmark = 4, K1 = 15,
                   model = model, n.chains1 = 1, n.iter1 = 5000, n.burnin1 = 2000, n.thin1 = 1,
                   Obstime = "time", ncl = 1, Limp = 100,
                   DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
  )
  
  
  # Merge survival data for modeling
  surv_data <- survival::tmerge(dataSurv_t, dataSurv_t, id = Id, endpt = event(deathTimes, Event))
  
  head(dataSurv_t)
  # Merge longitudinal data with survival data
  long.data1 <- survival::tmerge(surv_data, dataLong_t,
                                 id = Id, lY1 = tdc(time, TS2[[ij]]$lPredY[, 1]),
                                 lY2 = tdc(time, TS2[[ij]]$lPredY[, 2]), lY3 = tdc(time, TS2[[ij]]$lPredY[, 3]),
                                 lY4 = tdc(time, TS2[[ij]]$lPredY[, 4])
  )
  
  # Fit a Cox proportional hazards model using the joint model's longitudinal predictions
  cox_model_tsjm2[[ij]] <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
                                 data = long.data1, id = Id
  )
  
  
  # Prepare data for multiple imputation
  surv_data_new <- survival::tmerge(dataSurv_t, dataSurv_t, id = Id, endpt = event(deathTimes, Event))
  
  
  Limp <- 100
  SD2 <- matrix(0, Limp, length(cox_model_tsjm2[[ij]]$coefficients))
  TEM <- matrix(0, Limp, length(cox_model_tsjm2[[ij]]$coefficients))
  
  # Perform multiple imputation to estimate variability of coefficients
  for (k in 1:Limp) {
    lPredY <- TS2[[ij]]$LPredY[[k]]
    long.data1 <- survival::tmerge(surv_data_new, dataLong_t,
                                   id = Id,
                                   lY1 = tdc(time, TS2[[ij]]$lPredY[, 1]),
                                   lY2 = tdc(time, TS2[[ij]]$lPredY[, 2]),
                                   lY3 = tdc(time, TS2[[ij]]$lPredY[, 3]),
                                   lY4 = tdc(time, TS2[[ij]]$lPredY[, 4])
    )
    
    cox_model <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3 + lY4,
                       data = long.data1, id = Id
    )
    
    TEM[k, ] <- summary(cox_model)$coefficients[, 1]
    SD2[k, ] <- summary(cox_model)$coefficients[, 3]
  }
  
  # Compute the within-imputation variance (Wv) and between-imputation variance (Bv)
  Wv <- apply(SD2^2, 2, mean)
  Bv <- apply(TEM, 2, var)
  sdnew <- sqrt(Wv + (Limp + 1) / Limp * Bv)
  
  # Combine results into a single summary table
  Res <- cbind(
    cox_model_tsjm2[[ij]]$coefficients, sdnew,
    cox_model_tsjm2[[ij]]$coefficients - qnorm(.95) * sdnew,
    cox_model_tsjm2[[ij]]$coefficients + qnorm(.95) * sdnew
  )
  
  colnames(Res) <- c("coefficients", "sd", "L_CI", "U_CI")
  print(Res)
  
  end2 <- Sys.time()
  Timetsjm2[ij] <- difftime(end2, start2, units = "mins")
  TS2
  ##################################################
  
 
  ##################################################
  print(ij)
}



mean(Timetsjm0);sd(Timetsjm0); 
mean(Timetsjm1);sd(Timetsjm1); 
mean(Timetsjm2);sd(Timetsjm2); 
```

