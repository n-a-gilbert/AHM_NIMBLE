# 10.5 Counts + binned counts

library(AHMbook)
library(berryFunctions) # wtf weird package name?
library(nimble)
library(parallel)
library(MCMCvis)
library(tidyverse)

M1 <- 500 # number of sites with abundance class surveys
J1 <- 5 # number of surveys in abundance class surveys
M2 <- 100 # number of sites with counts
J2 <- 2 # number of replicate surveys for counts
mean.lam <- 50 # abundance intercept
beta.lam <- 0.5 # coefficient of site covariate, abundance
mean.p1 <- 0.4 # detection intercept in abundance-class survey
mean.p2 <- 0.5 # detection intercept in full-count survey
beta.p <- -1 # coefficient for site covariate in detection

set.seed(1)

# simulate data for abundance class surveys
data1 <- simNmix(
  nsite = M1, 
  nvisit = J1, 
  mean.lam = mean.lam, 
  mean.p = mean.p1, 
  beta2.lam = beta.lam, 
  beta3.p = beta.p, 
  show.plots = FALSE
)

# simulate data for straight up counts
data2 <- simNmix(
  nsite = M2, 
  nvisit = J2,
  mean.lam = mean.lam, 
  mean.p = mean.p2, 
  beta2.lam = beta.lam, 
  beta3.p = beta.p, 
  show.plots = FALSE
)

breaks <- c(0, 10, 25, 50, 100, 200)
Cclass <- berryFunctions::classify(c(data1$C), 
                                   method = 'custom',
                                   breaks = breaks)$index
Aclass <- matrix(Cclass, byrow = FALSE, ncol = data1$nvisits)

table(Aclass)
breaks
n <- length(Cclass)
limits <- array(NA, c(n, 2), list(1:n, c("Lower", "Upper")))
for(i in 1:n){
  limits[i, 1:2] <- c(breaks[Cclass[i]], breaks[Cclass[i] + 1])
}

head(cbind("True" = c(data1$C), "Class" = c(Aclass), limits))

# vectorize the environmental covariate and get site index
X3vec <- rep(data1$site.cov[,3], 5)
sitevec <- rep(1:500, 5) # 500 sites with 5 reps each

# response is simply a vector of 1's!
y <- rep(1, 2500)

data <- list(
  y = y, 
  X2 = data1$site.cov[,2],
  limits = limits,
  X3vec = X3vec
)

constants <- list(
  M = nrow(Aclass), 
  n = length(Cclass),
  sitevec = sitevec
  )

model <- nimble::nimbleCode({
  
  #priors
  alpha.p <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  beta.p ~ dnorm(0, 1.5)
  alpha.lam ~ dnorm(0, 10)
  mean.lam <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 1.5)
  
  # likelihood
  # model for abundance
  for(i in 1:M){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha.lam + beta.lam*X2[i]
  }
  
  # observation model for observed counts and for detection
  for(i in 1:n){
    y[i] ~ dinterval(C[i], limits[i,1:2]) # specify interval censoring
    C[i] ~ dbin(p[i], N[sitevec[i]]) # count becomes estimated quantity
    logit(p[i]) <- alpha.p + beta.p*X3vec[i]
  }
  
  # derived quantity
  Ntotal <- sum(N[1:M])
})

tmp <- matrix(limits[,2], ncol = 5)
Nst <- round(apply(tmp, 1, mean)) + 20
Cst <- limits[,1] + 1
inits <- function(){
  list(
    N = Nst, 
    C = Cst, 
    mean.p = rbeta(1, 1, 1), 
    beta.p = rnorm(1, 0, 1), 
    alpha.lam = rnorm(1, 0, 10), 
    beta.lam = rnorm(1, 0, 1)
  )}

keepers <- c("alpha.lam", "beta.lam", "mean.lam", "alpha.p", "beta.p", "mean.p", "Ntotal", "N", "C")

# Took 15 minutes to run with these settings, not quite converged (some rhats ~ 1.2)
# probably would have to crank up to 100k iterations / more aggressive thinning to reach convergence
ni <- 50000
nb <- 25000
nt <- 10
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = model,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# m$initializeInfo()

# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("model",
                              "inits",
                              "data",
                              "constants",
                              "keepers",
                              "nb",
                              "ni",
                              "nt"))

for(j in seq_along(cl)){
  set.seed(j)
  init <- inits()
  clusterExport(cl[j], "init")
}

out <- clusterEvalQ(cl, {
  library(nimble)
  
  m <- nimbleModel(code = model,
                   name = "model", 
                   constants = constants,
                   data = data,
                   inits = init)
  
  Cmodel <- compileNimble(m)
  modelConf <- configureMCMC(m)
  modelConf$addMonitors(keepers)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = m)
  out1 <- runMCMC(CmodelMCMC,
                  nburnin = nb,
                  niter = ni,
                  thin = nt)
  
  
})

# turn off your cluster
stopCluster(cl)
end <- Sys.time()
end - start

MCMCvis::MCMCsummary(out, params = c("alpha.lam", "beta.lam", "mean.lam",
                                     "alpha.p", "beta.p", "mean.p", "Ntotal"))

n_estimated <- MCMCsummary(out, params = c("N"))

n_estimated %>% 
  add_column(truth = data1$N) %>% 
  ggplot() +
  geom_point(aes(x = truth, y = mean), alpha = 0.2) + 
  geom_abline(slope = 1, color = "red", size = 2) + 
  labs(x = "True N", y = "Estimated N")

data <- list(
  y = y, 
  X21 = data1$site.cov[,2],
  limits = limits,
  X3vec1 = X3vec,
  C2 = data2$C, 
  X22 = data2$site.cov[,2],
  X32 = data2$site.cov[,3]
)

constants <- list(
  M1 = nrow(Aclass), 
  n = length(Cclass),
  sitevec = sitevec, 
  M2 = nrow(data2$C),
  J2 = ncol(data2$C)
)

imodel <- nimble::nimbleCode({
  
  # Priors
  # for dataset 1 (binned counts)
  alpha.p1 <- logit(mean.p1)
  mean.p1 ~ dbeta(1, 1)
  # for dataset 2 (counts)
  alpha.p2 <- logit(mean.p2)
  mean.p2 ~ dbeta(1, 1)
  # Priors for shared parameters
  alpha.lam ~ dnorm(0, 10)
  mean.lam <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 2)
  beta.p ~ dnorm(0, 2)
  
  # Likelihood
  for(i in 1:M1){
    N1[i] ~ dpois(lambda1[i])
    log(lambda1[i]) <- alpha.lam + beta.lam*X21[i]
  }
  for(i in 1:M2){
    N2[i] ~ dpois(lambda2[i])
    log(lambda2[i]) <- alpha.lam + beta.lam*X22[i]
  }
  
  # Observation model for binned counts
  for(i in 1:n){
    y[i] ~ dinterval(C1[i], limits[i, 1:2])
    C1[i] ~ dbin(p1[i], N1[sitevec[i]])
    logit(p1[i]) <- alpha.p1 + beta.p*X3vec1[i]
  }
  
  # Observation model for dataset 2
  for(i in 1:M2){
    logit(p2[i]) <- alpha.p2 + beta.p*X32[i]
    for(j in 1:J2){
      C2[i, j] ~ dbin(p2[i], N2[i])
    }
  }
  
  # Derived quantities
  Ntotal1 <- sum(N1[])
  Ntotal2 <- sum(N2[])
  GTotalN <- Ntotal1 + Ntotal2
  
})

tmp <- matrix(limits[,2], ncol = 5)
Nst1 <- round(apply(tmp, 1, mean)) + 20
Cst1 <- limits[,1] + 1
Nst2 <- apply(data$C2, 1, max) + 10
inits <- function(){
  list(
    N1 = Nst1, 
    C1 = Cst, 
    N2 = Nst2, 
    mean.p1 = rbeta(1, 1, 1), 
    mean.p2 = rbeta(1, 1, 1), 
    alpha.lam = rnorm(1, 0, 1), 
    beta.lam = rnorm(1, 0, 1), 
    beta.p = rnorm(1, 0, 1)
  )
}

# took 17 minutes with these settings & mixing was fine
ni <- 50000
nb <- 25000
nt <- 10
nc <- 3

keepers <- c("mean.lam", "beta.lam", "mean.p1", "mean.p2", "beta.p",
             "Ntotal1", "Ntotal2", "GTotalN", "N1", "C1", "N2")

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = imodel,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# m$initializeInfo()

# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("imodel",
                              "inits",
                              "data",
                              "constants",
                              "keepers",
                              "nb",
                              "ni",
                              "nt"))

for(j in seq_along(cl)){
  set.seed(j)
  init <- inits()
  clusterExport(cl[j], "init")
}

out2 <- clusterEvalQ(cl, {
  library(nimble)
  
  m <- nimbleModel(code = imodel,
                   name = "imodel", 
                   constants = constants,
                   data = data,
                   inits = init)
  
  Cmodel <- compileNimble(m)
  modelConf <- configureMCMC(m)
  modelConf$addMonitors(keepers)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = m)
  out1 <- runMCMC(CmodelMCMC,
                  nburnin = nb,
                  niter = ni,
                  thin = nt)
  
  
})

# turn off your cluster
stopCluster(cl)
end <- Sys.time()
end - start

MCMCsummary(out2, 
            params = c("mean.lam", "beta.lam", "mean.p1", "mean.p2", "beta.p",
                       "Ntotal1", "Ntotal2", "GTotalN"))
