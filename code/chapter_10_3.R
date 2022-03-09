# Section 10.3: raw and aggregated occupancy data

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)
library(tidyverse)


nsites1 <- 267
nsites2 <- 2000
nsurveys <- 3
mean.occ <- 0.4
beta1 <- -3
mean.det <- 0.4
alpha2 <- -3

set.seed(1)

str(data1 <- simOcc(
  M = nsites1, 
  J = nsurveys, 
  mean.occ = mean.occ, 
  beta1 = beta1, 
  beta2 = 0, 
  beta3 = 0, 
  mean.det = mean.det, 
  time.effects = c(0, 0), 
  alpha1 = 0, 
  alpha2 = alpha2, 
  alpha3 = 0, 
  sd.lp = 0, b = 0
))

set.seed(24)
str(data2 <- simOcc(
  M = nsites2,
  J = nsurveys, 
  mean.occ = mean.occ, 
  beta1 = beta1, 
  beta2 = 0, 
  beta3 = 0, 
  mean.det = mean.det, 
  time.effects = c(0, 0), 
  alpha1 = 0, 
  alpha2 = alpha2, 
  alpha3 = 0, 
  sd.lp = 0,
  b = 0
))

head(y1 <- data1$y)
head(y2agg <- apply(data2$y, 1, max))
elev1 <- data1$elev
elev2 <- data2$elev
wind1 <- data1$wind

# 10.3.1 - Standard occupancy model to fit raw occupancy data

str(data <- list(
  y1 = y1, 
  elev1 = elev1, 
  wind1 = wind1))
))

str(constants <- list(
  nsite1 = nrow(y1), 
  nrep1 = ncol(y1)
))

model1 <- nimble::nimbleCode({
  
  # Priors
  alpha.lpsi <- logit(mean.psi)
  mean.psi ~ dbeta(1, 1)
  beta.lpsi ~ dnorm(0, 1.5)
  alpha.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  beta.lp ~ dnorm(0, 1.5)
  
  # Likelihood dataset 1
  for(i in 1:nsite1){
    z1[i] ~ dbern(psi1[i])
    logit(psi1[i]) <- alpha.lpsi + beta.lpsi*elev1[i]
    for(j in 1:nrep1){
      y1[i, j] ~ dbern(z1[i]*p1[i, j])
      logit(p1[i, j]) <- alpha.lp + beta.lp*wind1[i, j]
    }
  }
  
  # Derived quantity dataset 1
  Nocc1 <- sum(z1[1:nsite1])
  
})

# initial values
zst1 <- apply(y1, 1, max)
inits <- function(){list(
  z1 = zst1, 
  mean.psi = rbeta(1, 1, 1), 
  beta.lpsi = rnorm(1, 0, 1.5), 
  mean.p = rbeta(1, 1, 1), 
  beta.lp = rnorm(1, 0, 1.5)
)}

keepers <- c("mean.psi", 
            "alpha.lpsi", 
            "beta.lpsi", 
            "mean.p", 
            "alpha.lp", 
            "beta.lp", 
            "Nocc1")

ni <- 4000
nt <- 2
nb <- 2000
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = model1, 
#                  constants = constants, 
#                  data = data, 
#                  inits = inits())

# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("model1",
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
  
  m <- nimbleModel(code = model1,
                   name = "model1", 
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


model1_summary <- MCMCsummary(out)

# 10.3.2 - Integrated model

data <- list(
  y1 = y1, 
  elev1 = elev1, 
  wind1 = wind1, 
  y2agg = y2agg, 
  elev2 = elev2
)

constants <- list(
  nsite1 = nrow(y1), 
  nrep1 = ncol(y1),
  nsite2 = length(y2agg)
)

model2 <- nimbleCode({
  
  # Priors
  alpha.lpsi <- logit(mean.psi) # occupancy intercecpt
  mean.psi ~ dbeta(1, 1)
  beta.lpsi ~ dnorm(0, 1.5)
  alpha.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1) 
  beta.lp ~ dnorm(0, 1.5)
  
  # likelihood for dataset 1 (raw data)
  for(i in 1:nsite1){
    z1[i] ~ dbern(psi1[i])
    # identical parameters in state model for dataset 2
    logit(psi1[i]) <- alpha.lpsi + beta.lpsi*elev1[i]
    for(j in 1:nrep1){
      y1[i, j] ~ dbern(z1[i]*p1[i, j])
      logit(p1[i, j]) <- alpha.lp + beta.lp*wind1[i, j]
    }
  }
  
  # likelihood for dataset 2 (aggregated data)
  for(i in 1:nsite2){
    z2[i] ~ dbern(psi2[i])
    logit(psi2[i]) <- alpha.lpsi + beta.lpsi*elev2[i]
    y2agg[i] ~ dbern(z2[i]*Pstar2)
  }
  
  Pstar2 <- 1 - pow((1 - mean.p), 3) # prob detected at least once

  # derived quantities
  Nocc1 <- sum(z1[1:nsite1])
  Nocc2 <- sum(z2[1:nsite2])

  })

zst1 <- apply(y1, 1, max)
zst2 <- rep(1, length(y2agg))
inits <- function(){list(
  z1 = zst1, 
  z2 = zst2,
  mean.psi = rbeta(1, 1, 1), 
  beta.lpsi = rnorm(1, 0, 1.5), 
  mean.p = rbeta(1, 1, 1), 
  beta.lp = rnorm(1, 0, 1.5)
)}

keepers <- c("mean.psi", 
             "alpha.lpsi", 
             "beta.lpsi", 
             "mean.p", 
             "alpha.lp", 
             "beta.lp", 
             "Nocc1",
             "Nocc2")

ni <- 4000
nt <- 2
nb <- 2000
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = model2,
#                  constants = constants,
#                  data = data,
#                  inits = inits())

# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("model2",
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
  
  m <- nimbleModel(code = model2,
                   name = "model2", 
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

model2_summary <- MCMCsummary(out2)

model2_summary %>% 
  tibble::as_tibble(rownames = "parameter") %>% 
  dplyr::filter(parameter == "Nocc1" | parameter == "Nocc2") %>% 
  tibble::add_column(truth = c(data1$sumZ, data2$sumZ),
                     observed = c(data1$sumZ.obs, data2$sumZ.obs)) %>% 
  dplyr::select(parameter, truth, observed, mean, lower = `2.5%`, upper = `97.5%`)

im_results <- model2_summary %>% 
  tibble::as_tibble(rownames = "parameter") %>% 
  dplyr::filter(!parameter %in% c("Nocc1", "Nocc2")) %>% 
  tibble::add_column(model = "int") %>% 
  dplyr::select(model, parameter, mean, lower = `2.5%`, upper = `97.5%`)

model1_summary %>% 
  tibble::as_tibble(rownames = "parameter") %>% 
  dplyr::filter(!parameter %in% c("Nocc1", "Nocc2")) %>% 
  tibble::add_column(model = "occ") %>% 
  dplyr::select(model, parameter, mean, lower = `2.5%`, upper = `97.5%`) %>% 
  dplyr::full_join(im_results) %>% 
  tidyr::pivot_wider(names_from = model, values_from = mean:upper) %>% 
  dplyr::filter(!grepl("alpha", parameter)) %>% 
  tibble::add_column(truth = c(alpha2,
                       beta1, 
                       mean.det, 
                       mean.occ)) %>% 
  dplyr::select(parameter, truth, mean_occ, lower_occ, upper_occ, mean_int, lower_int, upper_int)

  
