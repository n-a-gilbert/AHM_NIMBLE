# 10.4 counts + deteciton/nondetection

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

nsites1 <- 267       # number of sites with count data
nsites2 <- 520      # ... with detection/nondetection
nsurveys <- 3          
mean.lam <- 3
beta3.lam <- -1
mean.p <- 0.4
beta5.p <- 1         # effect of site covariate on p
beta.p.survey <- -2  # ...survey covariate on p

set.seed(1)
str(data1 <- simNmix(
  nsite = nsites1,
  nvisit = nsurveys, 
  mean.lam = mean.lam, 
  beta3.lam = beta3.lam, 
  mean.p = mean.p, 
  beta5.p = beta5.p, 
  beta.p.survey = beta.p.survey,
  show.plots = FALSE
))

# create dataset 2

set.seed(24)
str(data2 <- simNmix(
  nsite = nsites2, 
  nvisit = nsurveys, 
  mean.lam = mean.lam, 
  beta3.lam = beta3.lam, 
  mean.p = mean.p, 
  beta5.p = beta5.p, 
  beta.p.survey = beta.p.survey,
  show.plots = FALSE
))

head(y1 <- data2$C)

# degrade count data into deteciton/nondetection
head(y2 <- data2$C)
y2[y2 > 1] <- 1
head(y2)

elev1 <- data1$site.cov[,3]
elev2 <- data2$site.cov[,3]
hcov1 <- data1$site.cov[,5]
hcov2 <- data2$site.cov[,5]
wind1 <- data1$survey.cov
wind2 <- data2$survey.cov

# skipping fitting a N-mixture model to dataset1 only; moving on to integrated model

data <- list(
  y1 = data1$C, 
  y2 = y2, 
  elev1 = elev1, 
  elev2 = elev2, 
  hcov1 = hcov1, 
  hcov2 = hcov2, 
  wind1 = wind1, 
  wind2 = wind2
)

constants <- list(
  nsites1 = nsites1, 
  nsites2 = nsites2, 
  nsurveys = nsurveys
)

model <- nimble::nimbleCode({
  
  # Priors
  alpha.lam ~ dnorm(0, 2)
  beta.lam  ~ dnorm(0, 1)
  alpha.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  beta.lp1 ~ dnorm(0, 1)
  beta.lp2 ~ dnorm(0, 1)
  
  # Process model with shared parameters for datasets 1 and 2
  # Note that alpha.lam and beta.lam are shared
  
  for(i in 1:nsites1){
    N1[i] ~ dpois(lambda1[i])
    lambda1[i] <- exp(alpha.lam + beta.lam*elev1[i])
  }
  for(i in 1:nsites2){
    N2[i] ~ dpois(lambda2[i])
    lambda2[i] <- exp(alpha.lam + beta.lam*elev2[i])
  }

  # observational submodel for counts
  for(i in 1:nsites1){
    for(j in 1:nsurveys){
      y1[i, j] ~ dbin(p1[i, j], N1[i])
      logit(p1[i, j]) <- alpha.lp + beta.lp1*hcov1[i] + beta.lp2*wind1[i, j]
    }
  }
  
  # observational submodel for detection/nondetection
  for(i in 1:nsites2){
    for(j in 1:nsurveys){
      y2[i, j] ~ dbern(Pstar2[i, j])
      Pstar2[i, j] <- 1 - pow((1 - p2[i, j]), N2[i])
      logit(p2[i, j]) <- alpha.lp + beta.lp1*hcov2[i] + beta.lp2*wind2[i, j]
    }
  }
  
  # derived quantities
  Ntotal1 <- sum(N1[1:nsites1])
  Ntotal2 <- sum(N2[1:nsites2])
  
  
  })

# Initial values
Nst1 <- apply(data$y1, 1, max)
Nst2 <- rep(round(mean(data$y1)), nsites2)
inits <- function(){
  list(
    N1 = Nst1, 
    N2 = Nst2, 
    alpha.lam = rnorm(1, 0, 2), 
    beta.lam = rnorm(1, 0, 1), 
    mean.p = rbeta(1, 1, 1), 
    beta.lp1 = rnorm(1, 0, 1), 
    beta.lp2 = rnorm(1, 0, 1)
  )
}

keepers <- c("alpha.lam", "beta.lam", "alpha.lp", "mean.p", "beta.lp1", "beta.lp2", "Ntotal1", "Ntotal2")

ni <- 11000
nt <- 8
nb <- 3000
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = model,
#                  constants = constants,
#                  data = data,
#                  inits = inits())

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

library(tidyverse)
MCMCsummary(out) %>% 
  tibble::add_column(truth = c(data1$Ntotal, data2$Ntotal,
                               NA, NA, beta3.lam, beta5.p, beta.p.survey, mean.p))
