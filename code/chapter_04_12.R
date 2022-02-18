# 4.12 a demographic dynamic occupancy model

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

simDemoDynocc(psi1 = 1)
simDemoDynocc(nsites = 1000)
simDemoDynocc(nyears = 100)
simDemoDynocc(nvisit = 20)
simDemoDynocc(range.phi = c(0.8, 0.8))
simDemoDynocc(range.phi = c(0.2, 0.3), range.r = c(0, 0.2))
simDemoDynocc(range.phi = c(0.8, 1), range.r = c(0.5, 0.7))
simDemoDynocc(nvisit = 1)
simDemoDynocc(range.p = c(1, 1))

set.seed(24)
d<- simDemoDynocc(
  psi1 = 0.6, 
  nsites = 100, 
  nyears = 20, 
  nvisit = 5, 
  range.phi = c(0.1, 0.9), 
  range.r = c(0, 0.5), 
  range.p = c(0.1, 0.9)
)

data <- list(
  y = d$y
)

constants <- list(
  nsites = d$nsites, 
  nyears = d$nyears, 
  nvisit = d$nvisits, 
  first = d$f
)

demo_dynocc_1 <- nimbleCode({
  
  # priors
  for(t in 1:(nyears - 1)){
    phi[t] ~ dbeta(1, 1)
    r[t] ~ dbeta(1, 1)
    p[t] ~ dbeta(1, 1)
  }
  
  # alive / dead process specified conditional on first detection
  for(i in 1:nsites){
    z[i, first[i]] ~ dbern(1) # condition on year for first detection
    for(t in (first[i] + 1):nyears){
      z[i, t] ~ dbern(z[i, t-1] * ( phi[t-1] +  (1 - phi[t-1]) * r[t-1]  ) + (1 - z[i, t-1]) * r[t-1])
    }
  }
  
  # observations conditional on alive/ dead process
  for(i in 1:nsites){
    for(t in (first[i] + 1):nyears){
      for(j in 1:nvisit){
        y[i, j, t] ~ dbern(z[i, t] * p[t-1])
      }
    }
  }
})

zst <- zinit(apply(data$y, c(1, 3), max))      
inits <- function(){
  list(
    z = zst, 
    phi = rbeta(constants$nyears - 1, 1, 1), 
    r = rbeta(constants$nyears - 1, 1, 1), 
    p = rbeta(constants$nyears - 1, 1, 1)
  )
}

keepers <- c("phi", "r", "p")

ni <- 6000
nb <- 3000
nt <- 3
nc <- 3

# ART 1.5 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("demo_dynocc_1",
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
  
  m <- nimbleModel(code = demo_dynocc_1,
                   name = "demo_dynocc_1",
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

MCMCvis::MCMCsummary(out, params = c("phi", "r", "p"))

# fit model for french peregrine dataset

data("FrenchPeregrines")

fp <- FrenchPeregrines

y <- as.matrix(fp[,4:56])
f <- apply(y, 1, function(x) min(which(x!=0)))
f[f == "Inf"] <- ncol(y)

data <- list(
  y = y)

constants <- list(
  nsites = nrow(y), 
  nyears = ncol(y), 
  first = f
)

demo_dynocc_2 <- nimble::nimbleCode({
  
  phi[1] ~ dbeta(1, 1)
  r[1]   ~ dbeta(1, 1)
  lphi[1] <- logit(phi[1])
  lr[1] <- logit(r[1])
  for(t in 2:(nyears - 1)){
    logit(phi[t]) <- lphi[t]
    logit(r[t]) <- lr[t]
  }
  
  # random walk smoother
  for(t in 2:(nyears-1)){
    lphi[t] ~ dnorm(lphi[t-1], sd.lphi)
    lr[t]  ~ dnorm(lr[t-1], sd.lr)
  }
  
  sd.lphi ~ dexp(1)
  sd.lr ~ dexp(1)
  
  # likelihood
  # alive/dead process specified conditional on first detection
  for(i in 1:nsites){
    y[i, first[i]] ~ dbern(1) # condition on year of first detection
    for(t in (first[i] + 1):nyears){
      y[i, t] ~ dbern(y[i, t-1] * (phi[t-1] + (1 - phi[t-1]) * r[t-1]) + (1 - y[i, t-1]) * r[t-1])
    }
  }
  
})

y[is.na(y)] <- rbinom(1, 1, 0.5)
yst <- y

inits <- function(){
  list(
    y = yst,
    lphi = rnorm(52, 0, 2), 
    lr = rnorm(52, 0, 2),
    sd.lphi = rexp(1, 1), 
    sd.lr = rexp(1, 1), 
    phi = rbeta(52, 1, 1), 
    r = rbeta(52, 1, 1)
  )
}

keepers <-  c("phi", "r", "sd.lphi", "sd.lr")

ni <- 12000
nb <- 6000
nt <- 10
nc <- 3


m <- nimble::nimbleModel(code = demo_dynocc_2,
                         constants = constants,
                         data = data,
                         inits = inits())

m$initializeInfo()

# ART 1.75 min
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("demo_dynocc_2",
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
  
  m <- nimbleModel(code = demo_dynocc_2,
                   name = "demo_dynocc_2",
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

MCMCvis::MCMCsummary(out)
