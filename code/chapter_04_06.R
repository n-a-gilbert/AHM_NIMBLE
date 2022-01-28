# 4.6 trend estimation with a dynocc model
# note, this is simplest trend possible (linear)

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

set.seed(24)

gam <- 0.1
eps <- 0.2

d <- AHMbook::simDynocc(
  mean.psi1 = 0.8, 
  range.phi = c(1 - eps, 1-eps), 
  range.gamma = c(gam, gam), 
  range.p = c(0.4, 0.4)
  )

data <- list(
  y = d$y
)

constants <- list(
  nsites = d$nsites, 
  nsurveys = d$nsurveys, 
  nyears = d$nyears
)

dynocc3 <- nimble::nimbleCode({
  
  # priors and model for occupancy
  for(t in 1:nyears){
    logit(psi[t]) <- alpha + beta.trend * (t - 5.5) + eps.year[t]
    p[t] ~ dbeta(1, 1)
    eps.year[t] ~ dnorm(0, sd.lpsi)
  }
  
  alpha <- logit(mean.psi)
  mean.psi ~ dbeta(1, 1)
  beta.trend ~ dnorm(0, 1)
  sd.lpsi ~ dexp(1)
  
  #ecological submodel
  for(i in 1:nsites){
    for(t in 1:nyears){
      z[i, t] ~ dbern(psi[t])
      # observation model
      for(j in 1:nsurveys){
        y[i, j, t] ~ dbern(z[i, t] * p[t])
      }
    }
  }
  
  # derived parameters
  for(t in 1:nyears){
    n.occ[t] <- sum(z[1:nsites, t]) # finite sample ocupancy
    logit(psi.trend[t]) <- alpha + beta.trend*(t - 5.5)
  }
})

inits <- function(){
  list(
    z = apply(data$y, c(1, 3), max),
    p = rbeta(constants$nyears, 1, 1), 
    eps.year = rnorm(constants$nyears, 0, 1), 
    mean.psi = rbeta(1, 1, 1),
    beta.trend = rnorm(1, 0, 1), 
    sd.lpsi = rexp(1, 1)
  )
}

keepers <- c("psi", "psi.trend", "mean.psi", "alpha", "beta.trend", "sd.lpsi", "p", "n.occ")  

ni <- 5000
nb <- 2500
nt <- 2
nc <- 3

# m <- nimbleModel(code = dynocc3,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 1.5 mins
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("dynocc3",
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
  
  m <- nimbleModel(code = dynocc3,
                   name = "dynocc3",
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
