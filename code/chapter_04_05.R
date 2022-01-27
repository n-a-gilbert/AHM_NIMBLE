# 4.5 time-dependent dynocc

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

set.seed(1)

d <- AHMbook::simDynocc(
  nsites = 250, 
  nyears = 10, 
  nsurveys = 3, 
  mean.psi1 = 0.6, 
  range.phi = c(0.5, 1), 
  range.gamma = c(0, 0.5), 
  range.p = c(0.1, 0.5)
)

data <- list(
  y = d$y
)

constants <- list(
  nsites = d$nsites, 
  nsurveys = d$nsurveys, 
  nyears = d$nyears
)

dynocc2 <- nimble::nimbleCode({
  
  # priors
  psi1 ~ dbeta(1, 1)
  for(t in 1:(nyears - 1)){
    phi[t] ~ dbeta(1, 1)
    gamma[t] ~ dbeta(1, 1)
    p[t] ~ dbeta(1, 1)
  }
  
  p[nyears] ~ dbeta(1, 1)
  
  # ecological submodel
  for(i in 1:nsites){
    z[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      z[i, t] ~ dbern(z[i, t-1]*phi[t-1] + (1 - z[i, t-1])*gamma[t-1])
    }
  }
  
  # observation model
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      for(t in 1:nyears){
        y[i, j, t] ~ dbern(z[i, t] * p[t])
      }
    }
  }
  
  # derived params
  lpsi1 <- logit(psi1)
  lp[1] <- logit(p[1])
  psi[1] <- psi1
  n.occ[1] <- sum(z[1:nsites, 1])
  for(t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1 - psi[t-1])*gamma[t-1]
    n.occ[t] <- sum(z[1:nsites, t])
    growthr[t-1] <- psi[t] / psi[t-1]
    turnover[t-1] <- (1 - psi[t-1])*gamma[t-1]/psi[t]
    lgamma[t-1] <- logit(gamma[t-1])
    leps[t-1] <- logit(1 - phi[t-1])
    lp[t] <- logit(p[t])
  }
})

zobs <- apply(data$y, c(1, 3), function(x) max(x, na.rm = TRUE))
n_na <- length(zobs[zobs == "-Inf"])
zobs[zobs == "-Inf"] <- rbinom(n_na, 1, 0.5)
zst <- zobs
inits <- function(){
  list(
    z = zst,
    psi1 = rbeta(1, 1, 1),
    phi = rbeta(constants$nyears-1, 1, 1),
    gamma = rbeta(constants$nyears-1, 1, 1),
    p = rbeta(constants$nyears, 1, 1)
  )
}

keepers <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover", "lpsi1", "lgamma", "leps", "lp")
ni <- 20000
nb <- 10000
nt <- 10
nc <- 3

# m <- nimbleModel(code = dynocc2,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 3.5 mins
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("dynocc2",
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
  
  m <- nimbleModel(code = dynocc2,
                   name = "dynocc2",
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

# probably could have run those chains a bit longer
# (some kind of high rhats, low n.eff)
MCMCvis::MCMCsummary(out)
