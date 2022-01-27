# 4.3 Simulation and analysis of simplest possible dynocc

library(nimble)
library(parallel)
library(MCMCvis)

nsites <- 100
nyears <- 12
nsurveys <- 2
z <- array(NA, dim = c(nsites, nyears))
y <- array(NA, dim = c(nsites, nsurveys, nyears))

# set parameter values
psi1 <- 0.7
phi <- 0.9
gamma <- 0.05
p <- 0.25
( psi.eq <- gamma / (gamma + (1 - phi)))

set.seed(1)

# year 1
z[,1] <- rbinom(n = nsites, 
                size = 1, 
                prob = psi1)

# generate presence/absence in subsequent years
for(t in 2:nyears){
  exp.z <- z[,t-1]*phi + (1 - z[,t-1])*gamma
  z[,t] <- rbinom(n = nsites, 
                  size = 1, 
                  prob = exp.z)
}

# simulate measurement with false-negatives only
for(t in 1:nyears){
  for(j in 1:nsurveys){
    y[, j, t] <- rbinom(n = nsites, 
                        size = 1, 
                        prob = z[,t]*p)
  }
}

str(y)

# generate missing values

prob.missing <- 0.2
y2 <- y
for(i in 1:nsites){
  for(j in 1:nsurveys){
    for(t in 1:nyears){
      turnNA <- rbinom(1, 1, prob.missing)
      y2[i, j, t] <- ifelse(turnNA == 1, NA, y2[i, j, t])
    }
  }
}

str(y2)

table(nvisits <- apply(y2, c(1, 3), function(x) sum(!is.na(x))))

)
zobs <- apply(y2, c(1, 3), function(x) max(x, na.rm = TRUE))
n_na <- length(zobs[zobs == "-Inf"])
zobs[zobs == "-Inf"] <- rbinom(n_na, 1, 0.5)

data <- list(
  y = y2
)

constants <- list(
  nsites = dim(y2)[1],
  nsurveys = dim(y2)[2], 
  nyears = dim(y2)[3]
)

dynocc <- nimble::nimbleCode({
  
  # priors
  psi1  ~ dbeta(1, 1)
  phi   ~ dbeta(1, 1)
  gamma ~ dbeta(1, 1)
  p     ~ dbeta(1, 1)
  
  # ecological submodel
  for(i in 1:nsites){
    # initial state
    z[i, 1] ~ dbern(psi1)
    # state transitions
    for(t in 2:nyears){
      z[i, t] ~ dbern(z[i, t-1]*phi + (1 - z[i, t-1])*gamma)
    }
  }
  
  # observation model
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      for(t in 1:nyears){
        y[i, j, t] ~ dbern(z[i, t]*p)
      }
    }
  }
  
  # derived parameters
  eps <- 1 - phi # extinction probability
  psi[1] <- psi1 # population occupancy
  for(t in 2:nyears){
    psi[t] <- psi[t-1]*phi + (1 - psi[t-1])*gamma
  }
})

zst <- zobs
inits <- function(){
  list(
    z = zst, 
    psi1 = rbeta(1, 1, 1), 
    phi = rbeta(1, 1, 1), 
    gamma = rbeta(1, 1, 1), 
    p = rbeta(1, 1, 1)
  )
}

keepers <- c("psi1", "phi", "eps", "gamma", "p", "psi")

ni <- 15000
nb <- 10000
nt <- 5
nc <- 3

# m <- nimbleModel(code = dynocc,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 1.4 minus
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("dynocc",
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
  
  m <- nimbleModel(code = dynocc,
                   name = "dynocc",
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
