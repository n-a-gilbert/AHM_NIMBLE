# 2.9.2 Multistate Dail-Madsen models

library(AHMbook)
library(nimble)
library(parallel)

data("duskySalamanders")

n <- duskySalamanders

nsites <- dim(n)[1]
nyears <- dim(n)[2]
nObsAges <- dim(n)[3]
nsurveys <- dim(n)[4]

data <- list(
  n = n
)

constants <- list(
  nsites = nsites, 
  nyears = nyears, 
  nsurveys = nsurveys
)

msdm <- nimbleCode({
  
  # priors
  for(c in 1:3){
    lambda[c] ~ dexp(0.1)
    gamma[c]  ~ dexp(0.1)
  }
  
  phi[1] ~ dbeta(1, 1)
  phi[2] ~ dbeta(1, 1)
  p[1]   ~ dbeta(1, 1)
  p[2]   ~ dbeta(1, 1)
  
  # process model
  for(i in 1:nsites){
    
    # Initial abundance
    N[i, 1, 1] ~ dpois(lambda[1]) # year 1 juveniles
    N[i, 1, 2] ~ dpois(lambda[2]) # year 2 juveniles
    N[i, 1, 3] ~ dpois(lambda[3]) # adults
    
    for(t in 2:nyears){
      # number of survivors in each age class
      S[i, t, 1] ~ dbin(phi[1], N[i, t-1, 1])
      S[i, t, 2] ~ dbin(phi[1], N[i, t-1, 2])
      S[i, t, 3] ~ dbin(phi[2], N[i, t-1, 3])
      # number of recruits and immigrants
      G[i, t, 1] ~ dpois(gamma[1]*N[i, t-1, 3] + gamma[2])
      G[i, t, 2] ~ dpois(gamma[2])
      G[i, t, 3] ~ dpois(gamma[3])
      # sum all stages to get total N at each site i in year t
      N[i, t, 1] <- G[i, t, 1]
      N[i, t, 2] <- S[i, t, 1] + G[i, t, 2]
      N[i, t, 3] <- S[i, t, 2] + S[i, t, 3] + G[i, t, 3]
    }
    
    # observation model
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        # Stages S1 and S2 are indistinguisable
        n[i, t, 1, j] ~ dbin(p[1], (N[i, t, 1] + N[i, t, 2]))
        # detection probability is the same for adults
        n[i, t, 2, j] ~ dbin(p[2], N[i, t, 3])
      }
    }
  }
  
  
  # derived quantities: total N for each stage and year
  for(t in 1:nyears){
    Ntotal[1, t] <- sum(N[1:nsites, t, 1])
    Ntotal[2, t] <- sum(N[1:nsites, t, 2])
    Ntotal[3, t] <- sum(N[1:nsites, t, 3])
  }
  
  
})

lamNew <- NA
phiNew <- NA
gammaNew <- NA
pNew <- NA

lamNew[1] <- 50
lamNew[2] <- 50
phiNew[1] <- 0.9
phiNew[2] <- 0.9
gammaNew[1] <- 50
gammaNew[2] <- 50
gammaNew[3] <- 50
pNew[1] <- 0.8
pNew[2] <- 0.8

N <- N1 <- array(NA, dim = c(nsites, nyears, 3) )
for(i in 1:nsites){
  for(t in 1:nyears){
    max.it1 <- max(n[i, t, 1, ]) + 5
    max.it2 <- max(n[i, t, 2, ]) + 5
    N[i, t, 1] <- max.it1
    N[i, t, 2] <- max.it1
    N[i, t, 3] <- max.it2 + 2
  }
}

N1[,1,] <- N[,1,]

inits <- function(){
  list(
    phi = phiNew, 
    gamma = rexp(3, 0.1), 
    lambda = rexp(3, 0.1),
    p = pNew, 
    N = N,
    S = array(sample(1:3), dim = c(nsites, nyears, 3)),
    G = array(sample(1:3), dim = c(nsites, nyears, 3))
  )
}


keepers <- c("lambda", "phi", "gamma", "p", "Ntotal")

ni <- 100000
nb <- 50000
nt <- 50
nc <- 3

# m <- nimbleModel(code = msdm,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# 52 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("msdm",
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
  
  m <- nimbleModel(code = msdm,
                   name = "msdm",
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

# catastrophic failure - bigtime convergence probs
MCMCvis::MCMCsummary(out, params = "lambda")
