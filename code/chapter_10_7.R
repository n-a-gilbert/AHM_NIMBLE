# Section 10.7 - a simple IPM

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

data("greenWoodpecker")

C <- array(as.matrix(greenWoodpecker[,7:48]), dim = c(267, 3, 14))

data <- list(
  C = C
)

constants <- list(
  nsites = dim(C)[1],
  nsurveys = dim(C)[2], 
  nyears = dim(C)[3]
)

dm <- nimble::nimbleCode({
  
  # priors  
  lambda ~ dgamma(5, 2)
  phi ~ dbeta(1, 1)
  gamma ~ dgamma(3, 3)
  p ~ dbeta(1, 1)
  
  # likelihood for counts
  for(i in 1:nsites) {
    N[i, 1] ~ dpois(lambda)
    for(t in 1:(nyears-1)){
      S[i, t + 1] ~ dbin(phi, N[i, t])
      R[i, t + 1] ~ dpois(gamma)
      N[i, t + 1] <- S[i, t+1] + R[i, t+1]
    }
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i, j, t] ~ dbin(p, N[i, t])
      }
    }
  }
})

Rst <- apply(C, c(1, 3), max, na.rm = TRUE)
Rst[Rst == "-Inf"] <- 1
Rst[,1] <- NA
N1 <- apply(apply(C, c(1, 3), max, na.rm = TRUE), 1, max, na.rm = TRUE)
N1[N1 == "-Inf"] <- 2
Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N1 
Sst <- array(1, dim = dim(Rst))
inits <- function(){
  list(
    lambda = rgamma(1, 5, 2), 
    phi = rbeta(1, 1, 1), 
    gamma = rgamma(1, 3, 3), 
    p = rbeta(1, 1, 1), 
    R = Rst + 1, 
    S = Sst,
    N = Nst + 2
  )
}

keepers <- c("lambda", "phi", "gamma", "p")
ni <- 30000
nb <- 20000
nt <- 10
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = dm,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# #
# m$initializeInfo()

# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("dm",
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
  
  m <- nimbleModel(code = dm,
                   name = "dm", 
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

MCMCsummary(out)

set.seed(24)
data1 <- simCJS(
  n.occ = 6,
  n.marked = 20, 
  phi = 0.58,
  p = 0.4,
  show.plot = FALSE
  
)

( m_array1 <- ch2marray(data1$ch))
(r1 <- rowSums(m_array1))

data <- list(C = C,
             m_array = m_array1, 
             r = r1)

constants <- list(
  nsites = dim(C)[1],
  nsurveys = dim(C)[2], 
  nyears = dim(C)[3],
  n.occ = data1$n.occ
)

ipm <- nimble::nimbleCode({
  
  lambda ~ dgamma(2, 0.5)
  phiJoint ~ dbeta(1, 1)
  gamma ~ dgamma(3, 3)
  pDM ~ dbeta(1, 1)
  pCJS ~ dbeta(1, 1)
  
  # (1) Likelihood for the counts
  
  # make one parameter identical in the two submodels
  phiDM <- phiJoint
  
  for(i in 1:nsites){
    N[i, 1] ~ dpois(lambda)
    # state process - transition model
    for(t in 1:(nyears - 1)){
      S[i, t+1] ~ dbin(phiDM, N[i, t])
      R[i, t+1] ~ dpois(gamma)
      N[i, t+1] <- S[i, t+1] + R[i, t+1]
    }
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i, j, t] ~ dbin(phiDM, N[i, t])
      }
    }
  }
  
  # (2) Multinomial likelihood for the m-array format of the CR data 
  
  # constraints for CJS parameters
  for(t in 1:(n.occ-1)){
    phiCJS[t] <- phiJoint
    pCJS_t[t] <- pCJS
  }
  
  #multinomial likelhiood
  for(t in 1:(n.occ-1)){
    m_array[t, 1:n.occ] ~ dmulti(pr[t, 1:n.occ], r[t])
  }
  # define cell probabilities of the m-array
  for(t in 1:(n.occ-1)){
    q[t] <- 1 - pCJS_t[t] # probability of non-recacpture
    pr[t,t] <- phiCJS[t]*pCJS_t[t]
    for(j in (t + 1):(n.occ - 1)){
      pr[t, j] <- prod(phiCJS[t:j])*prod(q[t:(j-1)])*pCJS_t[t]
    }
    for(j in 1:(t-1)){
      pr[t,j] <- 0
    }
  }
  for(t in 1:(n.occ-1)){
    pr[t, n.occ] <- 1 - sum(pr[t,1:(n.occ-1)])
  }
  
})

Rst <- apply(C, c(1, 3), max, na.rm = TRUE)
Rst[Rst == "-Inf"] <- 1
Rst[,1] <- NA
N1 <- apply(apply(C, c(1, 3), max, na.rm = TRUE), 1, max, na.rm = TRUE)
N1[N1 == "-Inf"] <- 2
Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N1 
Sst <- array(1, dim = dim(Rst))
inits <- function(){
  list(
    lambda = rgamma(1, 5, 2), 
    phiJoint = rbeta(1, 1, 1), 
    gamma = rgamma(1, 3, 3), 
    pDM = rbeta(1, 1, 1),
    pCJS = rbeta(1, 1, 1),
    R = Rst + 1, 
    S = Sst,
    N = Nst + 2
  )
}

keepers <- c("lambda", "phiJoint", "gamma", "pDM", "pCJS")

ni <- 30000
nb <- 20000
nt <- 10
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = ipm,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# took ~ 20 mins with these settings
# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("ipm",
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
  
  m <- nimbleModel(code = ipm,
                   name = "ipm", 
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

MCMCsummary(out2)
