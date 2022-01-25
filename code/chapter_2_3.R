library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

data("greenWoodpecker")

d <- greenWoodpecker

counts <- as.matrix(d[,7:48])
dateso <- as.matrix(d[,49:90])
timeso <- as.matrix(d[, 91:132])
into <- timeso / d[,'route.length']


dates <- standardize(dateso)
int <- standardize(into)
times <- standardize(timeso)

C <- array(counts, dim = c(267, 3, 14))
DATE <- array(dates, dim = c(267, 3, 14))
DUR <- array(times, dim = c(267, 3, 14))
INT <- array(int, dim = c(267, 3, 14))

elev <- standardize(d$elevation)
forest <- standardize(d$forest)

DATE[is.na(DATE)] <- 0
INT[is.na(INT)] <- 0

# 2.3 YEAR-STRATIFIED N-MIXTURE MODEL
data <- list(
  C = C
)

constants <- list(
  nsites = dim(C)[1],
  nsurveys = dim(C)[2],
  nyears = dim(C)[3]
)

nmix1 <- nimbleCode({
  
  for(t in 1:nyears){
    lambda[t] ~ dgamma(1, 0.05)
    p[t] ~ dbeta(1, 1)
  }
  
  # ecological model for abundance
  for(i in 1:nsites){
    for(t in 1:nyears){
      N[i, t] ~ dpois(lambda[t])
      # observation model for replicated counts
      for(j in 1:nsurveys){
        C[i, j, t] ~ dbin(p[t], N[i, t])
      }
    }
  }
  
  # derived quantity: total abundance
  for(t in 1:nyears){
    totalN[t] <- sum(N[1:nsites, t])
  }
  
})

Nst <- apply(C, c(1, 3), max, na.rm = TRUE) + 1
Nst[Nst == "-Inf"] <- 1
inits <- function(){
  list(
    N = Nst, 
    lambda = runif(dim(C)[3]),
    p = rbeta(constants$nyears, 1, 1)
  )
}

keepers <- c("lambda", "p", "totalN")

ni <- 5000
nb <- 2500
nt <- 2
nc <- 3

# m <- nimbleModel(code = nmix1,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 3.4 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("nmix1",
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
  
  m <- nimbleModel(code = nmix1,
                   name = "nmix1",
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

# 2.3.1 - adding covariates and estimating a common trend over time

data <- list(
  C = C, 
  elev = elev, 
  forest = forest, 
  DATE = DATE,
  INT = INT
)

constants <- list(
  nsites = dim(C)[1], 
  nsurveys = dim(C)[2],
  nyears = dim(C)[3]
)

nmix2 <- nimbleCode({
  
  # priors
  for(t in 1:nyears){
    alpha0[t] <- logit(mean.p[t])
    mean.p[t] ~ dbeta(1, 1)
  }
  
  for(i in 1:2){
    alpha[i] ~ dnorm(0, sd = 0.5)
  }
  alpha[3]   ~ dlnorm(0, 1)        # increasing itensity of effort can only have a positive effect
  
  beta0 ~ dnorm(0, 2)
  
  for(i in 1:4){
    beta[i] ~ dnorm(0, sd = 0.5)
  }
  for(i in 1:nsites){
    for(t in 1:nyears){
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- beta0 + beta[1]*elev[i] + beta[2]*pow(elev[i], 2) + beta[3]*forest[i] + beta[4]*(t - 7.5)
      # observation model for replicated counts
      for(j in 1:nsurveys){
        C[i, j, t] ~ dbin(p[i, j, t], N[i, t])
        logit(p[i, j, t]) <- alpha0[t] + alpha[1]*DATE[i, j, t] + alpha[2]*pow(DATE[i, j, t], 2) + alpha[3]*INT[i, j, t]
      }
    }
  }
  
  # Derived quantity: total abundance across all surveyed sites
  for(t in 1:nyears){
    totalN[t] <- sum(N[1:nsites, t])
  }
  
})


Nst <- apply(C, c(1, 3), max, na.rm = TRUE) + 1
Nst[Nst == "-Inf"] <- 1
inits <- function(){
  list(
    N = Nst,
    alpha = c(rnorm(2, 0, 0.5), rlnorm(1, 0, 1)),
    mean.p = rbeta(constants$nyears, 1, 1),
    beta0 = rnorm(1, 0, 2),
    beta = rnorm(4, 0, 0.5)
  )
}

keepers <- c("alpha0", "alpha", "beta0", "beta","totalN")

ni <- 5000
nb <- 2500
nt <- 2
nc <- 3

m <- nimbleModel(code = nmix2,
                 constants = constants,
                 data = data,
                 inits = inits())
# 
m$initializeInfo()

# ART 5 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("nmix2",
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
  
  m <- nimbleModel(code = nmix2,
                   name = "nmix2",
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
