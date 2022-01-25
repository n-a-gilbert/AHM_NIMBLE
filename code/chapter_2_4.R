# 2.4 Modeling temporary emigration 
# pp 85-89

library(AHMbook)
library(nimble)
library(MCMCvis)
library(parallel)

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

data <- list(
  C = C, 
  DATE = DATE, 
  INT = INT
)

constants <- list(
  nsites = dim(C)[1], 
  nsurveys = dim(C)[2], 
  nyears = dim(C)[3]
)

nmix3 <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(1, 0.05)
  theta  ~ dbeta(1, 1)
  
  for(t in 1:nyears){
    alpha0[t] <- logit(mean.p[t])
    mean.p[t] ~ dbeta(1, 1)
  }
  
  for(k in 1:2){
    alpha[k] ~ dnorm(0, sd = 0.5)
  }
  
  alpha[3] ~ dlnorm(0, 1)
  
  # ecological model
  for(i in 1:nsites){
    M[i] ~ dpois(lambda) # superpopulation size M
    for(t in 1:nyears){
      N[i, t] ~ dbin(theta, M[i])
      for(j in 1:nsurveys){
        C[i, j, t] ~ dbin(p[i, j, t], N[i, t])
        logit(p[i, j, t]) <- alpha0[t] + alpha[1]*DATE[i, j, t] + alpha[2]*pow(DATE[i, j, t], 2) + alpha[3]*INT[i, j, t]
      }
    }
  }
  
  # derived quantity: Total M and total N across surveyed sites
  totalM <- sum(M[1:nsites])
  for(t in 1:nyears){
    totalN[t] <- sum(N[1:nsites, t])
  }
  
})


Nst <- apply(C, c(1, 3), max, na.rm = TRUE) + 1
Nst[Nst == "-Inf"] <- 1
Mst <- apply(Nst, 1, max)
inits <- function(){
  list(
    N = Nst,
    M = Mst, 
    lambda = runif(1, 0, 1), 
    theta = rbeta(1, 1, 1),
    alpha = c(rnorm(2, 0, 0.5), rlnorm(1, 0, 1)),
    mean.p = rbeta(constants$nyears, 1, 1)
    )
}

keepers <- c("alpha0", "alpha", "lambda", "theta", "totalM", "totalN")

ni <- 15000
nb <- 5000
nt <- 10
nc <- 3

# m <- nimbleModel(code = nmix3,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 12 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("nmix3",
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
  
  m <- nimbleModel(code = nmix3,
                   name = "nmix3",
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

# time-varying theta

nmix4 <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(1, 0.05)
  
  for(t in 1:nyears){
    theta[t] ~ dbeta(1, 1)
    alpha0[t] <- logit(mean.p[t])
    mean.p[t] ~ dbeta(1, 1)
  }
  
  for(k in 1:2){
    alpha[k] ~ dnorm(0, sd = 0.5)
  }
  
  alpha[3] ~ dlnorm(0, 1)
  
  # ecological model
  for(i in 1:nsites){
    M[i] ~ dpois(lambda)
    for(t in 1:nyears){
      N[i, t] ~ dbin(theta[t], M[i])
      for(j in 1:nsurveys){
        C[i, j, t] ~ dbin(p[i, j, t], N[i, t])
        logit(p[i, j, t]) <- alpha0[t] + alpha[1]*DATE[i, j, t] + alpha[2]*pow(DATE[i, j, t], 2) + alpha[3]*INT[i, j, t]
      }
    }
  }
  
  # derived quantities
  totalM <- sum(M[1:nsites])
  for(t in 1:nyears){
    totalN[t] <- sum(N[1:nsites, t])
  }
})



Nst <- apply(C, c(1, 3), max, na.rm = TRUE) + 1
Nst[Nst == "-Inf"] <- 1
Mst <- apply(Nst, 1, max)
inits <- function(){
  list(
    N = Nst,
    M = Mst, 
    lambda = runif(1, 0, 1), 
    theta = rbeta(1, 1, 1),
    alpha = c(rnorm(2, 0, 0.5), rlnorm(1, 0, 1)),
    mean.p = rbeta(constants$nyears, 1, 1)
  )
}

keepers <- c("alpha0", "alpha", "lambda", "theta", "totalM", "totalN")

ni <- 15000
nb <- 5000
nt <- 10
nc <- 3

# m <- nimbleModel(code = nmix4,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 12.5 mins
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("nmix4",
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
  
  m <- nimbleModel(code = nmix4,
                   name = "nmix4",
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
