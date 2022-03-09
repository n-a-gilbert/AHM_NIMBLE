library(unmarked)
library(AHMbook)
library(spdep)
library(nimble)
# 
# str(dat <- simDynoccSpatial(beta.Xpsi1 = 1, 
#                             beta.Xphi = 1, 
#                             beta.Xgamma = 1, 
#                             beta.Xp = 1) )

str(dat <- simDynoccSpatial(side = 30, 
                            nyears = 10, 
                            nsurveys = 3, 
                            mean.psi1 = 0.1, 
                            range.phi = c(0.5, 0.5), 
                            range.gamma = c(0.2, 0.2), 
                            range.p = c(0.4, 0.4), 
                            beta.Xautolog = c(1, 1), 
                            seed.XAC = 1, 
                            seed = 24, 
                            show.plots = FALSE))

# summary(dat$umf)

# fm <- colext(~1, ~Xauto, ~Xauto, ~1, dat$umf)
# summary(fm)
# confint(fm, type = "col")[2,]

nsites <- dat$side^2
nsurveys <- dat$nsurveys
nyears <- dat$nyears
y <- array(NA, dim = c(nsites, nsurveys, nyears))
for(i in 1:nyears){
  y[,,i] <- dat$umf@y[,(3*i-2):(3*i)]
}

amat <- dat$amatrix
table(apply(amat, 2, sum))

neigh <- dnearneigh(dat$grid, d1 = 0, d2 = sqrt(2) + 0.1)
str(winnb <- nb2WB(neigh))
numN <- winnb$num

neighID <- array(0, dim = c(nsites, 8))
for(i in 1:nsites){
  neighID[i, 1:numN[i]] <- unlist(neigh[i])
}

data <- list(y = y
             )

constants <- list(
  nsites = nsites, 
  nsurveys = nsurveys, 
  nyears = nyears,
  neighID = neighID,
  numN = numN)

# can't index with a Z variable
# Workaround: https://groups.google.com/g/nimble-users/c/TUs5RCw193M

# calc_autcov <- nimbleFunction(
#   run = function(z = double(2), neighID = double(2), numN = double(1)){
#     returnType(double(2))
#     nsites <- dim(z)[1]
#     nyears <- dim(z)[2]
#     result <- matrix(nrow = nsites,
#                      ncol = nyears-1,
#                      init = FALSE)
#     for(i in 1:nsites){
#       for(t in 2:nyears){
#         result[i, t-1] <- sum(z[neighID[i, 1:numN[i]], t-1]) / numN[i]
#       }
#     }
#     return(result)
#   })

calc_autcov <- nimbleFunction(
  run = function(z = double(2), 
                 neighID = double(2),
                 numN = double(1),
                 i = integer(),
                 t = integer()){
    returnType(double(0))
    result <- sum(z[neighID[i, 1:numN[i]], t-1]) / numN[i]
    return(result)
  })



z <- array(rpois(nsites*nyears, 1.5), dim = c(nsites, nyears))
neighID <- data$neighID
numN <- data$numN

# works
result <- array(NA, dim = c(nsites, nyears - 1))
for(i in 1:nsites){
  for(t in 2:nyears){
    result[i, t-1] <- calc_autcov(z, neighID, numN, i, t)
  }
}

autologistic1 <- nimble::nimbleCode({
  
  psi1 ~ dbeta(1, 1)
  
  alpha.lphi ~ dnorm(0, sd = 3)
  beta.lphi  ~ dnorm(0, sd = 1.5)
  alpha.lgamma ~ dnorm(0, sd = 3)
  beta.lgamma ~ dnorm(0, sd = 1.5)
  p ~ dbeta(1, 1)
  
  
  # likelihood
  # ecological submodel
  for(i in 1:nsites){
    z[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      z[i, t] ~ dbern(z[i, t-1]*phi[i, t-1] + (1 - z[i, t-1])*gamma[i, t-1])
      # compute autocovariate and specificy its effects on phi and gamma
      autocov[i, t-1] <- calc_autcov(z[i, t-1],
                                     neighID[i, 1:numN[i]],
                                     numN[i], 
                                     i = i, 
                                     t = t)
      logit(phi[i, t-1]) <- alpha.lphi + beta.lphi*autocov[i, t-1]
      logit(gamma[i, t-1]) <- alpha.lgamma + beta.lgamma*autocov[i, t-1]
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
})

zst <- array(1, dim = c(nsites, nyears))

inits <- function(){
  list(
    z = zst, 
    psi1 = rbeta(1, 1, 1), 
    alpha.lphi = rnorm(1, 0, 1), 
    beta.lphi = rnorm(1, 0, 1), 
    alpha.lgamma = rnorm(1, 0, 1), 
    beta.lgamma = rnorm(1, 0, 1), 
    p = rbeta(1, 1, 1)
  )
}

keepers <- c("psi1", "alpha.lphi", "beta.lphi", "alpha.lgamma", "beta.lgamma", "p")
ni <- 6000
nt <- 3
nb <- 3000
nc <- 1

m <- nimbleModel(code = autologistic1, 
                 constants = constants, 
                 data = data, 
                 inits = inits())

Cmodel <- compileNimble(m)
modelConf <- configureMCMC(m)
modelConf$addMonitors(keepers)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = m)
out1 <- runMCMC(CmodelMCMC,
                nburnin = nb,
                niter = ni,
                thin = nt)


# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("splines_nmix",
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
  
  m <- nimbleModel(code = splines_nmix,
                   name = "splines_nmix", 
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

str(bdata <- list(y = y, nsites = nsites, nsurveys = nsurveys, nyears = nyears,
                  neighID = neighID, numN = numN))
# List of 6
# $ y       : int [1:900, 1:3, 1:10] 0 0 0 0 1 0 0 0 0 0 ...
# $ nsites  : num 900
# $ nsurveys: num 3
# $ nyears  : num 10
# $ neighID : int [1:900, 1:8] 2 1 2 3 4 5 6 7 8 9 ...
# $ numN    : int [1:900] 3 5 5 5 5 5 5 5 5 5 ...


# Specify model in BUGS language
cat(file = "autologistic1.txt","
model {
  # Priors
  psi1 ~ dunif(0, 1)                   # Initial occupancy
  phi.int ~ dunif(0, 1)                # Persistence
  alpha.lphi <- logit(phi.int)
  beta.lphi ~ dnorm(0, 0.01)
  gamma.int ~ dunif(0, 1)              # Colonization
  alpha.lgamma <- logit(gamma.int)
  beta.lgamma ~ dnorm(0, 0.01)
  p ~ dunif(0, 1)                      # Detection
  # Likelihood
  # Ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1] + (1-z[i,t-1]) * gamma[i,t-1])
      # Compute autocovariate and specify its effects on phi and gamma
      autocov[i,t-1] <- sum(z[neighID[i,1:numN[i]], t-1]) / numN[i]
      logit(phi[i,t-1]) <- alpha.lphi + beta.lphi * autocov[i,t-1]
      logit(gamma[i,t-1]) <- alpha.lgamma + beta.lgamma * autocov[i,t-1]
    }
  }
  # Observation model
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      for (t in 1:nyears){
        y[i,j,t] ~ dbern(z[i,t] * p)
      }
    }
  }
}
")

# Initial values
zst <- array(1, dim = c(nsites, nyears))
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi1", "phi.int", "alpha.lphi", "beta.lphi", "gamma.int",
            "alpha.lgamma", "beta.lgamma", "p") # could also monitor "z"

# MCMC settings
# na <- 1000 ; ni <- 6000 ; nt <- 3 ; nb <- 3000 ; nc <- 3
na <- 1000 ; ni <- 600 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~ for testing, 8 mins

# Call JAGS (ART 30 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "autologistic1.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

str(out1)
