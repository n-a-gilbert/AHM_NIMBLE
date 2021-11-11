# Chapter 9.4.1
# CAR N-mixture model
# pages 540-545 in AHM2
# translated to NIMBLE
# 11 November 2021

library(AHMbook)
library(spdep)
library(nimble)
library(parallel)
library(MCMCvis)

data("BerneseOberland")
head(bo <- BerneseOberland)

set.seed(10)
s <- simExpCorrRF(theta = 10, size = 50)

nsites <- 2500
nreps <- 3

elev <- standardize(bo$elevation)
forest <- standardize(bo$forest)

# ecological process
beta0 <- 2
beta1 <- 2
beta2 <- -2
loglam0 <- beta0 + beta1*elev + beta2*elev*elev
loglam <- beta0 + beta1*elev + beta2*elev*elev + c(s$field)
lam0 <- exp(loglam0)
lam <- exp(loglam)

N <- rpois(n = nsites, lambda = lam)
table(N)
sum(N > 0) / nsites
(totalN <- sum(N))

wind <- matrix(rnorm(nsites*nreps), nrow = nsites, ncol = nreps)

#observation 
alpha0 <- 0
alpha1 <- -1
alpha2 <- -1

p <- array(NA, dim = c(nsites, nreps))
for(j in 1:nreps){
  p[, j] <- plogis(alpha0 + alpha1*forest + alpha2*wind[,j])
}

# count things
y <- array(dim = c(nsites, nreps))
for(j in 1:nreps){
  y[, j] <- rbinom(n = nsites, size = N, prob = p[, j])
}

set.seed(100, sample.kind = "Rounding")
sample.size <- 500
sample.sites <- sort(sample(1:nsites, size = sample.size))

yobs <- y

yobs[-sample.sites, ] <- NA
head(sample.sites)
head(yobs)

coordgrid <- cbind(bo$x, bo$y)
neigh <- spdep::dnearneigh(coordgrid, d1 = 0, d2 = sqrt(2)*1000 + 1)
winnb <- spdep::nb2WB(neigh)

# For NIMBLE, separate data and constants (used for indexing)
data <- list(y       = yobs, 
             adj     = winnb$adj, 
             weights = winnb$weights, 
             num     = winnb$num, 
             elev    = as.numeric(elev),
             elev2   = as.numeric(elev*elev),
             forest  = as.numeric(forest), 
             wind    = wind)

constants <- list(nsites = dim(data$y)[1], 
                  nrep   = dim(data$y)[2],
                  neigh  = length(data$adj))

# code is slightly modified from book, mostly just prior choices
# note, I wonder if beta0 should be dropped for this thing to converge?
car_nmix <- nimble::nimbleCode({
  
  # Priors
  # Using weakly information normal priors for intercepts and coefficients
  # Slightly wider prior for intercepts
  beta0  ~ dnorm(0, sd = 4) 
  alpha0 ~ dnorm(0, sd = 4)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, sd = 2)
    beta[v]  ~ dnorm(0, sd = 2)
  }
  
  # CAR prior for spatial random effects
  # adding zero_mean = 1 improves mixing :)
  eta[1:nsites] ~ dcar_normal(adj[1:neigh], weights[1:neigh], num[1:nsites], tau, zero_mean = 1)
  
  tau ~ dexp(1)
  
  # Model for abundance
  for(i in 1:nsites){
    log(lam[i]) <- beta0 + beta[1]*elev[i] + beta[2]*elev2[i] + eta[i]
    N[i] ~ dpois(lam[i])
  }
  
  # Observation model
  for(i in 1:nsites){
    for(j in 1:nrep){
      y[i, j] ~ dbin(p[i, j], N[i])
      logit(p[i, j]) <- alpha0 + alpha[1]*forest[i] + alpha[2]*wind[i, j]
    }
  }
  
  # Derived parameter: total abundance of the grid
  Ntotal <- sum(N[1:nsites])
  
})

# Initial values
Nst <- apply(yobs, 1, max)
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2

inits <- function(){
  list(
    N = Nst, 
    beta0 = rnorm(1, 0, 3), 
    alpha0 = rnorm(1, 0, 3), 
    alpha = rnorm(2, 0, 1.5), 
    beta = rnorm(2, 0, 1.5), 
    eta = rep(0, length(data$num)),
    tau = rexp(1, 1)
  )
}

keepers <- c("beta0", "beta", "alpha0", "alpha", "Ntotal", "tau")

# MCMC settings
# With these settings, it ran for ~ 1 hours (parallelized)
nc <- 3
ni <- 30000
nb <- 20000
nt <- 2  

# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("car_nmix",
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
  
  m <- nimbleModel(code = car_nmix,
                   name = "car_nmix", 
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

# Use the excellent MCMCvis package to summarise/visualize results
# Looking fairly good except for the abundance intercept,
# which is not quite converged and moreover pretty far from truth
# similar problem is noted in the text (and authors used much longer chains)
# Possible solution might be to simply drop the intercept and allow the 
# spatial random effect to act as the intercept?
MCMCvis::MCMCsummary(out)

# traceplots
MCMCtrace(out, 
          params = c("beta0", "alpha0", "beta", "alpha"),
          pdf = FALSE,
          type = "trace")

# visualization of estimates
# looking more or less good, except for abundance intercept :(
MCMCplot(out, 
         params = c("beta0", "alpha0", "beta", "alpha"), 
         ci = c(50, 95))
