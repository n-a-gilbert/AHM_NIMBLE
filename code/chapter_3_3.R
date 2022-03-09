library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

# Section 3.3

# 3.3.1
nspec <- 50
mu.lphi <- 0
sigma.lphi <- 0.5
mu.lp <- -1
sigma.lp <- 0.5

set.seed(123)
lphi <- rnorm(n = nspec, mean = mu.lphi, sd = sigma.lphi)
lp <- rnorm(n = nspec, mean = mu.lp, sd = sigma.lp)

phi <- plogis(lphi)
p <- plogis(lp)
round(sort(phi), 2)
round(sort(p), 2)

ch <- array(NA, dim = c(100, 6, nspec))

# simulate data for 50 spp
for(s in 1:nspec){
  data <- simCJS(phi = phi[s], p = p[s], show.plot = F)
  ch[,,s] <- data$ch
}

ch
ch[1:5,,]
head(tmp <- apply(ch, c(1, 3), sum))

constants <- list(
  f = data$f, 
  n.ind = data$n.ind, 
  n.occ = data$n.occ, 
  nspec = dim(ch)[3])

data <- list(
  y = ch )

cjs3 <- nimbleCode({
  
  # Priors and hyperpriors
  # define phi and p as random effects from a prior distributio
  for(s in 1:nspec){
    phi[s] <- ilogit(lphi[s])
    lphi[s] ~ dnorm(mu.lphi, sd = sd.lphi)
    p[s] <- ilogit(lp[s])
    lp[s] ~ dnorm(mu.lp, sd = sd.lp)
  }
  
  # hyperpriors
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dbeta(1, 1)
  sd.lphi ~ dexp(1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lp ~ dexp(1)
  
  # 'likelihood'
  for(s in 1:nspec){
    for(i in 1:n.ind){
      # define the latent state at first capture
      z[i, f[i], s] <- 1
      for(t in (f[i] + 1):n.occ){
        # state process: the latent alive/dead state
        z[i, t, s] ~ dbern(z[i, t-1, s] * phi[s] ) # phi indexed by species now
        # observation process
        y[i, t, s] ~ dbern(z[i, t, s] * p[s]) # p also indexed by species
      }
    }
  }
})

zst <- ch
for(s in 1:50){
  zst[,,s] <- zinit(ch[,,s])
}

inits <- function(){
  list(
    z = zst, 
    mean.phi = rbeta(1, 1, 1), 
    lphi = rnorm(constants$nspec, 0, sd = 2),
    lp = rnorm(constants$nspec, 0, sd = 2),
    sd.lphi = rexp(1, 1), 
    mean.p = rbeta(1, 1, 1), 
    sd.lp = rexp(1, 1)
  )
}

keepers <- c("mean.phi", "mu.lphi", "sd.lphi", "mean.p", "mu.lp", "sd.lp", "phi", "p")

ni <- 10000
nb <- 5000
nt <- 5
nc <- 3

# m <- nimbleModel(code = cjs3,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ~ 8 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs3",
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
  
  m <- nimbleModel(code = cjs3,
                   name = "cjs3", 
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

# 3.3.3 Community CJS models for exploration of patterns among species

set.seed(24)
nspec <- 50
beija <- c(rep("yes", nspec/2), rep("no", nspec/2))
mass <- round(runif(nspec, 10, 100))
massC <- mass - mean(mass)

# get the design matric and inspect
( DM <- model.matrix(~ beija*massC-1-massC) )

# pick values for parameter vector PV
mu.alpha.nonbeija <- 0
mu.alpha.beija <- 0.7
beta.nonbeija <- 0.03
beta.beija <- 0.01
( PV <- c(mu.alpha.nonbeija, mu.alpha.beija, beta.nonbeija, beta.beija) )

# expected survival at link scale
( LinP <- DM %*% PV )

# add random species effects and pick value of hyerparameter sigma_phi
sigma.phi <- 0.4
alpha0 <- rnorm(50, 0, sigma.phi)
lphi <- LinP + alpha0

# inverse-logit transformation yields survival for every species
exp.phi <- plogis(LinP) # expected survival wo species RE
phi <- plogis(lphi) # realized survival (with species effects)

data.frame(beija, mass, "expected.survival" = round(exp.phi, 3), "realized.survival" = round(phi, 3))

mu.lp <- -1
sigma.lp <- 0.5
lp <- rnorm(n = nspec, mean = mu.lp, sd = sigma.lp)
p <- plogis(p)
round(p, 3)

ch <- array(NA, dim = c(100, 6, 50))
for(s in 1:nspec){
  data <- simCJS(phi = phi[s], p = p[s], show.plot = F)
  ch[,,s] <- data$ch
}

new.beija <- rep(c(1, 2), each = 25)
constants <- list(
  f = data$f, 
  n.ind = data$n.ind,
  n.occ = data$n.occ,
  n.spec = dim(ch)[3],
  beija = new.beija)

data <- list(
  y = ch,
  massC = massC)

cjs5 <- nimbleCode({
  
  # Priors and hyperprios
  # submodel for how species vary in the community
  for(s in 1:n.spec){
    phi[s] <- ilogit(lphi[s])
    lphi[s] ~ dnorm(mu.lphi[s], sd = sd.lphi)
    mu.lphi[s] <- alpha.lphi[beija[s]] + beta.lphi[beija[s]]*massC[s]
    p[s] <- ilogit(lp[s])
    lp[s] ~ dnorm(mu.lp, sd = sd.lp)
  }

  # submodel for community hyperparameters
  for(k in 1:2){
    alpha.lphi[k] <- logit(mean.phi[k])
    mean.phi[k] ~ dbeta(1, 1)
    beta.lphi[k] ~ dnorm(0, sd = 1)
  }
  
  sd.lphi ~ dexp(1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lp ~ dexp(1)
  
  # likelihood
  for(s in 1:n.spec){
    for(i in 1:n.ind) {
      # define latent state at first capture
      z[i, f[i], s] <- 1
      for(t in (f[i] + 1):n.occ){
        # state process - dead/alive
        z[i, t, s] ~ dbern(z[i, t - 1, s] * phi[s] )
        # observation process
        y[i, t, s] ~ dbern(z[i, t, s] * p[s])
      }
    }
  }
})

zst <- ch

for(s in 1:50){
  zst[,,s] <- zinit(ch[,,s])
}

inits <- function(){
  list(
    z = zst, 
    mean.phi = rbeta(2, 1, 1),
    beta.lphi = rnorm(2, 0, 1),
    sd.lphi = rexp(1, 1), 
    mean.p = rbeta(1, 1, 1),
    sd.lp = rexp(1, 1), 
    lphi = rnorm(constants$n.spec, 0, sd = 2),
    lp = rnorm(constants$n.spec, 0, sd = 2)
    
  )
}

keepers <- c("mean.phi", "alpha.lphi", "beta.lphi", "sd.lphi", "mean.p", "mu.lp", "sd.lp", "mu.lphi", "phi", "p")

ni <- 20000
nb <- 10000
nt <- 10
nc <- 3

# m <- nimbleModel(code = cjs5,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ~ 8 minutes
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs5",
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
  
  m <- nimbleModel(code = cjs5,
                   name = "cjs5", 
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