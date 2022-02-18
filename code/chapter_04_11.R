# 4.11 accounting for preferential sampling

library(AHMbook)
library(nimble)
library(parallel)

data("FrenchPeregrines")

dat <- FrenchPeregrines

ain <- which(dat$department == "Ain")
jura <- which(dat$department == "Jura")
doubs <- which(dat$department == "Doubs")

ht <- as.numeric(dat$height)
y <- as.matrix(dat[,4:56])
nsites <- nrow(y)
nyears <- ncol(y)
year <- 1964:2016

( n.occ.obs <- apply(y, 2, sum, na.rm = TRUE) )
( n.visited <- apply(y, 2, function(x) sum(!is.na(x))) )
( n.ratio <- n.occ.obs / (n.visited / 284) )

# model 1 - no preferential sampling

data <- list(
  y = unname(y)
)

constants <- list(
  nsites = nsites, 
  nyears = nyears,
  height = ht#, 
  # ain = ain,
  # jura = jura,
  # doubs = doubs
)

dynocc1 <- nimble::nimbleCode({
  
  psi1 ~ dbeta(1, 1)
  # model for phi and gamma - cliff hiehgt + site + smooth random year effects
  for(i in 1:nsites){
    for(t in 1:(nyears - 1)){
      logit(phi[i, t]) <- lphi[i, t]
      lphi[i, t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i, t]) <- lgamma[i, t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], sd.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], sd.lgamma.site)
  }
  
  # priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  
  sd.lphi.site ~ dexp(1)
  sd.lgamma.site ~ dexp(1)
  
  # priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0
  lgamma.year[1] <- 0
  
  for(t in 2:(nyears - 1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], sd.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], sd.eps.lgamma )
  }
  
  sd.eps.lphi ~ dexp(1)
  sd.eps.lgamma ~ dexp(1)
  
  
  # ecological and observation submodels confounded (no p)
  
  for(i in 1:nsites){
    y[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      y[i, t] ~ dbern(y[i, t-1]*phi[i, t-1] + (1 - y[i, t-1])*gamma[i, t-1])
    }
  }
  
  # derived parameters
  psi[1] <- psi1
  n.occ[1] <- sum(y[1:nsites, 1])
  for(t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites, t])
  }
  # year-specific average values of phi and gamma
  for( t in 1:(nyears - 1)){
    mean.phi.year[t] <- mean(phi[1:nsites, t])
    mean.gamma.year[t] <- mean(gamma[1:nsites, t])
  }
  
  # average gamma nd phi per cliff hiehgt category in central year
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }
  
  
  # pop size in each department
  # UGH, vector indexing is a pain in Nimble!!!
  # COuld probably accomplish this with a nimble funtion but I'm so over this
  # for(t in 1:nyears){
  #   n.ain[t] <- sum(y[ain, t])
  #   n.jura[t] <- sum(y[jura, t])
  #   n.doubs[t] <- sum(y[doubs, t])
  # }
  
  
})

keepers <- c("psi1", "alpha.lphi", "initial.phi", 
             "alpha.lgamma", "initial.gamma", "lphi.year", 
             "lgamma.year", "mean.phi.year", "mean.gamma.year", 
             "sd.eps.lphi", "sd.eps.lgamma", "gamma.cliff", 
             "phi.cliff", "n.occ",
             # "n.ain",
             # "n.jura", "n.doubs",
             "y")

ni <- 20000
nb <- 10000
nt <- 10
nc <- 3

yst <- data$y
yst[is.na(yst)] <- rbinom(sum(is.na(yst)), 1, 0.5)

inits <- function(){
  list(
    psi1 = rbeta(1, 1, 1),
    lphi.site = rnorm(constants$nsites, 0, 1), 
    lgamma.site = rnorm(constants$nsites, 0, 1), 
    initial.phi = rbeta(3, 1, 1), 
    initial.gamma = rbeta(3, 1, 1), 
    sd.lphi.site = rexp(1, 1), 
    sd.lgamma.site = rexp(1, 1),
    lphi.year = c(0, rnorm(constants$nyears-2, 0, 1)),
    lgamma.year = c(0, rnorm(constants$nyears-2, 0, 1)),
    sd.eps.lphi = rexp(1, 1), 
    sd.eps.lgamma = rexp(1, 1), 
    y = yst
  )
}

m <- nimble::nimbleModel(code = dynocc1,
                         constants = constants,
                         data = data,
                         inits = inits())

m$initializeInfo()

# ART 
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("dynocc1",
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
  
  m <- nimbleModel(code = dynocc1,
                   name = "dynocc1",
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

# should run chains longer...
MCMCvis::MCMCsummary(out, params = c("psi1", "initial.phi", "initial.gamma"))


# Model 2 - modeling preferential sampling bu joint modeling of occupancy and visitation


# compute visitation data V
V <- y
V[V == 0] <- 1
V[is.na(V)] <- 0 # not visited

# check out visitation data
y[1:5, 1:5]
V[1:5, 1:5]


data <- list(
  y = unname(y),
  V = unname(V), 
  year = ((1664:2016) - 1990)/26)

constants <- list(
  nsites = nsites, 
  nyears = nyears,
  height = ht#, 
  # ain = ain,
  # jura = jura,
  # doubs = doubs
)

dynocc2 <- nimble::nimbleCode({
  
  psi1 ~ dbeta(1, 1)
  # model for phi and gamma - cliff hiehgt + site + smooth random year effects
  for(i in 1:nsites){
    for(t in 1:(nyears - 1)){
      logit(phi[i, t]) <- lphi[i, t]
      lphi[i, t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i, t]) <- lgamma[i, t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], sd.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], sd.lgamma.site)
  }
  
  # priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  
  sd.lphi.site ~ dexp(1)
  sd.lgamma.site ~ dexp(1)
  
  # priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0
  lgamma.year[1] <- 0
  
  for(t in 2:(nyears - 1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], sd.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], sd.eps.lgamma )
  }
  
  sd.eps.lphi ~ dexp(1)
  sd.eps.lgamma ~ dexp(1)
  
  
  # ecological and observation submodels confounded (no p)
  
  for(i in 1:nsites){
    y[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      y[i, t] ~ dbern(y[i, t-1]*phi[i, t-1] + (1 - y[i, t-1])*gamma[i, t-1])
    }
  }
  
  # derived parameters
  psi[1] <- psi1
  n.occ[1] <- sum(y[1:nsites, 1])
  for(t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites, t])
  }
  # year-specific average values of phi and gamma
  for( t in 1:(nyears - 1)){
    mean.phi.year[t] <- mean(phi[1:nsites, t])
    mean.gamma.year[t] <- mean(gamma[1:nsites, t])
  }
  
  # average gamma nd phi per cliff hiehgt category in central year
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }
  
  
  # pop size in each department
  # UGH, vector indexing is a pain in Nimble!!!
  # COuld probably accomplish this with a nimble funtion but I'm so over this
  # for(t in 1:nyears){
  #   n.ain[t] <- sum(y[ain, t])
  #   n.jura[t] <- sum(y[jura, t])
  #   n.doubs[t] <- sum(y[doubs, t])
  # }
  
  # Submodel 2 - whehter or not site is visited
  
  for(i in 1:nsites){
    for(t in 1:nyears){
      theta[i, t] <- ilogit(alpha.visit + beta.visit[1] * year[t] + beta.visit[2]*pow(year[t], 2) + 
                              beta.visit[3] * pow(year[t], 3) + kappa.lphi * lphi.site[i] + kappa.lgamma * lgamma.site[i] )
    }
  }
  
  alpha.visit <- logit(theta.int)
  theta.int ~ dbeta(1, 1)
  for(v in 1:3){
    beta.visit[v] ~ dnorm(0, sd = 2)
  }
  
  kappa.lphi ~ dnorm(0, sd = 2)
  kappa.lgamma ~ dnorm(0, sd = 2)
  
  # logisit regression for visits
  for(i in 1:nsites){
    for(t in 1:nyears){
      V[i, t] ~ dbern(theta[i, t])
    }
  }
})

inits <- function(){
  list(
    psi1 = rbeta(1, 1, 1),
    lphi.site = rnorm(constants$nsites, 0, 1), 
    lgamma.site = rnorm(constants$nsites, 0, 1), 
    initial.phi = rbeta(3, 1, 1), 
    initial.gamma = rbeta(3, 1, 1), 
    sd.lphi.site = rexp(1, 1), 
    sd.lgamma.site = rexp(1, 1),
    lphi.year = c(0, rnorm(constants$nyears-2, 0, 1)),
    lgamma.year = c(0, rnorm(constants$nyears-2, 0, 1)),
    sd.eps.lphi = rexp(1, 1), 
    sd.eps.lgamma = rexp(1, 1), 
    y = yst,
    theta.int = rbeta(1, 1, 1), 
    beta.visit = rnorm(3, 0, 2), 
    kappa.lphi = rnorm(1, 0, 2), 
    kappa.lgamma = rnorm(1, 0, 2)
  )
}

keepers <- c("psi1", "alpha.lphi", "initial.phi", 
             "alpha.lgamma", "initial.gamma", "lphi.year", 
             "lgamma.year", "mean.phi.year", "mean.gamma.year", 
             "sd.eps.lphi", "sd.eps.lgamma", "gamma.cliff", 
             "phi.cliff", "n.occ",
             # "n.ain",
             # "n.jura", "n.doubs",
             "kappa.lphi", "kappa.lgamma", "alpha.visit", "beta.visit")


# m <- nimble::nimbleModel(code = dynocc2,
#                          constants = constants,
#                          data = data,
#                          inits = inits())
# 
# m$initializeInfo()

# ART 33 minutes
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

# rhat of 744, wonder if it's converged
MCMCsummary(out, params = c("alpha.visit", "beta.visit"))


# yeah looks pretty bad. would have to run for much longer...kind of doubt this would converge ultimately?
MCMCsummary(out)

dynocc3 <- nimble::nimbleCode({
  
  psi1 ~ dbeta(1, 1)
  # model for phi and gamma - cliff hiehgt + site + smooth random year effects
  for(i in 1:nsites){
    for(t in 1:(nyears - 1)){
      logit(phi[i, t]) <- lphi[i, t]
      lphi[i, t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i, t]) <- lgamma[i, t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], sd.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], sd.lgamma.site)
  }
  
  # priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  
  sd.lphi.site ~ dexp(1)
  sd.lgamma.site ~ dexp(1)
  
  # priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0
  lgamma.year[1] <- 0
  
  for(t in 2:(nyears - 1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], sd.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], sd.eps.lgamma )
  }
  
  sd.eps.lphi ~ dexp(1)
  sd.eps.lgamma ~ dexp(1)
  
  
  # ecological and observation submodels confounded (no p)
  
  for(i in 1:nsites){
    y[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      y[i, t] ~ dbern(y[i, t-1]*phi[i, t-1] + (1 - y[i, t-1])*gamma[i, t-1])
    }
  }
  
  # derived parameters
  psi[1] <- psi1
  n.occ[1] <- sum(y[1:nsites, 1])
  for(t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites, t])
  }
  # year-specific average values of phi and gamma
  for( t in 1:(nyears - 1)){
    mean.phi.year[t] <- mean(phi[1:nsites, t])
    mean.gamma.year[t] <- mean(gamma[1:nsites, t])
  }
  
  # average gamma nd phi per cliff hiehgt category in central year
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }
  
  
  # pop size in each department
  # UGH, vector indexing is a pain in Nimble!!!
  # COuld probably accomplish this with a nimble funtion but I'm so over this
  # for(t in 1:nyears){
  #   n.ain[t] <- sum(y[ain, t])
  #   n.jura[t] <- sum(y[jura, t])
  #   n.doubs[t] <- sum(y[doubs, t])
  # }
  
  # Submodel 2 - whehter or not site is visited
  
  # no effect of past occupancy history in visitation at t = 1
  for(i in 1:nsites){
    theta[i, 1] <- ilogit(alpha.visit + beta.visit[1] * year[1] + beta.visit[2]*pow(year[1], 2) + beta.visit[3] * pow(year[1], 3) )
  }
  
  # including an effect of past occupancy from t = 2
  for(i in 1:nsites){
    for(t in 2:nyears){
      theta[i, t] <- ilogit(alpha.visit + beta.visit[1] * year[t] + beta.visit[2] * pow(year[t], 2) + beta.visit[3] * pow(year[t], 3) + 
                              kappa * y[i, t - 1] )
    }
  }
  
  alpha.visit <- logit(theta.int)
  theta.int ~ dbeta(1, 1)
  for(v in 1:3){
    beta.visit[v] ~ dnorm(0, sd = 2)
  }
  
  # coefficient for preferential sampling
  kappa ~ dnorm(0, sd = 2)
  
  # logisit regression for visits
  for(i in 1:nsites){
    for(t in 1:nyears){
      V[i, t] ~ dbern(theta[i, t])
    }
  }
})


inits <- function(){
  list(
    psi1 = rbeta(1, 1, 1),
    lphi.site = rnorm(constants$nsites, 0, 1), 
    lgamma.site = rnorm(constants$nsites, 0, 1), 
    initial.phi = rbeta(3, 1, 1), 
    initial.gamma = rbeta(3, 1, 1), 
    sd.lphi.site = rexp(1, 1), 
    sd.lgamma.site = rexp(1, 1),
    lphi.year = c(0, rnorm(constants$nyears-2, 0, 1)),
    lgamma.year = c(0, rnorm(constants$nyears-2, 0, 1)),
    sd.eps.lphi = rexp(1, 1), 
    sd.eps.lgamma = rexp(1, 1), 
    y = yst,
    theta.int = rbeta(1, 1, 1), 
    beta.visit = rnorm(3, 0, 2), 
    kappa = rnorm(1, 0, 2)
  )
}

keepers <- c("psi1", "alpha.lphi", "initial.phi", 
             "alpha.lgamma", "initial.gamma", "lphi.year", 
             "lgamma.year", "mean.phi.year", "mean.gamma.year", 
             "sd.eps.lphi", "sd.eps.lgamma", "gamma.cliff", 
             "phi.cliff", "n.occ",
             # "n.ain",
             # "n.jura", "n.doubs",
             "kappa", "alpha.visit", "beta.visit")

ni <- 50000
nb <- 40000
nt <- 5
nc <- 3

# m <- nimble::nimbleModel(code = dynocc3,
#                          constants = constants,
#                          data = data,
#                          inits = inits())
# 
# m$initializeInfo()

# ART 1.2 hours
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

# rhat of 5,000...
MCMCsummary(out)
