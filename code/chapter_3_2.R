library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

set.seed(1)
str(data <- simCJS(n.occ = 6, n.marked = 20, phi = 0.7, p = 0.4, show.plot = TRUE))

simCJS()

simCJS(n.occ = 20, phi = rep(0.7, 19), p = rep(0.4, 19))
simCJS(n.occ = 6, n.marked = c(10, 20, 30, 40, 50))

simCJS(n.occ = 20, phi= 0.9, p = 0.4)

simCJS(n.occ = 20, phi = runif(19, 0.4, 0.8),
       p = runif(19, 0.2, 0.5))

phi <- plogis(rnorm(19, qlogis(0.8), 1))
simCJS(n.occ = 20, phi = phi)


set.seed(1)
data <- simCJS(show.plot = FALSE)

bdata <- list(y = data$ch)
constants <- list(f = data$f, 
                  n.ind = data$n.ind, 
                  n.occ = data$n.occ)

cjs <- nimbleCode({
  
  # priors
  phi ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
  
  # likelihood
  for(i in 1:n.ind){
    z[i, f[i]] <- 1
    for(t in (f[i] + 1):n.occ){
      z[i, t] ~ dbern(z[i, t-1]*phi)
      y[i, t] ~ dbern(z[i, t]*p)
    }
  }
  
})

inits <- function(){
  list(
    z = zinit(data$ch), 
    phi = rbeta(1, 1, 1), 
    p = rbeta(1, 1, 1)
  )
}

keepers <- c("phi", "p", "z")

ni <- 50000
nb <- 30000
nt <- 10
nc <- 3

# m <- nimbleModel(code = cjs,
#                  constants = constants,
#                  data = bdata,
#                  inits = inits())
# 
# m$initializeInfo()

# ran in <1 min, but should bump up iterations to get better mixing
# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs",
                              "inits",
                              "bdata",
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
  
  m <- nimbleModel(code = cjs,
                   name = "cjs", 
                   constants = constants,
                   data = bdata,
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

MCMCsummary(out, params = c("phi", "p"))

# 3.2.4 time-dependent CJS model

set.seed(24)
dat <- simCJS(n.occ = 6, n.marked = 100, phi = runif(5, 0.2, 0.8), p = runif(5, 0.2, 0.9))

data <- list(
  y = dat$ch
)
constants <- list(
  f = dat$f, 
  n.ind = dat$n.ind, 
  n.occ = dat$n.occ
)

cjs2 <- nimbleCode({
  
  # priors
  for(t in 1:(n.occ-1)){
    phi[t] ~ dbeta(1, 1) #apparent survival
    p[t] ~ dbeta(1, 1)   # recap
  }
  
  for(i in 1:n.ind){
    z[i, f[i]] <- 1
    for(t in (f[i] + 1):n.occ){
      z[i, t] ~ dbern(z[i, t-1]*phi[t-1])
      y[i, t] ~ dbern(z[i, t] * p[t-1])
    }
  }
  
})

inits <- function(){
  list(
    z = zinit(dat$ch),
    phi = rbeta(constants$n.occ-1, 1, 1), 
    p = rbeta(constants$n.occ-1, 1, 1)
  )
}

keepers <- c("phi", "p", "z")
ni <- 50000
nb <- 30000
nt <- 10
nc <- 3

# m <- nimbleModel(code = cjs2,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# <2 mins
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs2",
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
  
  m <- nimbleModel(code = cjs2,
                   name = "cjs2", 
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

MCMCsummary(out, params = c("phi", "p"))
