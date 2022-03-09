# 10.8 Integrated model with SCR data

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)

# make a trapping grid
traplocs <- expand.grid(1:10, 1:10)
ntraps <- nrow(traplocs)

# define state-space by buffering traps
delta <- 2 # buffer width
Xl <- min(traplocs[,1] - delta) # lower x
Xu <- max(traplocs[,1] + delta) # upper x
Yl <- min(traplocs[,2] - delta) # lower y
Yu <- max(traplocs[,2] + delta) # upper y

# Distribute pop size of 50 in state space
N <- 50
K <- 4 # number of weeks of hair collection

# simulate activity centers
set.seed(121, kind = "Mersenne-Twister")
sx <- runif(N, Xl, Xu)
sy <- runif(N, Yl, Yu)
smat <- cbind(sx, sy)

# only some traps collect data
inner36 <- traplocs[,1] >= 3 & traplocs[,1] <=8 & traplocs[,2] >= 3 & traplocs[,2] <= 8
# others only collect count data
outer64 <- 1:ntraps
outer64 <- outer64[!inner36]

# parameters for the SCR model simulation
lam0 <- 0.6 # baseline encounter rate
sigma <- 0.6 # half-normal detection scale

# generate the encounters of every individual in every trap
D <- e2dist(smat, traplocs)
muy <- lam0 * exp(-(D*D) / (2*sigma^2))
Y <- matrix(NA, nrow = N, ncol = ntraps)
for(i in 1:N){
  Y[i,] <- rpois(ntraps, K*muy[i,])
}

# grab scr data from the inner 36 traps
Yscr <- Y[,inner36]
# captured individuals appear in the dataset
totalcaps <- apply(Yscr, 1, sum)
Yscr <- Yscr[totalcaps > 0, ]

#only count data appear in other set
Yocc <- Y[,outer64]
n <- apply(Yocc, 2, sum)

scrtraps <- traplocs[inner36,]
occtraps <- traplocs[outer64,]


# 10.8.3 - just the scr data
# set up data augmentation
M <- 150
Yaug <- matrix(0, M, dim(scrtraps)[1])
Yaug[1:dim(Yscr)[1],] <- Yscr

data <- list(
  y = Yaug, 
  Xl = Xl, 
  Yl = Yl, 
  Xu = Xu, 
  Yu = Yu, 
  X = as.matrix(scrtraps))

constants <- list(
  M = M, 
  nsurveys = K,
  ntraps = nrow(scrtraps)
)

scr <- nimble::nimbleCode({
  
  # Priors
  lam0 ~ dunif(0, 5)
  sigma ~ dunif(0, 10)
  psi ~ dunif(0, 1)
  
  # "Likelihood"
  for(i in 1:M){
    # process model
    z[i] ~ dbern(psi) # existence of individual
    s[i, 1] ~ dunif(Xl, Xu) # and its location, x
    s[i, 2] ~ dunif(Yl, Yu) # and y coordinate of activity center
    # observation model
    for(j in 1:ntraps){
      d[i, j] <- pow(pow(s[i, 1] - X[j, 1], 2) + pow(s[i, 2] - X[j, 2], 2), 0.5)
      lambda[i, j] <- z[i]*lam0*exp(-(d[i, j]*d[i, j])/(2*sigma*sigma))
      y[i, j] ~ dpois(nsurveys*lambda[i, j])
    }
  }
  N <- sum(z[1:M])
  
})

SinX <- runif(M, Xl, Xu)
SinY <- runif(M, Yl, Yu)
for(i in 1:dim(Yscr)[1]){
  SinX[i] <- sum(Yscr[i,]*scrtraps[,1]) / (sum(Yscr[i,]))
  SinY[i] <- sum(Yscr[i,]*scrtraps[,2]) / (sum(Yscr[i,]))
}

inits <- function(){
  list(
    z = c(rep(1, dim(Y)[1]), rbinom(M - dim(Y)[1], 1, 0.5)), 
    psi = runif(1), 
    s = cbind(SinX, SinY), 
    lam0 = runif(1, 0.5, 1.5), 
    sigma = runif(1, 0.5, 3)
  )
}

keepers <- c("psi", "lam0", "sigma", "N")

ni <- 30000
nb <- 10000
nt <- 20
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = scr,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()


# took ~ 5 mins with these settings
# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("scr",
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
  
  m <- nimbleModel(code = scr,
                   name = "scr", 
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

data <- list(
  n = n,
  Xl = Xl, 
  Yl = Yl, 
  Xu = Xu, 
  Yu = Yu, 
  X = as.matrix(occtraps)
)

constants <- list(
  M = M, 
  nsurveys = K, 
  ntraps = nrow(occtraps)
)

cr <- nimbleCode({
  
  # priors
  lam0 ~ dunif(0, 5)
  sigma ~ dunif(0, 10)
  psi ~ dbeta(1, 1)
  
  # "likelihood"
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i, 1] ~ dunif(Xl, Xu)
    s[i, 2] ~ dunif(Yl, Yu)
  }
  
  # observation model
  for(j in 1:ntraps){
    for(i in 1:M){
      d2[i, j] <- pow(s[i, 1] - X[j, 1], 2) + pow(s[i, 2] - X[j, 2], 2)
      lam[i, j] <- lam0*exp(-(d2[i, j]) / (2*sigma*sigma))*z[i]
    }
    lambda[j] <- sum(lam[1:M, j])
    n[j] ~ dpois(nsurveys*lambda[j])
  }
  N <- sum(z[1:M])
  
  
})

Sx <- runif(M, Xl, Xu)
Sy <- runif(M, Yl, Yu)
inits <- function(){
  list(
    z = as.vector(rep(1, M)), 
    psi = runif(1), 
    s = cbind(Sx, Sy), 
    lam0 = runif(1, 0.5, 1.5), 
    sigma = runif(1, 0.5, 3)
  )
}

keepers <- c("psi", "lam0", "sigma", "N")

ni <- 20000
nt <- 2
nb <- 15000
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
# m <- nimbleModel(code = cr,
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

parallel::clusterExport(cl, c("cr",
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
  
  m <- nimbleModel(code = cr,
                   name = "cr", 
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

data <- list(
  yscr = Yaug, 
  n = n, 
  Xl = Xl, 
  Yl = Yl, 
  Xu = Xu, 
  Yu = Yu, 
  scrtraps = as.matrix(scrtraps),
  occtraps = as.matrix(occtraps))

constants <- list(
  M = M, 
  nsurveys = K, 
  n.scrtraps = nrow(scrtraps), 
  n.occtraps = length(n)
)

im_scr <- nimbleCode({
  
  # Priors
  psi ~ dunif(0, 1)
  lam0 ~ dunif(0, 5)
  sigma ~ dunif(0, 10)
  
  # "likelihood"
  for(i in 1:M){
    # shared process model for both datasets
    z[i] ~ dbern(psi)
    s[i, 1] ~ dunif(Xl, Xu)
    s[i, 2] ~ dunif(Yl, Yu)
    
    # observation model for the SCR data
    for(j in 1:n.scrtraps){
      d[i, j] <- pow(pow(s[i, 1] - scrtraps[j, 1], 2) + pow(s[i, 2] - scrtraps[j, 2], 2), 0.5)
      lambda1[i, j] <- z[i]*lam0*exp(-(d[i, j]*d[i, j]) / (2*sigma*sigma))
      yscr[i, j] ~ dpois(nsurveys*lambda1[i, j]) 
    }
  }
  
  # observation model for unmarked counts
  for(j in 1:n.occtraps){
    for(i in 1:M){
      d2[i, j] <- pow(s[i, 1] - occtraps[j, 1], 2) + pow(s[i, 2] - occtraps[j, 2], 2)
      lam[i, j] <- z[i]*lam0*exp(-(d2[i, j]) / (2*sigma*sigma))
    }
    lambda2[j] <- sum(lam[1:M, j]) # expected captures / occasion at trap t
    n[j] ~ dpois(nsurveys*lambda2[j])
  }
  
  N <- sum(z[1:M])
  
})

# inits for activity centers
SinX <- SinY <- rep(NA, M)
for(i in 1:dim(Yscr)[1]){
  SinX[i] <- sum(Yscr[i,]*scrtraps[,1]) / (sum(Yscr[i, ]))
  SinY[i] <- sum(Yscr[i,]*scrtraps[,2]) / (sum(Yscr[i, ]))
}

inits <- function(){
  list(
    z = c(rep(1, dim(Y)[1]), rbinom(M - dim(Y)[1], 1, 0.5)),
    psi = runif(1), 
    s = cbind(SinX, SinY), 
    sigma = runif(1, 0.5, 3),
    lam0 = runif(1, 0, 5)
  )
}

keepers <- c("psi", "lam0", "sigma", "N")

ni <- 10000
nb <- 6000
nt <- 4
nc <- 3

# Check that everything is in order prior to setting up things to go parallel
m <- nimbleModel(code = im_scr,
                 constants = constants,
                 data = data,
                 inits = inits())

m$initializeInfo()


# took ~ 5 mins with these settings
# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("im_scr",
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

out3 <- clusterEvalQ(cl, {
  library(nimble)
  
  m <- nimbleModel(code = im_scr,
                   name = "im_scr", 
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

MCMCsummary(out3)

# 10.8.6 - SCR + occupancy data

# generate the encounters of every individual in every trap
set.seed(121, kind = "Mersenne-Twister")

D <- e2dist(smat, traplocs)
muy <- lam0 *exp(-(D*D) / (2*sigma^2))
Y <- array(NA, dim = c(N, ntraps, K))
for(i in 1:N){
  for(k in 1:K){
    Y[i, , k] <- rpois(ntraps, muy[i,])
  }
}

Yscr <- Y[, inner36, ]
# captured individuals appear in the dataset
totalcaps <- apply(Yscr, 1, sum)
nind <- sum(totalcaps > 0)
Yscr <- Yscr[totalcaps > 0, , ]

# only count data appear in the other dataset
Yocc <- Y[, outer64, ]
# total number of detections observed in outer traps
Yocc <- apply(Yocc, c(2, 3), sum)

# set up data augmentation for the encounter histories
M <- 150
Yaug <- array(0, dim = c(M, dim(scrtraps)[1], K))
Yaug[1:dim(Yscr)[1],,] <- Yscr

# convert to binary presence/absencec data
Yaug[Yaug > 1] <- 1
Yocc[Yocc > 1] <- 1

data <- list(
  yscr = Yaug, 
  yocc = Yocc, 
  Xl = Xl, 
  Yl = Yl, 
  Xu = Xu, 
  Yu = Yu,
  scrtraps = as.matrix(scrtraps), 
  occtraps = as.matrix(occtraps)
)

constants <- list(
  M = M, 
  nsurveys = K, 
  n.scrtraps = nrow(scrtraps), 
  n.occtraps = length(n)
)

scr_occ <- nimbleCode({
  
  psi ~ dbeta(1, 1)
  lam0 ~ dunif(0, 5)
  lam0occ ~ dunif(0, 5)
  sigma ~ dunif(0, 10)
  N <- sum(z[1:M])
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i, 1] ~ dunif(Xl, Xu)
    s[i, 2] ~ dunif(Yl, Yu)
    # compute detection probabilyt for scr
    for(j in 1:n.scrtraps){
      d2scr[i, j] <- (s[i, 1] - scrtraps[j, 1])^2 + (s[i, 2] - scrtraps[j, 2])^2
      lam[i, j] <- lam0*exp(-d2scr[i, j] / (2*sigma^2))
      pscr[i, j] <- 1 - exp(-lam[i, j]) 
    }
    
    # compute detection probability for occupancy
    for(j in 1:n.occtraps){
      d2occ[i, j] <- (s[i, 1] - occtraps[j, 1])^2 + (s[i, 2] - occtraps[j, 2])^2
      lam2[i, j] <- lam0occ*exp(-d2occ[i, j] / (2*sigma^2))
      pocc[i, j] <- 1 - exp(-lam2[i, j])
    }
    for(j in 1:n.scrtraps){
      for(k in 1:nsurveys){
        yscr[i, j, k] ~ dbern(pscr[i, j]*z[i])
      }
    }
    for(j in 1:n.occtraps){
      # for PA data compute probability of not captured
      pn[i, j] <- (1 - (pocc[i, j]*z[i]))
    }
  }
  
  # model for presence / absence data
  for(j in 1:n.occtraps){
    for(k in 1:nsurveys){
      yocc[j, k] ~ dbern(1 - prod(pn[1:M, j]))
    }
  }
})

SinX <- runif(M, Xl, Xu)
SinY <- runif(M, Yl, Yu)
ncaps <- apply(Yscr, c(1, 2), sum)
ncaps[ncaps > 1] <- 1
for(i in 1:nind){
  SinX[i] <- sum(ncaps[i,]*scrtraps[,1]) / (sum(ncaps[i,]))
  SinY[i] <- sum(ncaps[i,]*scrtraps[,2]) / (sum(ncaps[i,]))
}

inits <- function(){
  list(
    z = c(rep(1, dim(Y)[1]), rbinom(M - dim(Y)[1], 1, 0.5)), 
    psi = runif(1), 
    s = cbind(SinX, SinY), 
    sigma = runif(1, 0.5, 3),
    lam0occ = runif(1, 0.2, 0.5),
    lam0 = runif(1, 0, 5)
    
  )
}

keepers <- c("psi", "lam0", "lam0occ", "sigma", "N")

ni <- 12000
nb <- 6000
nt <- 6
nc <- 3


# Check that everything is in order prior to setting up things to go parallel
m <- nimbleModel(code = scr_occ,
                 constants = constants,
                 data = data,
                 inits = inits())

m$initializeInfo()


# took ~ 5 mins with these settings
# code to run chains in parallel
# also note, this can be modified to extend chains if they are not yet converged
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("scr_occ",
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

out4 <- clusterEvalQ(cl, {
  library(nimble)
  
  m <- nimbleModel(code = scr_occ,
                   name = "scr_occ", 
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
