# Chapter 3.4 - Spatial CJS

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)
library(spdep)

data("willowWarbler")
str(willowWarbler)

attach(willowWarbler)

ch <- as.matrix(birds[,1:11])
sitevec <- birds$cesID

table(apply(ch, 1 , sum))
apply(ch, 2, sum)
plot(table(table(sitevec)))
summary(as.numeric(table(sitevec)))

( nyear <- ncol(ch))
( nsite <- nrow(CES) )
( nblock <- nrow(blocks) )

( marr <- ch2marray(ch) )

(r <- apply(marr, 1, sum))

MARR <- array(NA, dim = c(10, nyear, nsite))
R <- array(NA, dim = c(10, nsite))
for(k in 1:nsite){
  sel.part <- ch[sitevec == k, ]
  ma <- ch2marray(sel.part)
  MARR[,,k] <- ma
  R[,k] <- apply(ma, 1, sum)
}
MARR

data <- list(
  MARR = MARR, 
  R = R
)

constants <- list(
  n.site = nsite, 
  n.occ = nyear
)

cjs6 <- nimbleCode({
  
  # priors and linear models
  for(s in 1:n.site){
    for(t in 1:(n.occ-1)){
      phi[t, s] <- ilogit(lphi[t, s]) # survival
      p[t, s] <- ilogit(lp[t, s]) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    
    # define random site effects
    alpha.lphi.site[s] ~ dnorm(mu.lphi, sd.lphi.site)
    alpha.lp.site[s]   ~ dnorm(mu.lp, sd.lp.site)
    mean.phi.site[s] <- ilogit(alpha.lphi.site[s])
    mean.p.site[s] <- ilogit(alpha.lp.site[s])
  }
  
  # define random year effects
  for(t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, sd.lphi.time)
    beta.lp.time[t]   ~ dnorm(0, sd.lp.time)
    mean.phi.time[t] <- ilogit(mu.lphi + beta.lphi.time[t])
    mean.p.time[t] <- ilogit(mu.lp + beta.lp.time[t])
  }
  
  
  # Hyperpriors
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dbeta(1, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lphi.site ~ dexp(1)
  sd.lp.site ~ dexp(1)
  sd.lphi.time ~ dexp(1)
  sd.lp.time ~ dexp(1)
  
  # multinomial likelihood for the m-array data
  for(s in 1:n.site){
    for(t in 1:(n.occ - 1)){
      MARR[t, 1:n.occ, s] ~ dmulti(pr[t, 1:n.occ, s], R[t, s])
    }
  }
  
  # define the cell probabilities of the m-array
  # main diagonal
  for(s in 1:n.site){
    for(t in 1:(n.occ - 1)){
      q[t, s] <- 1 - p[t, s]
      pr[t, t, s] <- phi[t, s]*p[t, s]
      # above the main diagona;
      for(j in (t+1):(n.occ-1)){
        pr[t, j, s] <- prod(phi[t:j, s])*prod(q[t:(j-1), s])*p[j, s]
      }
      # below the main diagonal
      for(j in 1:(t-1)){
        pr[t, j, s] <- 0
      }
    }
  }
  
  # last column: probabbility of non-recapture
  for(s in 1:n.site){
    for(t in 1:(n.occ-1)){
      pr[t, n.occ, s] <- 1 - sum(pr[t, t:(n.occ-1), s])
    }
  }
  
})

inits <- function(){
  list(
    alpha.lphi.site = rnorm(constants$n.site, 0, 2),
    alpha.lp.site = rnorm(constants$n.site, 0, 2),
    beta.lphi.time = rnorm(constants$n.occ - 1, 0, 2), 
    beta.lp.time = rnorm(constants$n.occ - 1, 0, 2),
    mean.phi = rbeta(1, 1, 1), 
    mean.p = rbeta(1, 1, 1), 
    sd.lphi.site = rexp(1, 1), 
    sd.lp.site = rexp(1, 1), 
    sd.lphi.time = rexp(1, 1), 
    sd.lp.time = rexp(1, 1)
  )
}

keepers <- c("mean.phi", "mean.p", "sd.lphi.site", "sd.lp.site", "sd.lphi.time", "sd.lp.time", "mean.phi.site", "mean.p.site", "mean.phi.time", "mean.p.time")

ni <- 30000
nb <- 20000
nt <- 10
nc <- 3

# m <- nimbleModel(code = cjs6,
# constants = constants,
# data = data,
# inits = inits())

# m$initializeInfo()
#35 mins
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs6",
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
  
  m <- nimbleModel(code = cjs6,
                   name = "cjs6", 
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


# 3.4.3 adding spatial covariates

scaled.gdd1 <- standardize(cells$gdd)
scaled.gdd2 <- standardize(cells$gdd^2)
scaled.lat <- standardize(cells$lat)

gdd1.site <- scaled.gdd1[CES$CellID]
gdd2.site <- scaled.gdd2[CES$CellID]
lat.site <- scaled.lat[CES$CellID]

data <- list(
  MARR = MARR, 
  R = R, 
  gdd1.site = gdd1.site, 
  gdd2.site = gdd2.site, 
  lat.site = lat.site
)

constants <- list(
  n.site = nsite, 
  n.occ = nyear
)

cjs7 <- nimbleCode({
  
  for(s in 1:n.site){
    for(t in 1:(n.occ-1)){
      phi[t, s] <- ilogit(lphi[t, s]) # survival
      p[t, s] <- ilogit(lp[t, s]) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lphi.time[t]
    }
    
    # linear model for site-level effects: add covariates
    alpha.lphi.site[s] ~ dnorm(mu.lphi.site[s], sd.lphi.site)
    mu.lphi.site[s] <- alpha.mu.lphi + beta1*gdd1.site[s] + beta2*gdd2.site[s] + beta3*lat.site[s]
    alpha.lp.site[s] ~ dnorm(mu.lp, sd.lp.site)
    
    # backtransform site means
    mean.phi.site[s] <- ilogit(alpha.lphi.site[s])
    mean.p.site[s] <- ilogit(alpha.lp.site[s])
  }
  
  # define year random effects
  for(t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, sd.lphi.time)
    beta.lp.time[t] ~ dnorm(0, sd.lp.time)
    
    # backtransform time means
    mean.phi.time[t] <- ilogit(alpha.mu.lphi + beta.lphi.time[t])
    mean.p.time[t] <- ilogit(mu.lp + beta.lp.time[t])
  }
  
  # hyperpriors
  alpha.mu.lphi <- logit(mean.phi)
  mean.phi ~ dbeta(1, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lphi.site ~ dexp(1)
  sd.lp.site ~ dexp(1)
  sd.lphi.time ~ dexp(1)
  sd.lp.time ~ dexp(1)
  
  # coefficients for gdd1, gdd2, and lat
  beta1 ~ dnorm(0, 1.5)
  beta2 ~ dnorm(0, 1.5)
  beta3 ~ dnorm(0, 1.5)
  
  # multinomial likelihood for the m-array
  for(s in 1:n.site){
    for(t in 1:(n.occ-1)){
      MARR[t, 1:n.occ, s] ~ dmulti(pr[t, 1:n.occ, s], R[t, s])
    }
  }
  
  # define the cell probabilities
  for(s in 1:n.site){
    for(t in 1:(n.occ-1)){
      q[t, s] <- 1 - p[t, s] # prob of non-recapture
      pr[t, t, s] <- phi[t, s]*p[t, s]
      # above diagonal
      for(j in (t+1):(n.occ-1)){
        pr[t, j, s] <- prod(phi[t:j, s])*prod(q[t:(j-1), s])*p[j, s]
      }
      # below the main diagonal
      for(j in 1:(t - 1)){
        pr[t, j, s] <- 0
      }
    }
  }
  
  # last column of m-array: probability of non-recapture
  for(s in 1:n.site){
    for(t in 1:(n.occ-1)){
      pr[t, n.occ, s] <- 1 - sum(pr[t, 1:(n.occ-1), s])
    }
  }
  
})

inits <- function(){
  list(
    alpha.lphi.site = rnorm(constants$n.site, 0, 1.5), 
    alpha.lp.site = rnorm(constants$n.site, 0, 1.5), 
    beta.lphi.time = rnorm(constants$n.occ-1, 0, 1.5),
    beta.lp.time = rnorm(constants$n.occ-1, 0, 1.5),
    mean.phi = rbeta(1, 1, 1),
    mean.p = rbeta(1, 1, 1), 
    sd.lphi.site = rexp(1, 1), 
    sd.lp.site = rexp(1, 1), 
    sd.lphi.time = rexp(1, 1), 
    sd.lp.time = rexp(1, 1), 
    beta1 = rnorm(1, 0, 1.5),
    beta2 = rnorm(1, 0, 1.5), 
    beta3 = rnorm(1, 0, 1.5)
  )
}

keepers <- c("mean.phi", "mean.p", 
             "alpha.mu.lphi", "mu.lp", 
             "sd.lphi.site", "sd.lphi.time", 
             "sd.lp.site", "sd.lp.time", 
             "mean.phi.site", "mean.p.site", 
             "mean.phi.time", "mean.p.time", "beta1", 
             "beta2", "beta3")

ni <- 50000
nb <- 40000
nt <- 5
nc <- 3

# m <- nimbleModel(code = cjs7,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# 35 mins
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs7",
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
  
  m <- nimbleModel(code = cjs7,
                   name = "cjs7", 
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

MCMCvis::MCMCsummary(out, params = c("mean.phi", "mean.p", "beta1", "beta2", "beta3"))

# 3.4.4 Spatial CJS with spatial autocorrelation only

blockcoordgrid <- cbind(as.matrix(blocks))
neigh <- spdep::dnearneigh(blockcoordgrid, d1 = 0, d2 = sqrt(2 * 25^2) + 1)
winnb <- spdep::nb2WB(neigh)

dim(MARR)

# this is the form following the book with the re-arranged m-array
# gonna stick with the JAGS form since Nimble can handle it apparently?

# dim( MARRWB <- aperm(MARR, c(3, 1, 2)) )
# 
# data <- list(
#   MARRWB = MARRWB,
#   R = R,
#   adj = winnb$adj,
#   weights = winnb$weights,
#   num = winnb$num)
# 
# constants <- list(
#   n.sites = nsite,
#   n.occ = nyear,
#   n.block = nblock,
#   BlockID = CES$BlockID,
#   neigh = length(data$adj)
# )
# 
# cjs8 <- nimbleCode({
# 
#   for(s in 1:n.sites){
#     for(t in 1:(n.occ - 1)){
#       phi[t, s] <- 1 / ( 1 + exp(-lphi[t, s])) # survival
#       p[t, s] <- 1 / (1 + exp(-lp[t, s])) # recapture
#       lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
#       lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
#     }
# 
#     # eta is spatial random effect at the block level
#     alpha.lphi.site[s] <- mu.lphi + eta[BlockID[s]]
#     alpha.lp.site[s] ~ dnorm(mu.lp, sd.lp.site)
# 
#     # backtransform site means
#     mean.p.site[s] <- 1 / (1 + exp(-alpha.lp.site[s]))
#   }
# 
#   for(t in 1:(n.occ-1)){
#     beta.lphi.time[t] ~ dnorm(0, sd.lphi.time)
#     beta.lp.time[t]   ~ dnorm(0, sd.lp.time)
# 
#     # backtransform time means
#     mean.phi.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lphi.time[t]))
#     mean.p.time[t] <- 1 / (1 + exp(-mu.lp + beta.lp.time[t]))
# 
#   }
# 
#   # hyperpriors
#   mu.lphi <- logit(mean.phi)
#   mean.phi ~ dbeta(1, 1)
#   mu.lp <- logit(mean.p)
#   mean.p ~ dbeta(1, 1)
#   sd.lp.site ~ dexp(1)
#   sd.lphi.time ~ dexp(1)
#   sd.lp.time ~ dexp(1)
# 
#   # CAR prior distribution
#   eta[1:n.block] ~ dcar_normal(adj[1:neigh], weights[1:neigh], num[1:n.block], tau, zero_mean = 1)
#   tau ~ dgamma(0.5, 0.0005)
#   veta <- 1 / tau
#   sdeta <- sqrt(veta)
# 
#   # multinomial likelihood for m-array
#   # open index comes last
#   # I bet I don't have to do this in Nimble, will try later
#   for(s in 1:n.sites){
#     for(t in 1:(n.occ-1)){
#       MARRWB[s, t, 1:n.occ] ~ dmulti(pr[t, s, 1:n.occ], R[t, s])
#     }
#   }
# 
#   # define the cell probabilitie sof the m-array
#   # main diagonal
#   for(s in 1:n.sites){
#     for(t in 1:(n.occ-1)){
#       q[t, s] <- 1 - p[t, s]
#       pr[t, s, t] <- phi[t, s]*p[t, s]
#       # above diagonal
#       for(j in (t+1):(n.occ-1)){
#         pr[t, s, j] <- prod(phi[t:j,s])*prod(q[t:(j-1), s])*p[j, s]
#       }
#       # below the diagonal
#       for(j in 1:(t-1)){
#         pr[t, s, j] <- 0
#       }
# 
#     }
#   }
# 
#   # last column: probability of non-recapture
#   for(s in 1:n.sites){
#     for(t in 1:(n.occ-1)){
#       pr[t, s, n.occ] <- 1 - sum(pr[t, s, 1:(n.occ-1)])
#     }
#   }
# 
# })
# 
# inits <- function(){
#   list(
#     alpha.lp.site = rnorm(constants$n.site, 0, 1),
#     beta.lphi.time = rnorm(constants$n.occ-1, 0, 1),
#     beta.lp.time = rnorm(constants$n.occ-1, 0, 1),
#     mean.phi = rbeta(1, 1, 1),
#     mean.p = rbeta(1, 1, 1),
#     sd.lp.site = rexp(1, 1),
#     sd.lphi.time = rexp(1, 1),
#     sd.lp.time = rexp(1, 1),
#     eta = rep(0, length(data$num)),
#     tau = rgamma(1, 0.5, 0.0005)
#   )
# }
# 
# keepers <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp", "sd.lp.site", "sd.lphi.time",
#              "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time", "veta", "sdeta", "eta")
# 
# ni <- 100000
# nb <- 50000
# nt <- 50
# nc <- 3
# 
# m <- nimbleModel(code = cjs8,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# start <- Sys.time()
# cl <- parallel::makeCluster(nc)
# 
# parallel::clusterExport(cl, c("cjs7",
#                               "inits",
#                               "data",
#                               "constants",
#                               "keepers",
#                               "nb",
#                               "ni",
#                               "nt"))
# 
# for(j in seq_along(cl)){
#   set.seed(j)
#   init <- inits()
#   clusterExport(cl[j], "init")
# }
# 
# out <- clusterEvalQ(cl, {
#   library(nimble)
#   
#   m <- nimbleModel(code = cjs7,
#                    name = "cjs7", 
#                    constants = constants,
#                    data = data,
#                    inits = init)
#   
#   Cmodel <- compileNimble(m)
#   modelConf <- configureMCMC(m)
#   modelConf$addMonitors(keepers)
#   modelMCMC <- buildMCMC(modelConf)
#   CmodelMCMC <- compileNimble(modelMCMC, project = m)
#   out1 <- runMCMC(CmodelMCMC,
#                   nburnin = nb,
#                   niter = ni,
#                   thin = nt)
#   
#   
# })
# 
# # turn off your cluster
# stopCluster(cl)
# end <- Sys.time()
# end - start

# trying with the regular MARR form - this seems to work
# must be an idiosyncracy of WinnBUGS

data <- list(
  MARR = MARR,
  R = R,
  adj = winnb$adj,
  weights = winnb$weights,
  num = winnb$num)

constants <- list(
  n.sites = nsite,
  n.occ = nyear,
  n.block = nblock,
  BlockID = CES$BlockID,
  neigh = length(data$adj)
)

cjs8 <- nimbleCode({

  for(s in 1:n.sites){
    for(t in 1:(n.occ - 1)){
      phi[t, s] <- 1 / ( 1 + exp(-lphi[t, s])) # survival
      p[t, s] <- 1 / (1 + exp(-lp[t, s])) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }

    # eta is spatial random effect at the block level
    alpha.lphi.site[s] <- mu.lphi + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, sd.lp.site)

    # backtransform site means
    mean.p.site[s] <- 1 / (1 + exp(-alpha.lp.site[s]))
  }

  for(t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, sd.lphi.time)
    beta.lp.time[t]   ~ dnorm(0, sd.lp.time)

    # backtransform time means
    mean.phi.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lphi.time[t]))
    mean.p.time[t] <- 1 / (1 + exp(-mu.lp + beta.lp.time[t]))

  }

  # hyperpriors
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dbeta(1, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lp.site ~ dexp(1)
  sd.lphi.time ~ dexp(1)
  sd.lp.time ~ dexp(1)

  # CAR prior distribution
  eta[1:n.block] ~ dcar_normal(adj[1:neigh], weights[1:neigh], num[1:n.block], tau, zero_mean = 1)
  tau ~ dgamma(0.5, 0.0005)
  veta <- 1 / tau
  sdeta <- sqrt(veta)

  # multinomial likelihood for the m-array
  for(s in 1:n.sites){
    for(t in 1:(n.occ-1)){
      MARR[t, 1:n.occ, s] ~ dmulti(pr[t, 1:n.occ, s], R[t, s])
    }
  }
  
  # define the cell probabilities
  for(s in 1:n.sites){
    for(t in 1:(n.occ-1)){
      q[t, s] <- 1 - p[t, s] # prob of non-recapture
      pr[t, t, s] <- phi[t, s]*p[t, s]
      # above diagonal
      for(j in (t+1):(n.occ-1)){
        pr[t, j, s] <- prod(phi[t:j, s])*prod(q[t:(j-1), s])*p[j, s]
      }
      # below the main diagonal
      for(j in 1:(t - 1)){
        pr[t, j, s] <- 0
      }
    }
  }
  
  # last column of m-array: probability of non-recapture
  for(s in 1:n.sites){
    for(t in 1:(n.occ-1)){
      pr[t, n.occ, s] <- 1 - sum(pr[t, 1:(n.occ-1), s])
    }
  }

})

inits <- function(){
  list(
    alpha.lp.site = rnorm(constants$n.site, 0, 1),
    beta.lphi.time = rnorm(constants$n.occ-1, 0, 1),
    beta.lp.time = rnorm(constants$n.occ-1, 0, 1),
    mean.phi = rbeta(1, 1, 1),
    mean.p = rbeta(1, 1, 1),
    sd.lp.site = rexp(1, 1),
    sd.lphi.time = rexp(1, 1),
    sd.lp.time = rexp(1, 1),
    eta = rep(0, length(data$num)),
    tau = rgamma(1, 0.5, 0.0005)
  )
}

keepers <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp", "sd.lp.site", "sd.lphi.time",
             "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time", "veta", "sdeta", "eta")

ni <- 100000
nb <- 50000
nt <- 50
nc <- 3

# m <- nimbleModel(code = cjs8,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
#  
# m$initializeInfo()

start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs8",
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

  m <- nimbleModel(code = cjs8,
                   name = "cjs8",
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

MCMCsummary(out, params = c("mean.phi", "mean.p"))

# 3.4.5 - Spatial hierarchical CJS with covariates!

data <- list(
  MARR = MARR, 
  R = R, 
  adj = winnb$adj, 
  weights = winnb$weights, 
  num = winnb$num,
  gdd1.site = gdd1.site, 
  gdd2.site = gdd2.site, 
  lat.site = lat.site
)

constants <- list(
  n.sites = nsite, 
  n.occ = nyear, 
  n.block = nblock,
  BlockID = CES$BlockID,
  neigh = length(data$adj)
)

cjs9 <- nimbleCode({
  
  for(s in 1:n.sites){
    for(t in 1:(n.occ - 1)){
      phi[t, s] <- 1 / ( 1 + exp(-lphi[t, s])) # survival
      p[t, s] <- 1 / (1 + exp(-lp[t, s])) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    
    # eta is spatial random effect at the block level
    alpha.lphi.site[s] <- mu.lphi + beta1*gdd1.site[s] + beta2*gdd2.site[s] + beta3*lat.site[s] + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, sd.lp.site)
    
    # backtransform site means
    mean.p.site[s] <- 1 / (1 + exp(-alpha.lp.site[s]))
  }
  
  for(t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, sd.lphi.time)
    beta.lp.time[t]   ~ dnorm(0, sd.lp.time)
    
    # backtransform time means
    mean.phi.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lphi.time[t]))
    mean.p.time[t] <- 1 / (1 + exp(-mu.lp + beta.lp.time[t]))
    
  }
  
  # hyperpriors
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dbeta(1, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lp.site ~ dexp(1)
  sd.lphi.time ~ dexp(1)
  sd.lp.time ~ dexp(1)
  
  beta1 ~ dnorm(0, 1.5)
  beta2 ~ dnorm(0, 1.5)
  beta3 ~ dnorm(0, 1.5)
  
  # CAR prior distribution
  eta[1:n.block] ~ dcar_normal(adj[1:neigh], weights[1:neigh], num[1:n.block], tau, zero_mean = 1)
  tau ~ dgamma(0.5, 0.0005)
  veta <- 1 / tau
  sdeta <- sqrt(veta)
  
  # multinomial likelihood for the m-array
  for(s in 1:n.sites){
    for(t in 1:(n.occ-1)){
      MARR[t, 1:n.occ, s] ~ dmulti(pr[t, 1:n.occ, s], R[t, s])
    }
  }
  
  # define the cell probabilities
  for(s in 1:n.sites){
    for(t in 1:(n.occ-1)){
      q[t, s] <- 1 - p[t, s] # prob of non-recapture
      pr[t, t, s] <- phi[t, s]*p[t, s]
      # above diagonal
      for(j in (t+1):(n.occ-1)){
        pr[t, j, s] <- prod(phi[t:j, s])*prod(q[t:(j-1), s])*p[j, s]
      }
      # below the main diagonal
      for(j in 1:(t - 1)){
        pr[t, j, s] <- 0
      }
    }
  }
  
  # last column of m-array: probability of non-recapture
  for(s in 1:n.sites){
    for(t in 1:(n.occ-1)){
      pr[t, n.occ, s] <- 1 - sum(pr[t, 1:(n.occ-1), s])
    }
  }
  
})

inits <- function(){
  list(
    alpha.lp.site = rnorm(constants$n.site, 0, 1),
    beta.lphi.time = rnorm(constants$n.occ-1, 0, 1),
    beta.lp.time = rnorm(constants$n.occ-1, 0, 1),
    mean.phi = rbeta(1, 1, 1),
    mean.p = rbeta(1, 1, 1),
    sd.lp.site = rexp(1, 1),
    sd.lphi.time = rexp(1, 1),
    sd.lp.time = rexp(1, 1),
    eta = rep(0, length(data$num)),
    tau = rgamma(1, 0.5, 0.0005),
    beta1 = rnorm(1, 0, 1.5), 
    beta2 = rnorm(1, 0, 1.5), 
    beta3 = rnorm(1, 0, 1.5)
  )
}

keepers <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp", "sd.lp.site", "sd.lphi.time",
             "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time", "veta", "sdeta", "eta", "beta1", "beta2", "beta3")

ni <- 100000
nb <- 50000
nt <- 50
nc <- 3

# m <- nimbleModel(code = cjs9,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

#2.7 hours
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("cjs9",
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
  
  m <- nimbleModel(code = cjs9,
                   name = "cjs9",
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

MCMCsummary(out, params = c("mean.phi", "mean.p", "beta1", "beta2", "beta3"))
