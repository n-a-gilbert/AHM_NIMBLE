# 4.10 - Analaysis of citizen science data with occupancy models

library(AHMbook)
library(nimble)
library(parallel)

# 4.10.1 - Temporal trends in detection heterogeneity

set.seed(1)

d <- AHMbook::simDynocc(
  nsites = 250, 
  nyears = 20, 
  nsurveys = 3, 
  mean.psi1 = 0.6, 
  range.p = c(0.5, 0.5), 
  range.phi = c(0.8, 0.8), 
  range.gamma = c(0.3, 0.3), 
  trend.sd.site = c(0, 2)
)

data <- list(
  y = d$y
)

constants <- list(
  nsites = dim(data$y)[1], 
  nsurveys = dim(data$y)[2], 
  nyears = dim(data$y)[3]
)

dynoccH <- nimble::nimbleCode({
  
  # Priors
  psi1 ~ dbeta(1, 1)
  for(t in 1:(nyears - 1)){
    phi[t] ~ dbeta(1, 1)
    gamma[t] ~ dbeta(1, 1)
  }
  for(t in 1:nyears){
    lp[t] <- logit(mean.p[t])
    mean.p[t] ~ dbeta(1, 1)
  }
  
  # random effects priors
  for(t in 1:nyears){
    for(i in 1:nsites){
      eps[i, t] ~ dnorm(0, sd.eps[t])
    }
    sd.eps[t] ~ dexp(1)
  }
  
  # ecological submodel
  for(i in 1:nsites){
    z[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      z[i, t] ~ dbern(z[i, t-1]*phi[t-1] + (1 - z[i, t-1])*gamma[t-1])
    }
  }
  
  # observation submodel
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      for(t in 1:nyears){
        logit(p[i, j, t]) <- lp[t] + eps[i, t] # time plus site-survey effects
        y[i, j, t] ~ dbern(z[i, t]*p[i, j, t])
      }
    }
  }
  
  # compute pop and sample occupancy
  psi[1] <- psi1
  psi.fs[1] <- sum(z[1:nsites, 1]) / 250 # sample occupancy
  for(t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1 - psi[t-1])*gamma[t-1]
    psi.fs[t] <- sum(z[1:nsites, t]) / 250
  }
})

inits <- function(){
  list(
    z = apply(data$y, c(1, 3), max), 
    psi1 = rbeta(1, 1, 1), 
    phi = rbeta(constants$nyears-1, 1, 1), 
    gamma = rbeta(constants$nyears-1, 1, 1), 
    mean.p = rbeta(constants$nyears, 1, 1), 
    eps = matrix(rnorm(constants$nsites*constants$nyears, 0, 1), nrow = constants$nsites, ncol = constants$nyears),
    sd.eps = rexp(constants$nyears, 1)
  )
}

keepers <- c("psi", "psi.fs", "phi", "gamma", "mean.p", "sd.eps")

ni <- 50000
nb <- 25000
nt <- 10
nc <- 3

# m <- nimbleModel(code = dynoccH,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()

# ART 
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("dynoccH",
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
  
  m <- nimbleModel(code = dynoccH,
                   name = "dynoccH",
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

MCMCvis::MCMCsummary(out)

# 4.10.2 analysis of CS data on Middle-spotted Woodpeckers

data("spottedWoodpecker")

dat <- spottedWoodpecker

# add to data scaled ate (such that 0 is 1 may and 1 unit is 1 month)
dat$date <- (dat$jdate - 121) / 30

# sample size in original data set
( nsites <- length(unique(dat$site)) )
( nyears <- length(unique(dat$year)) )
( ndays <- length(unique(dat$jdate)) )

dat.full <- dat
prop.data <- 0.3
ncase <- nrow(dat)

set.seed(1, sample.kind = "Rounding")
sel.cases <- sort(sample(1:ncase, ncase *prop.data ))

dat <- dat[sel.cases, ]

# have to renumber the sites since we lost some in subsampling
dat$site <- as.numeric(as.factor(dat[,"site"]))

#sample sites in new subsampled dataset
( nsites <- length(unique(dat$site)) )
( nyears <- length(unique(dat$ye)) )
( ndays <- length(unique(dat$jdate)) )

# compute detection frequency
table ( df <- tapply(dat$y, 
                     list(dat$site, dat$year), 
                     sum,
                     na.rm = TRUE)
)

# proportion of missing values per site x year combo
( prop.NA <- sum(is.na(df)) / prod(dim(df)) )

# proportion of missing value per site x year x dat combo
# lol most of them
( prop.NA <- 1 - (nrow(dat) / (nsites * nyears * ndays)) )

# compute observed occupancy
zobs <- tapply(dat$y, list(dat$site, dat$year), max, na.rm = TRUE)
zobs[zobs>1] <- 1
psiobs <- apply(zobs, 2, mean, na.rm = TRUE)
psiobs

data <- list(
  y = dat[, "y"],
  date = dat$date,
  nsurveys = dat[,"nsurvey"]
)

constants <- list(
  site = dat[,"site"],
  year = dat[, "year"]-1989, 
  nsites = nsites, 
  nyears = nyears, 
  nobs = nrow(dat)
)


# dynamic model with RE of year in phi, gamma, and p; plus, annual heterogeneity in site-level extra-dispersion in p
occmodel7 <- nimble::nimbleCode({
  
  # Priors
  psi1 ~ dbeta(1, 1)
  for(t in 1:(nyears - 1)){
    logit(phi[t]) <- lphi[t]
    lphi[t] ~ dnorm(mu.lphi, sd.lphi)
    logit(gamma[t]) <- lgamma[t]
    lgamma[t] ~ dnorm(mu.lgamma, sd.lgamma)
  }
  
  for(t in 1:nyears){
    alpha.lp[t] ~ dnorm(mu.alpha.lp,  sd.alpha.lp)
    beta.lp.1[t] ~ dnorm(mu.beta.lp1, sd.beta.lp1)
    beta.lp.2[t] ~ dnorm(mu.beta.lp2, sd.beta.lp2)
  }
  
  # HYPE
  
  mu.lphi <- logit(mean.phi)
  mu.lgamma <- logit(mean.gamma)
  mu.alpha.lp <- logit(mean.p)
  
  mean.phi    ~ dbeta(1, 1)
  mean.gamma  ~ dbeta(1, 1)
  mean.p      ~ dbeta(1, 1)
  mu.beta.lp1 ~ dnorm(0, sd = 0.5)
  mu.beta.lp2 ~ dnorm(0, sd = 0.5)
  sd.lphi     ~ dexp(1)
  sd.lgamma   ~ dexp(1)
  sd.alpha.lp ~ dexp(1)
  sd.beta.lp1 ~ dexp(1)
  sd.beta.lp2 ~ dexp(1)
  
  # anually varying site rnadom effects in detection
  for(t in 1:nyears){
    for(i in 1:nsites){
      eps.site[i, t] ~ dnorm(0, sd.p.site[year[i]])
    }
    sd.p.site[t] ~ dexp(1)
  }
  
  # Ecological submodel
  for(i in 1:nsites){
    z[i, 1] ~ dbern(psi1)
    for(t in 2:nyears){
      z[i, t] ~ dbern(z[i, t-1]*phi[t-1] + (1 - z[i, t-1])*gamma[t-1])
    }
  }
  
  # observation model
  for(i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]]*date[i] + beta.lp.2[year[i]]*pow(date[i], 2) + eps.site[site[i], year[i]]
    y[i] ~ dbin(z[site[i], year[i]] * p[i], nsurveys[i])
  }

  # derived parameters
  psi[1] <- psi1
  n.occ[1] <- sum(z[1:nsites, 1])
  for(t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1 - psi[t-1])*gamma[t-1]
    n.occ[t] <- sum(z[1:nsites, t])
  }
  
  
  })

keepers <- c(
  "psi", "phi", "gamma", 
  "n.occ", "mean.phi", "mu.lphi", 
  "sd.lphi", "mean.gamma", "mu.lgamma", 
  "sd.lgamma", "mean.p", "mu.alpha.lp", 
  "mu.beta.lp1", "mu.beta.lp2", "sd.alpha.lp", 
  "sd.beta.lp1", "sd.beta.lp2", "alpha.lp", 
  "beta.lp.1", "beta.lp.2", "sd.p.site"
)

zobs[is.na(zobs)] <- rbinom(sum(is.na(zobs)), 1, 0.5)
zst <- unname(zobs)

inits <- function(){
  list(
    z = zst,
    psi1= rbeta(1, 1, 1),
    lphi = rnorm(constants$nyears - 1, 0, 1),
    lgamma = rnorm(constants$nyears - 1, 0, 1), 
    alpha.lp = rnorm(constants$nyears, 0, 1), 
    beta.lp.1 = rnorm(constants$nyear, 0, 1), 
    beta.lp.2 = rnorm(constants$nyear, 0, 1), 
    mean.phi = rbeta(1, 1, 1), 
    mean.gamma = rbeta(1, 1, 1), 
    mean.p = rbeta(1, 1, 1), 
    mu.beta.lp1 = rnorm(1, 0, 0.5), 
    mu.beta.lp2 = rnorm(1, 0, 0.5), 
    sd.lphi = rexp(1, 1), 
    sd.lgamma = rexp(1, 1), 
    sd.alpha.lp = rexp(1, 1), 
    sd.beta.lp1 = rexp(1, 1), 
    sd.beta.lp2 = rexp(1, 1),
    eps.site = matrix(rnorm(constants$nsites * constants$nyears, 0, 1), nrow = constants$nsites, ncol = constants$nyears),
    sd.p.site = rexp(constants$nyears, 1) 
  )
}


ni <- 60000
nb <- 30000
nt <- 10
nc <- 3
 
# m <- nimbleModel(code = occmodel7,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
#  
# m$initializeInfo()

# ART 5 hours
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("occmodel7",
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
  
  m <- nimbleModel(code = occmodel7,
                   name = "occmodel7",
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

MCMCvis::MCMCsummary(out)
