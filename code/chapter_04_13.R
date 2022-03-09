library(AHMbook)
library(nimble)
library(parallel)

data("SwissMarbledWhite")

dat <- SwissMarbledWhite

y <- as.matrix(dat[,14:24])

DATE <- as.matrix(dat[,3:13]) # dates on which surveys were conducted

for(t in 1:11){
  DATE[is.na(DATE[, t]), t] <- mean(DATE[,t], na.rm = T)
}

year <- dat$year

nsites <- length(unique(dat$site))

nyears <- length(unique(dat$year))

nsurveys <- ncol(y)

nobs <- nrow(y)

data <- list(
  y = y, 
  DATE = DATE, 
  yr = year - 2004#,
  # site = dat$site
)

constants <- list(
  nobs = nobs, 
  nsites = nsites, 
  nyears = nyears, 
  year = year - 1997, 
  nsurveys = nsurveys
)

pheno_occ <- nimble::nimbleCode({
  
  # linear model for annual site occupancy with its priors
  for(i in 1:nobs){
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- min(10, max(-10, lpsi[i]))
    lpsi[i] <- beta.lpsi[1] + beta.lpsi[2]*yr[i] + eps.lpsi[year[i]]
  }
  
  for(t in 1:nyears){
    eps.lpsi[t] ~ dnorm(0, sd.lpsi)
  }
  
  # priors for occ
  beta.lpsi[1] <- logit(mean.psi)
  mean.psi ~ dbeta(1, 1)
  beta.lpsi[2] ~ dnorm(0, sd = 1)
  sd.lpsi ~ dexp(1)
  
  # logit-linear model of detection on year
  for(t in 1:nyears){
    for(j in 1:nsurveys){
      logit(p[t, j]) <- min(10, max(-10, lp[t, j]))
      lp[t, j] <- beta.lp[1, j] + beta.lp[2, j]*(t - 7) # year centered
    }
  }
  
  # priors for detection
  for(j in 1:nsurveys){
    beta.lp[1, j] ~ dnorm(mu.lp1, sd.lp1)
    beta.lp[2, j] ~ dnorm(mu.lp2, sd.lp2)
  }
  
  mu.lp1 ~ dnorm(0, sd = 0.5)
  mu.lp2 ~ dnorm(0, sd = 0.5)
  sd.lp1 ~ dexp(1)
  sd.lp2 ~ dexp(1)
  
  # linear regression of arrival date on year
  for(i in 1:nobs){
    arr[i] ~ dnorm(mu.arr1[i], sd.arr)
    mu.arr1[i] <- min(200, max(30, mu.arr[i]))
    mu.arr[i] <- beta.arr[1] + beta.arr[2]*yr[i]
  }
  
  # priors for arrival submodel
  beta.arr[1] ~ dnorm(90, sd = 10)
  beta.arr[2] ~ dnorm(0, sd = 5)
  
  sd.arr ~ dexp(1)
  
  # linear regression for departure date on year
  for(i in 1:nobs){
    dep[i] ~ dnorm(mu.dep1[i], sd.dep)
    mu.dep1[i] <- min(500, max(0, mu.dep[i]))
    mu.dep[i] <- beta.dep[1] + beta.dep[2] * yr[i]
  }
  
  # priors for departure model
  beta.dep[1] ~ dnorm(120, sd = 10)
  beta.dep[2] ~ dnorm(0, sd = 5)
  
  sd.dep ~ dexp(1)
  
  # model for the observed data
  for(i in 1:nobs){
    for(j in 1:nsurveys){
      y[i, j] ~ dbern(mu1[i, j])
      mu1[i, j] <- min(0.99, max(0.01, mu[i, j]))
      mu[i, j] <- z[i] * step(DATE[i, j] - arr[i]) * step(dep[i] - DATE[i, j]) * p[year[i], j]
    }
  }
  
  # derived quantities
  # average occupancy per year
  for(t in 1:nyears){
    for(s in 1:nsites){
      logit(tmp[s, t]) <- beta.lpsi[1] + beta.lpsi[2] * (t - 7) + eps.lpsi[t]
    }
    
    psi.pred[t] <- mean(tmp[1:nsites, t])
    
  }
  
  # average detection per year and visit
  for(t in 1:nyears){
    for(j in 1:nsurveys){
      logit(p.pred[t, j]) <- beta.lp[1, j] + beta.lp[2, j] * (t - 7)
    }
  }
  
  # average arrival and departure time per year and length of flight period
  for(t in 1:nyears){
    arr.pred[t] <- beta.arr[1] + beta.arr[2] * (t - 7)
    dep.pred[t] <- beta.dep[1] + beta.dep[2] * (t - 7)
    fp.pred[t] <- dep.pred[t] - arr.pred[t]
  }
  
})

zst <- apply(y, 1, max, na.rm = T)
inits <- function(){
  list(
    z            = zst, 
    eps.lpsi     = rnorm(constants$nyears, 0, 1), 
    mean.psi     = rbeta(1, 1, 1), 
    beta.lpsi    = rnorm(2, mean = 0, sd = 2),
    sd.lpsi      = rexp(1, 1), 
    beta.lp      = matrix(rnorm(2*constants$nsurveys),
                          ncol = constants$nsurveys,
                          nrow = 2),
    mu.lp1       = rnorm(1, 0, 1), 
    mu.lp2       = rnorm(1, 0, 1), 
    sd.lp1       = rexp(1, 1), 
    sd.lp2       = rexp(1, 1), 
    arr          = rnorm(constants$nobs, 90, 10), 
    beta.arr     = c(rnorm(1, 90, 10), rnorm(1, 0, 5)),
    sd.arr       = rexp(1, 1),
    dep          = rnorm(constants$nobs, 120, 10), 
    beta.dep     = c(rnorm(1, 120, 10), rnorm(1, 0, 5)),
    sd.dep       = rexp(1, 1)
  )
}

keepers <- c("mean.psi", "beta.lpsi", "sd.lpsi", "beta.lp", "mu.lp1", 
             "mu.lp2", "sd.lp1", "beta.arr", "sd.arr", "beta.dep", "sd.dep", 
             "psi.pred", "p.pred", "arr.pred", "dep.pred", "fp.pred")

nb <- 50000
ni <- nb + 50000
nt <- 50
nc <- 3


m <- nimble::nimbleModel(code = pheno_occ,
                         constants = constants,
                         data = data,
                         inits = inits())

m$initializeInfo()

# ART 
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("pheno_occ",
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
  
  m <- nimbleModel(code = pheno_occ,
                   name = "pheno_occ",
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
