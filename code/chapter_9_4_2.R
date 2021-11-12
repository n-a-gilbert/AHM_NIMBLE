# Chapter 9.4.2 - GAM for spatial autocorrelation
# 12 November 2021

library(AHMbook)
library(fields)
library(raster)
library(nimble)
library(parallel)
library(MCMCvis)

RNGversion("3.5.3")

str( dat <- AHMbook::simNmixSpatial(
  nsurveys = 3, 
  mean.lambda = exp(2), 
  beta = c(2, -2), 
  mean.p = 0.5, 
  alpha = c(-1, -1), 
  sample.size = 500, 
  variance.RF = 1, 
  theta.RF = 10, 
  seeds = c(10, 100), 
  truncN = 6, 
  show.plots = TRUE
) )

# scale both sets of coordinates
head( coordgrid <- scale(cbind(dat$xcoord, dat$ycoord)) ) # all 2500 cells
head( sitelocs <- coordgrid[dat$surveyed.sites, ]) # 500 surveyed cells

# takes a few minutes
system.time(
  knots <- fields::cover.design(R = coordgrid, 
                                nd = 125, 
                                nruns = 10, 
                                num.nn = 200, 
                                max.loop = 20)
)

# define the Z matrix for knot coefficients
knotlocs <- knots$design
omega <- ( e2dist(knotlocs, knotlocs)/10)^3
svd.omega <- svd(omega)
sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
Zk <- (e2dist(coordgrid, knotlocs)/10)^3
Zmat <- t(solve(sqrt.omega, t(Zk)))

head(Zmat)

n.knots <- 125
y <- dat$yobs

data <- list(
  y = y, 
  elev = dat$elevationS, 
  forest = dat$forestS, 
  wind = dat$wind, 
  Zmat = Zmat
  )

constants <- list(
  nsites = dim(y)[1], 
  nreps = dim(y)[2], 
  n.knots = n.knots
)

# code slightly modified (mostly slightly different priors)
# seem to be some problems as is -
# poor mixing of knot coefficients, and some trouble with abundance intercept/fixed effects
# possible solutions - 1) dropping beta0 intercept and let spatial field function as intercept
# and/or reduce the number of knots - not sure if that would help model along ?
splines_nmix <- nimble::nimbleCode({
  
  # intercepts
  beta0  ~ dnorm(0, 5)
  alpha0 ~ dnorm(0, 5)
  
  # fixed effects
  for(v in 1:2){
    alpha[v] ~ dnorm(0, sd = 1) 
    beta[v]  ~ dnorm(0, sd = 1)
  }
  
  # priors for random effects representing the splines
  for(k in 1:n.knots){
    b[k] ~ dnorm(0, tau.b)
  }
  
  # prior for random effects dispersion
  tau.b ~ dexp(1)
  
  # Model for abundance
  for(i in 1:nsites){
    N[i] ~ dpois(lam[i])
    log(lam[i]) <- beta0 + beta[1]*elev[i] + beta[2]*pow(elev[i], 2) + smooth[i]
    # authors note that subtracting mean of spatial field from spatial field
    # at every iteration to avoid confounding with the intercept
    smooth[i] <- smooth2[i] - mean(smooth2[1:nsites])
    smooth2[i] <- inprod(Zmat[i, 1:n.knots], b[1:n.knots])
  }
  
  # measurement error model
  for(i in 1:nsites){
    for(j in 1:nreps){
      y[i, j] ~ dbin(p[i, j], N[i])
      logit(p[i, j]) <- alpha0 + alpha[1]*forest[i] + alpha[2]*wind[i, j]
    }
  }

  # derived parameter: total population of grid
  Ntotal <- sum(N[1:nsites])
  
})

# initial values - start at max observed count (or 2 if unsurveyed or none detected)
Nst <- apply(dat$yobs, 1, max)
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2

inits <- function(){
  list(
    N = Nst, 
    beta0 = rnorm(1, 0, 5), 
    alpha0 = rnorm(1, 0, 5), 
    beta = rnorm(2, 0, 1), 
    alpha = rnorm(2, 0, 1), 
    b = rnorm(constants$n.knots, 0, 1), 
    tau.b = rexp(1, 1)
  )
}

keepers <- c("beta0", "beta", "alpha0", "alpha", "Ntotal", "tau.b", "b")

# MCMC settings
# With these settings, it ran for ~ 20 hours (parallelized)
# this is MUCH longer than the runtime noted by authors - not sure what's going on
# (JAGS / Nimble differences, or maybe I'm messing something up?)
nc <- 3
ni <- 50000
nb <- 40000
nt <- 2  

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

# mixing seems okay for fixed effects and intercepts
# abundance interept (beta0) having some trouble
# as well as coefficient for forest (alpha[1])
MCMCvis::MCMCsummary(out,
                     params = c("beta0", "alpha0", "beta", "alpha"))

# mixing problems with knot coefficients
# only ~10% have rhat < 1.1
MCMCvis::MCMCsummary(out, 
                     params = c("b"))

# traceplots
MCMCvis::MCMCtrace(out, 
                   params = c("beta0", "alpha0", "beta", "alpha"),
                   pdf = FALSE,
                   type = "trace")

# visualization of estimates
# looking okay-ish
MCMCvis::MCMCplot(out, 
                  params = c("beta0", "alpha0", "beta", "alpha"), 
                  ci = c(50, 95))

# yikes, crazy estimates for the knot coefficients
MCMCvis::MCMCplot(out, 
                  params = c("b"), 
                  ci = c(50, 95))