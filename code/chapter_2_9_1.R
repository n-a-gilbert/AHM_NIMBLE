# 2.9.1 Disease-structured Dail-Madsen model

library(AHMbook)
library(nimble)
library(parallel)
library(MCMCvis)


frogs <- nimbleCode({
  
  # Priors
  #-- Not infected
  for(t in 1:(nyears - 1)){
    gammaN[t] <- mean.gammaN
    omegaN[t] <- mean.omegaN
    psi_NI[t] <- mean.psi_NI
  }
  
  alpha.lamN  ~ dnorm(0, sd = 10)
  mean.pN     ~ dbeta(1, 1)
  mean.gammaN ~ dnorm(0, sd = 10)
  mean.omegaN ~ dbeta(1, 1)
  mean.psi_NI ~ dbeta(1, 1)
  
  #-- Infeted
  for(t in 1:(nyears-1)){
    gammaI[t] <- mean.gammaI
    omegaI[t] <- mean.omegaI
    psi_IN[t] <- mean.psi_IN
  }
  
  alpha.lamI  ~ dnorm(0, sd = 10)
  mean.pI     ~ dbeta(1, 1)
  mean.gammaI ~ dnorm(0, sd = 10)
  mean.omegaI ~ dbeta(1, 1)
  mean.psi_IN ~ dbeta(1, 1)
  
  #-- Detection
  for(t in 1:nyears){
    pN[t] <- mean.pN
    pI[t] <- mean.pI
  }
  
  # Ecologicacl state model
  #-- First year
  for(i in 1:nsites){
    #-- Not infected
    NN[i, 1] ~ dpois(lambdaN[i])
    log(lambdaN[i]) <- alpha.lamN
    #-- Infected
    NI[i, 1] ~ dpois(lambdaI[i])
    log(lambdaI[i]) <- alpha.lamI
  }
  
  #-- Second and subsequent years
  for(t in 2:nyears){
    for(i in 1:nsites){
      #-- Not infected
      SN[i, t] ~ dbin( omegaN[t-1], NN[i, t-1]) # total survivors
      TN[i, t] ~ dbin( psi_NI[t-1], SN[i, t]) # survive, become infected
      GN[i, t] ~ dpois(GaN[i, t])
      log(GaN[i, t]) <- gammaN[t-1]
      #-- Infected
      SI[i, t] ~ dbin( omegaI[t-1], NI[i, t-1] ) # infecteds who survive
      TI[i, t] ~ dbin( psi_NI[t-1], SI[i, t] )   # survive but become infected
      GI[i, t] ~ dpois( GaI[i, t] )
      log(GaI[i, t]) <- gammaI[t-1]
      # subsequent population size; add everything up and subtract out-transitions
      NN[i, t] <- SN[i, t] - TN[i, t] + GN[i, t] + TI[i, t]
      NI[i, t] <- SI[i, t] - TI[i, t] + GI[i, t] + TN[i, t]
    }
  }
  
  # Observation model
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      for(t in 1:nyears){
        yN[i, j, t] ~ dbin(pN[t], NN[i, t] )
        yI[i, j, t] ~ dbin(pI[t], NI[i, t] )
        # predictions of new data
        y.newN[i, j, t] ~ dbin(pN[t], NN[i, t])
        y.newI[i, j, t] ~ dbin(pI[t], NI[i, t])
      }
    }
  }

  # Bayesian p-values
  # for(t in 1:nyears){
  #   for(i in 1:nsites){
  #     evalN[i, t] <- pN[t] * NN[i, t]
  #     evalI[i, t] <- pI[t] * NI[i, t]
  #     for(j in 1:nsurveys){
  #       EN[i, j, t] <- pow(( yN[i, j, t] - evalN[i, t]), 2) / ( evalN[i, t] + 0.5 )
  #       EI[i, j, t] <- pow(( yI[i, j, t] - evalI[i, t]), 2) / ( evalI[i, t] + 0.5 )
  #       E.newN[i, j, t] <- pow(( y.newN[i, j, t] - evalN[i, t]), 2) / (evalN[i, t] + 0.5 )
  #       E.newI[i, j, t] <- pow(( y.newI[i, j, t] - evalI[i, t]), 2) / (evalI[i, t] + 0.5 )
  #     }
  #   }
  # }
  # 
  # fitN <- sum(EN[1:nsites, 1:nsurveys, 1:nyears])
  # fitN.new <- sum(E.newN[1:nsites, 1:nsurveys, 1:nyears])
  # fitI <- sum(EI[1:nsites, 1:nsurveys, 1:nyears])
  # fitI.new <- sum(E.newI[1:nsites, 1:nsurveys, 1:nyears])
  # 
  
  })

set.seed(2019)
sodata <- simFrogDisease(
  nsites = 100, 
  nyears = 3, 
  nsurveys = 3, 
  alpha.lam = 3, 
  omega = c(0.9, 0.7), 
  gamma = c(2, 1), 
  p = c(0.8, 0.8, 0.8), 
  recovery = 0.1, 
  infection = 0.1
)

data <- list(
  yN = sodata$yN,
  yI = sodata$yI
)

constants <- list(
  nsites = dim(sodata$yN)[1],
  nsurveys = dim(sodata$yN)[2],
  nyears = dim(sodata$yN)[3]
)

inits <- function(){
  list(
    alpha.lamN = runif(1, 2, 3), 
    mean.pN = runif(1, 0.9, 1),
    mean.omegaN = runif(1, 0.7, 1), 
    mean.gammaN = runif(1, 2, 3), 
    mean.psi_NI = runif(1, 0, 0.3), 
    alpha.lamI = runif(1, 2, 3), 
    mean.pI = runif(1, 0.9, 1), 
    mean.omegaI = runif(1, 0.7, 1), 
    mean.gammaI = runif(1, 2, 3), 
    mean.psi_IN = runif(1, 0, 0.3)
  )
}

keepers <- c(
  "alpha.lamN", 
  "alpha.lamI", 
  "mean.pN", 
  "mean.pI", 
  "mean.omegaN", 
  "mean.omegaI", 
  "mean.gammaN", 
  "mean.gammaI", 
  "mean.psi_NI", 
  "mean.psi_IN"#, 
  # "fitN", 
  # "fitN.new", 
  # "fitI", 
  # "fitI.new"
)


ni <- 60000
nb <- 40000
nt <- 10
nc <- 3

# set.seed(1

# m <- nimbleModel(code = frogs,
#                  constants = constants,
#                  data = data,
#                  inits = inits())
# 
# m$initializeInfo()


start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("frogs",
                              "inits",
                              "data",
                              "constants",
                              "keepers",
                              "nb",
                              "ni",
                              "nt"))

for(j in seq_along(cl)){
  set.seed(127)
  init <- inits()
  clusterExport(cl[j], "init")
}

out <- clusterEvalQ(cl, {
  library(nimble)
  
  m <- nimbleModel(code = frogs,
                   name = "frogs",
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

# major convergence issues 
# MCMCsummary(out, params = "alpha.lamN")
