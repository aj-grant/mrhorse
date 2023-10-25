libs = .libPaths()

library(MendelianRandomization, lib.loc = libs)
library(rjags, lib.loc = libs)
library(tidyverse, lib.loc = libs)
library(BWMR, lib.loc = libs)
library(parallel, lib.loc = libs)

source("sims_functions.R")

sim_dat = function(p, iv, n, theta, balanced_pl = T){
  bx = runif(p, 0.05, 0.15) * sample(c(-1, 1), p, replace = T)
  if (balanced_pl == T){
    a = runif(p, -0.2, 0.2)
  } else {
    a = runif(p, -0.1, 0.3)
  }
  d = runif(p, -0.1, 0.1)
  valid_gv = sample(p, (1-iv) * p)
  a[valid_gv] = 0
  d[valid_gv] = 0
  Z = rep(1, p)
  Z[valid_gv] = 0
  maf = runif(p, 0.1, 0.3)
  G = sapply(1:p, function(j){rbinom(2*n, 2, maf[j])})
  U = rnorm(2*n, 0, 1) + G %*% (d*Z)
  X = G %*% (bx*(1-Z)) + U + rnorm(2*n, 0, 1)
  Y = theta * X + G %*% a + U + rnorm(2*n, 0, 1)
  
  bxhat = byhat = sebx = seby = vector(length = p)
  sam1 = 1:n
  sam2 = (n+1):(2*n)
  
  for (j in 1:p){
    xmod = sstat(X[sam1], G[sam1, j])
    bxhat[j] = xmod$bhat
    sebx[j] = xmod$se
    ymod = sstat(Y[sam2], G[sam2, j])
    byhat[j] = ymod$bhat
    seby[j] = ymod$se
  }
  r2 = summary(lm(X[sam1]~G[sam1, ]))$r.squared
  return(list(bxhat = bxhat, sebx = sebx, byhat = byhat, seby = seby, r2 = r2, valid_gv = valid_gv))
}

sim_res = function(D, p, n, no_ini, variable.names){
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrob_oracle = mr_input(bx = D$bxhat[D$valid_gv], bxse = D$sebx[D$valid_gv], by = D$byhat[D$valid_gv], byse = D$seby[D$valid_gv])
  mrest_ivw = mr_ivw(mrob)
  mrest_ivw_oracle = mr_ivw(mrob_oracle)
  mrest_median = mr_median(mrob)
  
  A = BWMR2(gammahat = D$bxhat, Gammahat = D$byhat, sigmaX = D$sebx, sigmaY = D$seby)
  
  mrest_cml = mr_cML(mrob, n = n, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, n = n)
  
  mrjags = jags.model(file = "mrhorse-cor-model.txt",
                      data = list(by = D$byhat, bx = D$bxhat, sy = D$seby, sx = D$sebx, N = p),
                      n.chains = no_ini,
                      quiet = TRUE)
  update(mrjags, n.iter = 10000, progress.bar = "none")
  mrjags.coda <- coda.samples(mrjags,
                              variable.names=variable.names,
                              n.iter=10000,
                              quiet = TRUE, progress.bar = "none")
  
  c(mrest_ivw$Estimate, mrest_ivw$StdError, mrest_ivw$CILower, mrest_ivw$CIUpper,
    mrest_ivw_oracle$Estimate, mrest_ivw_oracle$StdError, mrest_ivw_oracle$CILower, mrest_ivw_oracle$CIUpper,
    mrest_median$Estimate, mrest_median$StdError, mrest_median$CILower, mrest_median$CIUpper,
    A$beta, A$se_beta, A$beta - qnorm(0.975, 0, 1) * A$se_beta, A$beta + qnorm(0.975, 0, 1) * A$se_beta,
    mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper,
    summary(mrjags.coda[, "theta"])$statistics[1], summary(mrjags.coda[, "theta"])$statistics[2],
    summary(mrjags.coda[, "theta"])$quantiles[1], summary(mrjags.coda[, "theta"])$quantiles[5]
  )
}

cl = makeCluster(24)
clusterExport(cl, c('libs'))
clusterEvalQ(cl, library(MendelianRandomization, lib.loc = libs))
clusterEvalQ(cl, library(BWMR, lib.loc = libs))
clusterEvalQ(cl, library(tidyverse, lib.loc = libs))
clusterEvalQ(cl, library(rjags, lib.loc = libs))

M = 1000
clusterSetRNGStream(cl, 20220305)
p = 100
n = 50000
clusterExport(cl, c('M', 'p', 'n', 'BWMR2', 'sstat', 'sim_dat'))

D_20 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.2, n, theta = 0.1)
  return(A)
})

D_20_null = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.2, n, theta = 0)
  return(A)
})

D_20_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.2, n, theta = 0.1, balanced_pl = F)
  return(A)
})

D_20_null_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.2, n, theta = 0, balanced_pl = F)
  return(A)
})

clusterSetRNGStream(cl, 20220305)
D_40 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.4, n, theta = 0.1)
  return(A)
})

D_40_null = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.4, n, theta = 0)
  return(A)
})

D_40_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.4, n, theta = 0.1, balanced_pl = F)
  return(A)
})

D_40_null_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.4, n, theta = 0, balanced_pl = F)
  return(A)
})

clusterSetRNGStream(cl, 20220305)
D_60 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.6, n, theta = 0.1)
  return(A)
})

D_60_null = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.6, n, theta = 0)
  return(A)
})

D_60_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.6, n, theta = 0.1, balanced_pl = F)
  return(A)
})

D_60_null_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.6, n, theta = 0, balanced_pl = F)
  return(A)
})

########################################################################################################################

no_ini = 3
variable.names = "theta"
clusterExport(cl, c('M', 'no_ini', 'variable.names', 'sim_res'))
clusterExport(cl, c('D_20', 'D_20_null', 'D_20_dir', 'D_20_null_dir', 'D_40', 'D_40_null', 'D_40_dir', 'D_40_null_dir',
                    'D_60', 'D_60_null', 'D_60_dir', 'D_60_null_dir'))

clusterSetRNGStream(cl, 20220306)
R_20 = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_20_null = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_20_dir = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_dir = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_40 = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_40_null = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_40_dir = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_dir = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_60 = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_60_null = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_60_dir = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_dir = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

stopCluster(cl)

#MR-Cue
sim_res_cue = function(D){
  mrest_cue = MRCUEIndep(gammah = D$bxhat, Gammah = D$byhat, se1 = D$sebx, se2 = D$seby, rho = 0)
  c(mrest_cue$beta.hat, mrest_cue$beta.se, mrest_cue$beta.hat - qnorm(0.975, 0, 1)*mrest_cue$beta.se,
    mrest_cue$beta.hat + qnorm(0.975, 0, 1)*mrest_cue$beta.se)
}

cl = makeCluster(8)
clusterEvalQ(cl, library(MR.CUE))

clusterExport(cl, c('M', 'sim_res_cue'))
clusterExport(cl, c('D_20', 'D_20_null', 'D_20_dir', 'D_20_null_dir', 'D_40', 'D_40_null', 'D_40_dir', 'D_40_null_dir',
                    'D_60', 'D_60_null', 'D_60_dir', 'D_60_null_dir'))

clusterSetRNGStream(cl, 20220306)
R_20_cue = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_cue = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_20_dir_cue = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_dir_cue = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_40_cue = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_cue = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_40_dir_cue = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_dir_cue = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_60_cue = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_cue = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_60_dir_cue = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  sim_res_cue(D)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_dir_cue = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  sim_res_cue(D)
})

stopCluster(cl)
