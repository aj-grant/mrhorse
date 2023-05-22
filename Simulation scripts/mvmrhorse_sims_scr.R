libs = .libPaths()

library(MendelianRandomization, lib.loc = libs)
library(rjags, lib.loc = libs)
library(tidyverse, lib.loc = libs)
library(GRAPPLE, lib.loc = libs)
library(parallel, lib.loc = libs)

sim_dat = function(p, K, iv, n, theta, balanced_pl = T){
  bx = sapply(1:K, function(k){runif(p, 0.05, 0.15) * sample(c(-1, 1), p, replace = T)})
  if (balanced_pl == T){
    a = runif(p, -0.2, 0.2)
  } else {
    a = runif(p, -0.1, 0.3)
  }
  d = runif(p, -0.1, 0.1)
  valid_gv = sample(p, (1-iv) * p)
  a[valid_gv] = 0
  d[valid_gv] = 0
  maf = runif(p, 0.1, 0.3)
  G = sapply(1:p, function(j){rbinom(2*n, 2, maf[j])})
  U = rnorm(2*n, 0, 1) + G %*% d
  ex = mv_rnorm(2*n, c(0, 0), matrix(c(1, 0.3, 0.3, 1), nrow = 2))
  X = G %*% bx + U %*% rep(1, K) + ex
  Y = X %*% theta + G %*% a + U + rnorm(2*n, 0, 1)
  
  byhat = seby = vector(length = p)
  bxhat = sebx = matrix(nrow = p, ncol = K)
  sam1 = 1:n
  sam2 = (n+1):(2*n)
  
  for (j in 1:p){
    for (k in 1:K){
      xmod = sstat(X[sam1, k], G[sam1, j])
      bxhat[j, k] = xmod$bhat
      sebx[j, k] = xmod$se
    }
    ymod = sstat(Y[sam2], G[sam2, j])
    byhat[j] = ymod$bhat
    seby[j] = ymod$se
  }
  return(list(bxhat = bxhat, sebx = sebx, byhat = byhat, seby = seby, valid_gv = valid_gv))
}

sim_res = function(D, no_ini, variable.names){
  p = dim(D$bxhat)[1]
  K = dim(D$bxhat)[2]
  mrob = mr_mvinput(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrob_oracle = mr_mvinput(bx = D$bxhat[D$valid_gv, ], bxse = D$sebx[D$valid_gv, ], by = D$byhat[D$valid_gv], byse = D$seby[D$valid_gv])
  mrest_ivw = mr_mvivw(mrob)
  mrest_ivw_oracle = mr_mvivw(mrob_oracle)
  mrest_median = mr_mvmedian(mrob)
  
  grapple_dat = data.frame("SNP" = seq(1:p), "gamma_exp1" = D$bxhat[, 1], "gamma_exp2" = D$bxhat[, 2],
                           "se_exp1" = D$sebx[, 1], "se_exp2" = D$sebx[, 2],
                           "gamma_out1" = D$byhat, "se_out1" = D$seby)
  grapple_est = grappleRobustEst(grapple_dat, plot.it = FALSE)
  
  Tx = matrix(nrow = K, ncol = p*K)
  for (j in 1:p){
    Tx[, ((j-1)*K+1):(j*K)] = diag(1 / D$sebx[j, ]^2)
  }
  mrjags = jags.model(file = "mvmrhorse-cor-model.txt",
                      data = list(by = D$byhat, bx = D$bxhat, sy = D$seby, Tx = Tx, N = p, K = K, R = diag(K)),
                      n.chains = no_ini,
                      quiet = TRUE)
  update(mrjags, n.iter = 10000, progress.bar = "none")
  mrjags.coda <- coda.samples(mrjags,
                              variable.names=variable.names,
                              n.iter=10000,
                              quiet = TRUE, progress.bar = "none")
  
  q = qnorm(0.975, 0, 1)
  c(mrest_ivw$Estimate, mrest_ivw$StdError, mrest_ivw$CILower, mrest_ivw$CIUpper,
    mrest_ivw_oracle$Estimate, mrest_ivw_oracle$StdError, mrest_ivw_oracle$CILower, mrest_ivw_oracle$CIUpper,
    mrest_median$Estimate, mrest_median$StdError, mrest_median$CILower, mrest_median$CIUpper,
    grapple_est$beta.hat, sqrt(diag(grapple_est$beta.var)),
    grapple_est$beta.hat - q * sqrt(diag(grapple_est$beta.var)),
    grapple_est$beta.hat + q * sqrt(diag(grapple_est$beta.var)),
    summary(mrjags.coda)$statistics[, 1], summary(mrjags.coda)$statistics[, 2],
    summary(mrjags.coda)$quantiles[, 1], summary(mrjags.coda)$quantiles[, 5]
  )
}

cl = makeCluster(24)
clusterExport(cl, c('libs'))
clusterEvalQ(cl, library(MendelianRandomization, lib.loc = libs))
clusterEvalQ(cl, library(GRAPPLE, lib.loc = libs))
clusterEvalQ(cl, library(rjags, lib.loc = libs))

tic()
M = 200
p = 100
K = 2
n = 50000
clusterExport(cl, c('M', 'p', 'K', 'n', 'sstat', 'mv_rnorm', 'sim_dat'))

clusterSetRNGStream(cl, 20220506)
D_20 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.2, n, theta = c(0.1, 0.1))
  return(A)
})

D_20_null = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.2, n, theta = c(0, 0.1))
  return(A)
})

D_20_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.2, n, theta = c(0.1, 0.1), balanced_pl = F)
  return(A)
})

D_20_null_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.2, n, theta = c(0, 0.1), balanced_pl = F)
  return(A)
})

clusterSetRNGStream(cl, 20220506)
D_40 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.4, n, theta = c(0.1, 0.1))
  return(A)
})

D_40_null = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.4, n, theta = c(0, 0.1))
  return(A)
})

D_40_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.4, n, theta = c(0.1, 0.1), balanced_pl = F)
  return(A)
})

D_40_null_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.4, n, theta = c(0, 0.1), balanced_pl = F)
  return(A)
})

clusterSetRNGStream(cl, 20220506)
D_60 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.6, n, theta = c(0.1, 0.1))
  return(A)
})

D_60_null = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.6, n, theta = c(0, 0.1))
  return(A)
})

D_60_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.6, n, theta = c(0.1, 0.1), balanced_pl = F)
  return(A)
})

D_60_null_dir = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, K, iv = 0.6, n, theta = c(0, 0.1), balanced_pl = F)
  return(A)
})

########################################################################################################################

M = 200
no_ini = 3
variable.names = "theta"
clusterExport(cl, c('M', 'no_ini', 'variable.names', 'sim_res'))
clusterExport(cl, c('D_20', 'D_20_null', 'D_20_dir', 'D_20_null_dir', 'D_40', 'D_40_null', 'D_40_dir', 'D_40_null_dir',
              'D_60', 'D_60_null', 'D_60_dir', 'D_60_null_dir'))

clusterSetRNGStream(cl, 20220506)
R_20 = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_20_null = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_20_dir = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_20_null_dir = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_40 = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_40_null = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_40_dir = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_40_null_dir = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_60 = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_60_null = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_60_dir = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  sim_res(D, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220506)
R_60_null_dir = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  sim_res(D, no_ini, variable.names)
})

stopCluster(cl)
