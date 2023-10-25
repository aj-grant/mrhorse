libs = .libPaths()

library(MendelianRandomization, lib.loc = libs)
library(rjags, lib.loc = libs)
library(tidyverse, lib.loc = libs)
library(BWMR, lib.loc = libs)
library(parallel, lib.loc = libs)

source("sims_functions.R")

n = 50000
p = 100000
m = 20
hx = 0.3
hy = 0.3
rho = 0.6

sim_dat = function(p, iv, n, m, theta, hx, hy){
  maf = runif(p, 0.1, 0.3)
  pk = c((1-iv)*m/p, iv*m/p, 1 - m/p)
  vx = c(hx/m, 2*hx/(3*m), 0)
  vd = c(0, hx/(3*m), 0)
  va = c(0, (hy - (theta^2 + 1/3)*hx)/m, 0)
  
  B = sapply(1:3, function(k){rnorm(p, 0, sqrt(vx[k]))})
  D = sapply(1:3, function(k){rnorm(p, 0, sqrt(vd[k]))})
  A = sapply(1:3, function(k){rnorm(p, 0, sqrt(va[k]))})
  P = rmultinom(p, 1, pk)
  
  bt = apply(B * t(P), 1, sum)
  dt = apply(D * t(P), 1, sum)
  at = apply(A * t(P), 1, sum)
  
  bx = (bt + dt) / sqrt(2*maf*(1-maf))
  a = (at + dt) / sqrt(2*maf*(1-maf))
  
  by = theta * bx + a
  s = 1 / sqrt(2 * n * maf * (1-maf))
  
  bxhat = rnorm(p, bx, s)
  byhat = rnorm(p, by, s)
  
  sig1 = which(2*(1-pnorm(abs(bxhat/s))) < 0.001)
  sig2 = which(2*(1-pnorm(abs(bxhat/s))) < 5e-8)
  return(list(bxhat_all = bxhat, byhat_all = byhat, se_all = s, sig1 = sig1, sig2 = sig2))
}

sim_res = function(D, p, n, no_ini, variable.names){
  mrob = mr_input(bx = D$bxhat_all[D$sig2], bxse = D$se_all[D$sig2], by = D$byhat_all[D$sig2], byse = D$se_all[D$sig2])
  mrest_ivw = mr_ivw(mrob)
  mrest_median = mr_median(mrob)
  
  A = BWMR2(gammahat = D$bxhat_all[D$sig2], Gammahat = D$byhat_all[D$sig2], sigmaX = D$se_all[D$sig2], sigmaY = D$se_all[D$sig2])
  
  mrest_cml = mr_cML(mrob, n = n, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, n = n)
  
  c(mrest_ivw$Estimate, mrest_ivw$StdError, mrest_ivw$CILower, mrest_ivw$CIUpper,
    mrest_median$Estimate, mrest_median$StdError, mrest_median$CILower, mrest_median$CIUpper,
    A$beta, A$se_beta, A$beta - qnorm(0.975, 0, 1) * A$se_beta, A$beta + qnorm(0.975, 0, 1) * A$se_beta,
    mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
  )
}

sim_res_horse = function(D, p, n, no_ini, variable.names){
  mrjags = jags.model(file = "mrbayes-cor-model.txt",
                      data = list(by = D$byhat_all[D$sig2], bx = D$bxhat_all[D$sig2], sy = D$se_all[D$sig2], sx = D$se_all[D$sig2], N = length(D$sig2)),
                      n.chains = no_ini,
                      quiet = TRUE)
  update(mrjags, n.iter = 10000, progress.bar = "none")
  mrjags.coda <- coda.samples(mrjags,
                              variable.names=variable.names,
                              n.iter=10000,
                              quiet = TRUE, progress.bar = "none")
  
  
  c(summary(mrjags.coda[, "theta"])$statistics[1], summary(mrjags.coda[, "theta"])$statistics[2],
    summary(mrjags.coda[, "theta"])$quantiles[1], summary(mrjags.coda[, "theta"])$quantiles[5]
  )
}

sim_res_cause = function(D, p){
  X1 = data.frame("SNPX" = sprintf("SNP_%i", 1:p), "bxhat" = D$bxhat_all, "sebx" = D$se_all, "EAX" = rep("C", p), "OAX" = rep("T", p))
  X2 = data.frame("SNPY" = sprintf("SNP_%i", 1:p), "byhat" = D$byhat_all, "seby" = D$se_all, "EAY" = rep("C", p), "OAY" = rep("T", p))
  X0 = gwas_merge(X1, X2, snp_name_cols = c("SNPX", "SNPY"), beta_hat_cols = c("bxhat", "byhat"), se_cols = c("sebx", "seby"),
                  A1_cols = c("EAX", "EAY"), A2_cols = c("OAX", "OAY"))
  params = est_cause_params(X0, X0$snp)
  topvars = sprintf("SNP_%i", 1:p)[D$sig1]
  mr_cause = cause(X = X0, variants = topvars, param_ests = params, force = TRUE)
  
  c(summary(mr_cause, ci_size = 0.95)[1]$quants[[2]][, 1], pnorm(mr_cause$elpd[3, 5])
  )
}

cl = makeCluster(16)
clusterExport(cl, c('libs'))
clusterEvalQ(cl, library(myfun, lib.loc = libs))
clusterEvalQ(cl, library(MendelianRandomization, lib.loc = libs))
clusterEvalQ(cl, library(tidyverse, lib.loc = libs))
clusterEvalQ(cl, library(BWMR, lib.loc = libs))
clusterEvalQ(cl, library(rjags, lib.loc = libs))
clusterEvalQ(cl, library(cause, lib.loc = libs))

M = 200
clusterExport(cl, c('M', 'p', 'n', 'm', 'hx', 'hy', 'BWMR2', 'sstat', 'sim_dat'))

clusterSetRNGStream(cl, 20220616)
D_gwas_20 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.2, n, m, theta = 0.1, hx, hy)
  return(A)
})

clusterSetRNGStream(cl, 20220616)
D_gwas_40 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.4, n, m, theta = 0.1, hx, hy)
  return(A)
})

clusterSetRNGStream(cl, 20220616)
D_gwas_60 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.6, n, m, theta = 0.1, hx, hy)
  return(A)
})

clusterSetRNGStream(cl, 20220616)
D_gwas_null_20 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.2, n, m, theta = 0, hx, hy)
  return(A)
})

clusterSetRNGStream(cl, 20220616)
D_gwas_null_40 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.4, n, m, theta = 0, hx, hy)
  return(A)
})

clusterSetRNGStream(cl, 20220616)
D_gwas_null_60 = parLapply(cl, 1:M, function(i){
  A = sim_dat(p, iv = 0.6, n, m, theta = 0, hx, hy)
  return(A)
})
save('D_gwas_20', 'D_gwas_40', 'D_gwas_60', 'D_gwas_null_20', 'D_gwas_null_40', 'D_gwas_null_60', file = 'sims_dat_gwas.R')
no_ini = 3
variable.names = "theta"
clusterExport(cl, c('M', 'no_ini', 'variable.names', 'sim_res', 'sim_res_horse', 'sim_res_cause'))
clusterExport(cl, c('D_gwas_20', 'D_gwas_40', 'D_gwas_60', 'D_gwas_null_20', 'D_gwas_null_40', 'D_gwas_null_60'))

clusterSetRNGStream(cl, 20220617)
R_gwas_20 = parSapply(cl, 1:M, function(i){
  D = D_gwas_20[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_40 = parSapply(cl, 1:M, function(i){
  D = D_gwas_40[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_60 = parSapply(cl, 1:M, function(i){
  D = D_gwas_60[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_20 = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_20[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_40 = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_40[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_60 = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_60[[i]]
  sim_res(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_20_horse = parSapply(cl, 1:M, function(i){
  D = D_gwas_20[[i]]
  sim_res_horse(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_40_horse = parSapply(cl, 1:M, function(i){
  D = D_gwas_40[[i]]
  sim_res_horse(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_60_horse = parSapply(cl, 1:M, function(i){
  D = D_gwas_60[[i]]
  sim_res_horse(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_20_horse = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_20[[i]]
  sim_res_horse(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_40_horse = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_40[[i]]
  sim_res_horse(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_60_horse = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_60[[i]]
  sim_res_horse(D, p, n, no_ini, variable.names)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_20_cause = parSapply(cl, 1:M, function(i){
  D = D_gwas_20[[i]]
  sim_res_cause(D, p)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_40_cause = parSapply(cl, 1:M, function(i){
  D = D_gwas_40[[i]]
  sim_res_cause(D, p)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_60_cause = parSapply(cl, 1:M, function(i){
  D = D_gwas_60[[i]]
  sim_res_cause(D, p)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_20_cause = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_20[[i]]
  sim_res_cause(D, p)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_40_cause = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_40[[i]]
  sim_res_cause(D, p)
})

clusterSetRNGStream(cl, 20220617)
R_gwas_null_60_cause = parSapply(cl, 1:M, function(i){
  D = D_gwas_null_60[[i]]
  sim_res_cause(D, p)
})
stopCluster(cl)
