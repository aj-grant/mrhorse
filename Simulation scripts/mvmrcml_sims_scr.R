libs = .libPaths()

library(parallel, lib.loc = libs)
library(MVMRcML, lib.loc = libs)

sim_res = function(D, Kmax){
  p = dim(D$bxhat)[1]
  K = dim(D$bxhat)[2]
  
  rho_mat = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)
  Sig_inv_l = invcov_mvmr(se_bx = D$sebx, se_by = D$seby, rho_mat = rho_mat)
  mrest_cml = MVmr_cML_DP(b_exp = D$bxhat,
                          b_out = as.matrix(D$byhat),
                          se_bx = D$sebx,
                          Sig_inv_l = Sig_inv_l,
                          n = 50000,
                          num_pert = 100,
                          K_vec = 0:p)
  c(mrest_cml$BIC_DP_theta, mrest_cml$BIC_DP_se,
    mrest_cml$BIC_theta-qnorm(0.975)*mrest_cml$BIC_DP_se, mrest_cml$BIC_theta+qnorm(0.975)*mrest_cml$BIC_DP_se,
    mrest_cml$eff_DP_B)
}

cl = makeCluster(24)
clusterExport(cl, c('libs'))
clusterEvalQ(cl, library(MVMRcML, lib.loc = libs))

########################################################################################################################

M = 200
clusterExport(cl, c('M', 'sim_res'))
clusterExport(cl, c('D_20', 'D_20_null', 'D_20_dir', 'D_20_null_dir', 'D_40', 'D_40_null', 'D_40_dir', 'D_40_null_dir',
              'D_60', 'D_60_null', 'D_60_dir', 'D_60_null_dir'))

clusterSetRNGStream(cl, 20220506)
R_20 = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_20_null = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_20_dir = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_20_null_dir = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_40 = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_40_null = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_40_dir = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_40_null_dir = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_60 = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_60_null = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_60_dir = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  sim_res(D)
})

clusterSetRNGStream(cl, 20220506)
R_60_null_dir = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  sim_res(D)
})
stopCluster(cl)
