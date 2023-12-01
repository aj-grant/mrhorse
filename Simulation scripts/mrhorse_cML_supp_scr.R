library(MendelianRandomization)
library(parallel)

p = 100
n = 50000
nx = 50000

M = 1000

cl = makeCluster(8)
clusterEvalQ(cl, library(MendelianRandomization))

clusterExport(cl, c('M', 'nx'))
clusterExport(cl, c('D_20', 'D_20_null', 'D_20_dir', 'D_20_null_dir',
                    'D_40', 'D_40_null', 'D_40_dir', 'D_40_null_dir',
                    'D_60', 'D_60_null', 'D_60_dir', 'D_60_null_dir'))
################################################################################
#Set K to its true value
################################################################################
K0 = 20
clusterExport(cl, c('K0'))
clusterSetRNGStream(cl, 20220306)
R_20_cML_K = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_cML_K = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_dir_cML_K = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_dir_cML_K = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

K0 = 40
clusterExport(cl, c('K0'))
clusterSetRNGStream(cl, 20220306)
R_40_cML_K = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_cML_K = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_dir_cML_K = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_dir_cML_K = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

K0 = 60
clusterExport(cl, c('K0'))
clusterSetRNGStream(cl, 20220306)
R_60_cML_K = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_cML_K = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_dir_cML_K = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_dir_cML_K = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

################################################################################
#Set K to 20
################################################################################
K0 = 20
clusterExport(cl, c('K0'))
clusterSetRNGStream(cl, 20220306)
R_20_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_dir_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_dir_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_dir_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_dir_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_dir_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_dir_cML_K20 = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

################################################################################
#Set K to 40
################################################################################
K0 = 40
clusterExport(cl, c('K0'))
clusterSetRNGStream(cl, 20220306)
R_20_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_dir_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_dir_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_dir_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_dir_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_dir_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_dir_cML_K40 = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})


################################################################################
#Set K to 60
################################################################################
K0 = 60
clusterExport(cl, c('K0'))
clusterSetRNGStream(cl, 20220306)
R_20_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_20[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_20_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_dir_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_20_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_20_null_dir_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_20_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_40[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_40_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_dir_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_40_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_40_null_dir_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_40_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_60[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_60_null[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_dir_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_60_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})

clusterSetRNGStream(cl, 20220306)
R_60_null_dir_cML_K60 = parSapply(cl, 1:M, function(i){
  D = D_60_null_dir[[i]]
  mrob = mr_input(bx = D$bxhat, bxse = D$sebx, by = D$byhat, byse = D$seby)
  mrest_cml = mr_cML(mrob, K = K0, n = nx, DP = FALSE)
  mrest_cml_DP = mr_cML(mrob, K = K0, n = nx, DP = TRUE)
  c(mrest_cml$Estimate, mrest_cml$StdError, mrest_cml$CILower, mrest_cml$CIUpper,
    mrest_cml_DP$Estimate, mrest_cml_DP$StdError, mrest_cml_DP$CILower, mrest_cml_DP$CIUpper)
})
