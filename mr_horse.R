library(R2jags)

mr_horse_model = function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1/(sy[i] * sy[i]))
    mu[i] = theta * bx0[i] + alpha[i]
    bx[i] ~ dnorm(bx0[i], 1 / (sx[i] * sx[i]))
    
    bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0))
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] -1
    
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }
  
  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)
  
  vx0 ~ dnorm(0, 1);T(0, )
  mx0 ~ dnorm(0, 1)
  
  theta ~ dunif(-10, 10)
}

mvmr_horse_model = function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1 / (sy[i] * sy[i]))
    mu[i] = inprod(bx0[i, 1:K], theta) + alpha[i]
    bx[i,1:K] ~ dmnorm(bx0[i,1:K], Tx[1:K, ((i-1)*K+1):(i*K)])
    
    kappa[i] = (rho[i]^2 / (1 + K*rho[i]^2))
    bx0[i,1:K] ~ dmnorm(mx + sx0 * rho[i] * alpha[i] / (phi[i] * tau), A - kappa[i] * B)
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] -1
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }
  
  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)
  
  mx ~ dmnorm(rep(0, K), R[,])
  
  for (k in 1:K){
    vx0[k] ~ dnorm(0, 1);T(0, )
    sx0[k] = sqrt(vx0[k])
    theta[k] ~ dunif(-10, 10)
    for (j in 1:K){
      A[j, k] = ifelse(j==k, 1/vx0[j], 0)
      B[j, k] = 1 / (sx0[j] * sx0[k])
    }
  }
}

mr_horse = function(D, no_ini = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000){
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }
  jags_fit = jags(data = list(by = D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse, N = length(D$betaY)),
                  parameters.to.save = variable.names,
                  n.chains = no_ini,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  model.file = mr_horse_model)
  mr.coda = as.mcmc(jags_fit)
  mr_estimate = data.frame("Estimate" = round(unname(summary(mr.coda[, "theta"])$statistics[1]), 3),
                           "SD" = round(unname(summary(mr.coda[, "theta"])$statistics[2]), 3),
                           "2.5% quantile" = round(unname(summary(mr.coda[, "theta"])$quantiles[1]), 3),
                           "97.5% quantile" = round(unname(summary(mr.coda[, "theta"])$quantiles[5]), 3),
                           "Rhat" = round(unname(gelman.diag(mr.coda)$psrf[1]), 3))
  names(mr_estimate) = c("Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}

mvmr_horse = function(D, no_ini = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000){
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }
  
  p = dim(D)[1]
  K = sum(sapply(1:dim(D)[2], function(j){substr(names(D)[j], 1, 5)=="betaX"}))/2
  
  Bx = D[, sprintf("betaX%i", 1:K)]
  Sx = D[, sprintf("betaX%ise", 1:K)]
  Tx = matrix(nrow = K, ncol = p*K)
  for (j in 1:p){
    Tx[, ((j-1)*K+1):(j*K)] = diag(1 / Sx[j, ]^2)
  }
  jags_fit = jags(data = list(by = D$betaY, bx = Bx, sy = D$betaYse, Tx = Tx, N = p, K = K, R = diag(K)),
                  parameters.to.save = variable.names,
                  n.chains = no_ini,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  model.file = mvmr_horse_model)
  mr.coda = as.mcmc(jags_fit)
  mr_estimate = data.frame("Parameter" = sprintf("theta[%i]", 1:K),
                           "Estimate" = round(unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 1]), 3),
                           "SD" = round(unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 2]), 3),
                           "2.5% quantile" = round(unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 1]), 3),
                           "97.5% quantile" = round(unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 5]), 3),
                           "Rhat" = round(unname(gelman.diag(mr.coda)$psrf[sprintf("theta[%i]", 1:K), 1]), 3))
  names(mr_estimate) = c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}
