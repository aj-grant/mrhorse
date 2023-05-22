mr_horse = function(D, no_ini, variable.names, n.iter = 10000, n.burnin = 10000){
  mrjags = jags.model(file = "mrhorse-cor-model.txt",
                      data = list(by = D@betaY, bx = D@betaX, sy = D@betaYse, sx = D@betaXse, N = length(D@betaY)),
                      n.chains = no_ini,
                      quiet = TRUE)
  update(mrjags, n.iter = n.burnin, progress.bar = "none")
  mrjags.coda <- coda.samples(mrjags,
                              variable.names=variable.names,
                              n.iter=n.iter,
                              quiet = TRUE, progress.bar = "none")
}

mr_horse_summary = function(mr_horse.coda){
  c(summary(mr_horse.coda[, "theta"])$statistics[1], summary(mr_horse.coda[, "theta"])$statistics[2],
    summary(mr_horse.coda[, "theta"])$quantiles[1], summary(mr_horse.coda[, "theta"])$quantiles[5],
    gelman.diag(mr_horse.coda, confidence = 0.8)$psrf[1])
}

mr_cause_run = function(Xgwas, Ygwas, sig_snps){
  X1 = data.frame("SNPX" = Xgwas$SNP, "bxhat" = Xgwas$beta, "sebx" = Xgwas$se, "EAX" = Xgwas$EA, "OAX" = Xgwas$NEA, "p" = Xgwas$p)
  X2 = data.frame("SNPY" = Ygwas$SNP, "byhat" = Ygwas$beta, "seby" = Ygwas$se, "EAY" = Ygwas$EA, "OAY" = Ygwas$NEA)
  X0 = gwas_merge(X1, X2, snp_name_cols = c("SNPX", "SNPY"), beta_hat_cols = c("bxhat", "byhat"), se_cols = c("sebx", "seby"),
                  A1_cols = c("EAX", "EAY"), A2_cols = c("OAX", "OAY"))
  set.seed(100)
  varlist = with(X0, sample(snp, size = 1000000, replace = FALSE))
  params = est_cause_params(X0, varlist)
  X0_p = inner_join(X0, data.frame("snp" = sig_snps$SNP, "chr" = sig_snps$chr, "pos" = sig_snps$pos, "EA" = sig_snps$EA,
                                   "NEA" = sig_snps$NEA, "p" = sig_snps$p), by = "snp")
  X0_clump = clump_data(format_data(X0_p, chr_col = "chr", pos_col = "pos", snp_col = "snp",
                                    effect_allele_col = "EA", other_allele_col = "NEA", pval_col = "p",
                                    type = "exposure"),
                        clump_r2 = 0.01, clump_p1 = 0.001)
  top_vars = X0_clump$SNP
  cause_est = cause(X = X0, variants = top_vars, param_ests = params, force = TRUE)
}
