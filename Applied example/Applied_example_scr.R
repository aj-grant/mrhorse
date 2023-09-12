library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(rjags)
library(cause)

var_rsid_map = read.table(file = 'variants_rsid_map.txt',
                          colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric"),
                          header = TRUE)
names(var_rsid_map) = c("variant", "chr", "pos", "NEA", "EA", "SNP", "minor_allele", "MAF")
var_rsid_map_chrpos = read.table(file = 'var_rsid_map_chrpos.tsv', sep = '\t', colClasses = c("character"),
                                 header = TRUE)
names(var_rsid_map_chrpos) = "CHRPOS"
var_rsid_map = cbind(var_rsid_map, var_rsid_map_chrpos)
rm(var_rsid_map_chrpos)

cad_dat = read.table(file = 'cad_out.txt', colClasses = c("character", "character", "character", "numeric", "numeric", "numeric"),
                     header = TRUE)
names(cad_dat) = c("SNP", "EA_cad", "NEA_cad", "EAF", "beta_cad", "se_cad")
cad_fm = format_data(cad_dat, snp_col = "SNP", beta_col = "beta_cad", se_col = "se_cad",
                    effect_allele_col = "EA_cad", other_allele_col = "NEA_cad", eaf = "EAF",
                    type = "outcome")

t2d_dat = read.table(file = 'T2D_out.txt', colClasses = c("character", "character", "character", "numeric", "numeric", "numeric"),
                     header = TRUE)
names(t2d_dat) = c("CHRPOS", "EA_t2d", "NEA_t2d", "EAF", "beta_t2d", "se_t2d")
t2d_dat = inner_join(var_rsid_map, t2d_dat, by = "CHRPOS")
t2d_dat = t2d_dat[, c("SNP", "EA_t2d", "NEA_t2d", "EAF", "beta_t2d", "se_t2d")]
t2d_fm = format_data(t2d_dat, snp_col = "SNP", beta_col = "beta_t2d", se_col = "se_t2d",
                     effect_allele_col = "EA_t2d", other_allele_col = "NEA_t2d", eaf = "EAF",
                     type = "outcome")

ast_dat = read.table(file = 'asthma_out.txt', colClasses = c("character", "character", "character", "numeric", "numeric", "numeric"),
                     header = TRUE)
names(ast_dat) = c("SNP", "EA_ast", "NEA_ast", "EAF", "beta_ast", "se_ast")
ast_fm = format_data(ast_dat, snp_col = "SNP", beta_col = "beta_ast", se_col = "se_ast",
                     effect_allele_col = "EA_ast", other_allele_col = "NEA_ast", eaf = "EAF",
                     type = "outcome")

str_dat = read.table(file = 'stroke_out.txt', colClasses = c("character", "character", "character", "numeric", "numeric", "numeric"),
                     header = TRUE)
names(str_dat) = c("SNP", "EA_str", "NEA_str", "EAF", "beta_str", "se_str")
str_fm = format_data(str_dat, snp_col = "SNP", beta_col = "beta_str", se_col = "se_str",
                     effect_allele_col = "EA_str", other_allele_col = "NEA_str", eaf = "EAF",
                     type = "outcome")

########################################################################################################################
#Birthweight from Horikoshi
########################################################################################################################
bw_all = read.table(file = "bw_sig_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                            "character", "numeric", "numeric", "numeric", "numeric"),
                    header = FALSE)
names(bw_all) = c("SNP", "chr", "pos", "EA", "NEA", "EAF", "beta_bw", "se_bw", "p_bw")

bw_all_fm = format_data(bw_all, chr_col = "chr", pos_col = "pos", snp_col = "SNP",
                        effect_allele_col = "EA", other_allele_col = "NEA", eaf_col = "EAF",
                        beta_col = "beta_bw", se_col = "se_bw", pval_col = "p_bw",
                        type = "exposure")
bw_clump = clump_data(bw_all_fm, clump_r2 = 0.001)

bw_cad_harmonise = harmonise_data(exposure_dat = bw_clump, outcome_dat = cad_fm, action = 2)
bw_cad = data.frame("SNP" = bw_cad_harmonise$SNP,
                    "beta_bw" = bw_cad_harmonise$beta.exposure, "se_bw" = bw_cad_harmonise$se.exposure,
                    "beta_cad" = bw_cad_harmonise$beta.outcome, "se_cad" = bw_cad_harmonise$se.outcome)

bw_t2d_harmonise = harmonise_data(exposure_dat = bw_clump, outcome_dat = t2d_fm, action = 2)
bw_t2d = data.frame("SNP" = bw_t2d_harmonise$SNP,
                    "beta_bw" = bw_t2d_harmonise$beta.exposure, "se_bw" = bw_t2d_harmonise$se.exposure,
                    "beta_t2d" = bw_t2d_harmonise$beta.outcome, "se_t2d" = bw_t2d_harmonise$se.outcome)

bw_ast_harmonise = harmonise_data(exposure_dat = bw_clump, outcome_dat = ast_fm, action = 2)
bw_ast = data.frame("SNP" = bw_ast_harmonise$SNP,
                    "beta_bw" = bw_ast_harmonise$beta.exposure, "se_bw" = bw_ast_harmonise$se.exposure,
                    "beta_ast" = bw_ast_harmonise$beta.outcome, "se_ast" = bw_ast_harmonise$se.outcome)

rm(bw_all_fm, bw_ast_harmonise, bw_cad_harmonise, bw_t2d_harmonise, bw_clump)

########################################################################################################################
#FG from MAGIC
########################################################################################################################
fg_all = read.table(file = "fg_sig_out.txt", colClasses = c("character", "numeric", "numeric", "character",
                                                              "character", "numeric", "numeric", "numeric", "numeric"),
                     header = FALSE)
names(fg_all) = c("SNP", "chr", "pos", "EA", "NEA", "EAF", "beta_fg", "se_fg", "p_fg")

fg_all_fm = format_data(fg_all, chr_col = "chr", pos_col = "pos", snp_col = "SNP",
                         effect_allele_col = "EA", other_allele_col = "NEA", eaf_col = "EAF",
                         beta_col = "beta_fg", se_col = "se_fg", pval_col = "p_fg",
                         type = "exposure")
fg_clump = clump_data(fg_all_fm, clump_r2 = 0.001)

fg_cad_harmonise = harmonise_data(exposure_dat = fg_clump, outcome_dat = cad_fm, action = 2)
fg_cad = data.frame("SNP" = fg_cad_harmonise$SNP,
                     "beta_fg" = fg_cad_harmonise$beta.exposure, "se_fg" = fg_cad_harmonise$se.exposure,
                     "beta_cad" = fg_cad_harmonise$beta.outcome, "se_cad" = fg_cad_harmonise$se.outcome)

fg_t2d_harmonise = harmonise_data(exposure_dat = fg_clump, outcome_dat = t2d_fm, action = 2)
fg_t2d = data.frame("SNP" = fg_t2d_harmonise$SNP,
                     "beta_fg" = fg_t2d_harmonise$beta.exposure, "se_fg" = fg_t2d_harmonise$se.exposure,
                     "beta_t2d" = fg_t2d_harmonise$beta.outcome, "se_t2d" = fg_t2d_harmonise$se.outcome)

fg_ast_harmonise = harmonise_data(exposure_dat = fg_clump, outcome_dat = ast_fm, action = 2)
fg_ast = data.frame("SNP" = fg_ast_harmonise$SNP,
                    "beta_fg" = fg_ast_harmonise$beta.exposure, "se_fg" = fg_ast_harmonise$se.exposure,
                    "beta_ast" = fg_ast_harmonise$beta.outcome, "se_ast" = fg_ast_harmonise$se.outcome)

rm(fg_all, fg_all_fm, fg_ast_harmonise, fg_cad_harmonise, fg_t2d_harmonise, fg_clump)

########################################################################################################################
#LDL from GLGC
########################################################################################################################
ldl_all = read.table(file = "ldl0_sig_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                              "character", "numeric", "numeric", "numeric", "numeric"),
                     header = FALSE)
names(ldl_all) = c("SNP", "chr", "pos", "NEA", "EA", "EAF", "beta_ldl", "se_ldl", "p_ldl")

ldl_all_fm = format_data(ldl_all, chr_col = "chr", pos_col = "pos", snp_col = "SNP",
                         effect_allele_col = "EA", other_allele_col = "NEA", eaf_col = "EAF",
                         beta_col = "beta_ldl", se_col = "se_ldl", pval_col = "p_ldl",
                         type = "exposure")
ldl_clump = clump_data(ldl_all_fm, clump_r2 = 0.001)

ldl_cad_harmonise = harmonise_data(exposure_dat = ldl_clump, outcome_dat = cad_fm, action = 2)
ldl_cad = data.frame("SNP" = ldl_cad_harmonise$SNP,
                     "beta_ldl" = ldl_cad_harmonise$beta.exposure, "se_ldl" = ldl_cad_harmonise$se.exposure,
                     "beta_cad" = ldl_cad_harmonise$beta.outcome, "se_cad" = ldl_cad_harmonise$se.outcome)

ldl_t2d_harmonise = harmonise_data(exposure_dat = ldl_clump, outcome_dat = t2d_fm, action = 2)
ldl_t2d = data.frame("SNP" = ldl_t2d_harmonise$SNP,
                     "beta_ldl" = ldl_t2d_harmonise$beta.exposure, "se_ldl" = ldl_t2d_harmonise$se.exposure,
                     "beta_t2d" = ldl_t2d_harmonise$beta.outcome, "se_t2d" = ldl_t2d_harmonise$se.outcome)

ldl_ast_harmonise = harmonise_data(exposure_dat = ldl_clump, outcome_dat = ast_fm, action = 2)
ldl_ast = data.frame("SNP" = ldl_ast_harmonise$SNP,
                    "beta_ldl" = ldl_ast_harmonise$beta.exposure, "se_ldl" = ldl_ast_harmonise$se.exposure,
                    "beta_ast" = ldl_ast_harmonise$beta.outcome, "se_ast" = ldl_ast_harmonise$se.outcome)

rm(ldl_all, ldl_all_fm, ldl_ast_harmonise, ldl_cad_harmonise, ldl_t2d_harmonise, ldl_clump)

########################################################################################################################
#TC from GLGC
########################################################################################################################
tc_all = read.table(file = "tc0_sig_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                              "character", "numeric", "numeric", "numeric", "numeric"),
                     header = FALSE)
names(tc_all) = c("SNP", "chr", "pos", "NEA", "EA", "EAF", "beta_tc", "se_tc", "p_tc")

tc_all_fm = format_data(tc_all, chr_col = "chr", pos_col = "pos", snp_col = "SNP",
                         effect_allele_col = "EA", other_allele_col = "NEA", eaf_col = "EAF",
                         beta_col = "beta_tc", se_col = "se_tc", pval_col = "p_tc",
                         type = "exposure")
tc_clump = clump_data(tc_all_fm, clump_r2 = 0.001)

tc_cad_harmonise = harmonise_data(exposure_dat = tc_clump, outcome_dat = cad_fm, action = 2)
tc_cad = data.frame("SNP" = tc_cad_harmonise$SNP,
                     "beta_tc" = tc_cad_harmonise$beta.exposure, "se_tc" = tc_cad_harmonise$se.exposure,
                     "beta_cad" = tc_cad_harmonise$beta.outcome, "se_cad" = tc_cad_harmonise$se.outcome)

tc_t2d_harmonise = harmonise_data(exposure_dat = tc_clump, outcome_dat = t2d_fm, action = 2)
tc_t2d = data.frame("SNP" = tc_t2d_harmonise$SNP,
                     "beta_tc" = tc_t2d_harmonise$beta.exposure, "se_tc" = tc_t2d_harmonise$se.exposure,
                     "beta_t2d" = tc_t2d_harmonise$beta.outcome, "se_t2d" = tc_t2d_harmonise$se.outcome)

tc_ast_harmonise = harmonise_data(exposure_dat = tc_clump, outcome_dat = ast_fm, action = 2)
tc_ast = data.frame("SNP" = tc_ast_harmonise$SNP,
                     "beta_tc" = tc_ast_harmonise$beta.exposure, "se_tc" = tc_ast_harmonise$se.exposure,
                     "beta_ast" = tc_ast_harmonise$beta.outcome, "se_ast" = tc_ast_harmonise$se.outcome)

rm(tc_all, tc_all_fm, tc_ast_harmonise, tc_cad_harmonise, tc_t2d_harmonise, tc_clump)

########################################################################################################################
#IVW estimates
########################################################################################################################
bw_cad_obj = mr_input(bx = bw_cad$beta_bw, bxse = bw_cad$se_bw, by = bw_cad$beta_cad, byse = bw_cad$se_cad)
fg_cad_obj = mr_input(bx = fg_cad$beta_fg, bxse = fg_cad$se_fg, by = fg_cad$beta_cad, byse = fg_cad$se_cad)
ldl_cad_obj = mr_input(bx = ldl_cad$beta_ldl, bxse = ldl_cad$se_ldl, by = ldl_cad$beta_cad, byse = ldl_cad$se_cad)
tc_cad_obj = mr_input(bx = tc_cad$beta_tc, bxse = tc_cad$se_tc, by = tc_cad$beta_cad, byse = tc_cad$se_cad)

bw_cad_ivw = mr_ivw(bw_cad_obj)
fg_cad_ivw = mr_ivw(fg_cad_obj)
ldl_cad_ivw = mr_ivw(ldl_cad_obj)
tc_cad_ivw = mr_ivw(tc_cad_obj)

bw_t2d_obj = mr_input(bx = bw_t2d$beta_bw, bxse = bw_t2d$se_bw, by = bw_t2d$beta_t2d, byse = bw_t2d$se_t2d)
fg_t2d_obj = mr_input(bx = fg_t2d$beta_fg, bxse = fg_t2d$se_fg, by = fg_t2d$beta_t2d, byse = fg_t2d$se_t2d)
ldl_t2d_obj = mr_input(bx = ldl_t2d$beta_ldl, bxse = ldl_t2d$se_ldl, by = ldl_t2d$beta_t2d, byse = ldl_t2d$se_t2d)
tc_t2d_obj = mr_input(bx = tc_t2d$beta_tc, bxse = tc_t2d$se_tc, by = tc_t2d$beta_t2d, byse = tc_t2d$se_t2d)

bw_t2d_ivw = mr_ivw(bw_t2d_obj)
fg_t2d_ivw = mr_ivw(fg_t2d_obj)
ldl_t2d_ivw = mr_ivw(ldl_t2d_obj)
tc_t2d_ivw = mr_ivw(tc_t2d_obj)

bw_ast_obj = mr_input(bx = bw_ast$beta_bw, bxse = bw_ast$se_bw, by = bw_ast$beta_ast, byse = bw_ast$se_ast)
fg_ast_obj = mr_input(bx = fg_ast$beta_fg, bxse = fg_ast$se_fg, by = fg_ast$beta_ast, byse = fg_ast$se_ast)
ldl_ast_obj = mr_input(bx = ldl_ast$beta_ldl, bxse = ldl_ast$se_ldl, by = ldl_ast$beta_ast, byse = ldl_ast$se_ast)
tc_ast_obj = mr_input(bx = tc_ast$beta_tc, bxse = tc_ast$se_tc, by = tc_ast$beta_ast, byse = tc_ast$se_ast)

bw_ast_ivw = mr_ivw(bw_ast_obj)
fg_ast_ivw = mr_ivw(fg_ast_obj)
ldl_ast_ivw = mr_ivw(ldl_ast_obj)
tc_ast_ivw = mr_ivw(tc_ast_obj)

########################################################################################################################
#MR-cML-DP estimates
########################################################################################################################
set.seed(20230412)
bw_cad_cml = mr_cML(bw_cad_obj, n = 153781)
fg_cad_cml = mr_cML(fg_cad_obj, n = 184305)
ldl_cad_cml = mr_cML(ldl_cad_obj, n = 184305)
tc_cad_cml = mr_cML(tc_cad_obj, n = 184305)

bw_t2d_cml = mr_cML(bw_t2d_obj, n = 153781)
fg_t2d_cml = mr_cML(fg_t2d_obj, n = 281416)
ldl_t2d_cml = mr_cML(ldl_t2d_obj, n = 188577)
tc_t2d_cml = mr_cML(tc_t2d_obj, n = 188577)

bw_ast_cml = mr_cML(bw_ast_obj, n = 142486)
fg_ast_cml = mr_cML(fg_ast_obj, n = 142486)
ldl_ast_cml = mr_cML(ldl_ast_obj, n = 142486)
tc_ast_cml = mr_cML(tc_ast_obj, n = 142486)

########################################################################################################################
#MR-cML estimates
########################################################################################################################
set.seed(20230412)
bw_cad_cml0 = mr_cML(bw_cad_obj, n = 153781, DP = FALSE)
fg_cad_cml0 = mr_cML(fg_cad_obj, n = 184305, DP = FALSE)
ldl_cad_cml0 = mr_cML(ldl_cad_obj, n = 184305, DP = FALSE)
tc_cad_cml0 = mr_cML(tc_cad_obj, n = 184305, DP = FALSE)

bw_t2d_cml0 = mr_cML(bw_t2d_obj, n = 153781, DP = FALSE)
fg_t2d_cml0 = mr_cML(fg_t2d_obj, n = 281416, DP = FALSE)
ldl_t2d_cml0 = mr_cML(ldl_t2d_obj, n = 188577, DP = FALSE)
tc_t2d_cml0 = mr_cML(tc_t2d_obj, n = 188577, DP = FALSE)

bw_ast_cml0 = mr_cML(bw_ast_obj, n = 142486, DP = FALSE)
fg_ast_cml0 = mr_cML(fg_ast_obj, n = 142486, DP = FALSE)
ldl_ast_cml0 = mr_cML(ldl_ast_obj, n = 142486, DP = FALSE)
tc_ast_cml0 = mr_cML(tc_ast_obj, n = 142486, DP = FALSE)

########################################################################################################################
#MR-Horse estimates
########################################################################################################################
set.seed(20230412)
bw_cad_horse = mr_horse(bw_cad_obj, no_ini = 3, variable.names = "theta")
fg_cad_horse = mr_horse(fg_cad_obj, no_ini = 3, variable.names = "theta")
ldl_cad_horse = mr_horse(ldl_cad_obj, no_ini = 3, variable.names = "theta")
tc_cad_horse = mr_horse(tc_cad_obj, no_ini = 3, variable.names = "theta")

bw_t2d_horse = mr_horse(bw_t2d_obj, no_ini = 3, variable.names = "theta")
fg_t2d_horse = mr_horse(fg_t2d_obj, no_ini = 3, variable.names = "theta")
ldl_t2d_horse = mr_horse(ldl_t2d_obj, no_ini = 3, variable.names = "theta")
tc_t2d_horse = mr_horse(tc_t2d_obj, no_ini = 3, variable.names = "theta")

bw_ast_horse = mr_horse(bw_ast_obj, no_ini = 3, variable.names = "theta")
fg_ast_horse = mr_horse(fg_ast_obj, no_ini = 3, variable.names = "theta")
ldl_ast_horse = mr_horse(ldl_ast_obj, no_ini = 3, variable.names = "theta")
tc_ast_horse = mr_horse(tc_ast_obj, no_ini = 3, variable.names = "theta")

########################################################################################################################
#CAUSE
########################################################################################################################
bw_all = read.table(file = "bw_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                        "character", "numeric", "numeric", "numeric", "numeric"),
                    header = TRUE)
names(bw_all) = c("SNP", "chr", "pos", "EA", "NEA", "EAF", "beta_bw", "se_bw", "p_bw")
bw_sig = read.table('bw_cause_sig_out.txt', colClasses = c("character", "character", "numeric", "character", "character",
                                                           "numeric", "numeric", "numeric"), header = FALSE)
names(bw_sig) = c("SNP", "chr", "pos", "EA", "NEA", "EAF", "beta", "se", "p")

fg_all = read.table(file = "fg_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                            "character", "numeric", "numeric", "numeric", "numeric"),
                    header = TRUE)
names(fg_all) = c("SNP", "chr", "pos", "EA", "NEA", "EAF", "beta_fg", "se_fg", "p_fg")
fg_sig = read.table('fg_cause_sig_out.txt', colClasses = c("character", "character", "numeric", "character", "character",
                                                           "numeric", "numeric", "numeric"), header = FALSE)

names(fg_sig) = c("SNP", "chr", "pos", "EA", "NEA", "EAF", "beta", "se", "p")

ldl_all = read.table(file = "ldl0_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                              "character", "numeric", "numeric", "numeric", "numeric"),
                     header = TRUE)
names(ldl_all) = c("SNP", "chr", "pos", "NEA", "EA", "EAF", "beta_ldl", "se_ldl", "p_ldl")
ldl_sig = read.table('ldl0_cause_sig_out.txt', colClasses = c("character", "character", "numeric", "character", "character",
                                                           "numeric", "numeric", "numeric"), header = FALSE)
names(ldl_sig) = c("SNP", "chr", "pos", "NEA", "EA", "EAF", "beta", "se", "p")

tc_all = read.table(file = "tc0_out.txt", colClasses = c("character", "character", "numeric", "character",
                                                            "character", "numeric", "numeric", "numeric", "numeric"),
                    header = TRUE)
names(tc_all) = c("SNP", "chr", "pos", "NEA", "EA", "EAF", "beta_tc", "se_tc", "p_tc")
tc_sig = read.table('tc0_cause_sig_out.txt', colClasses = c("character", "character", "numeric", "character", "character",
                                                           "numeric", "numeric", "numeric"), header = FALSE)
names(tc_sig) = c("SNP", "chr", "pos", "NEA", "EA", "EAF", "beta", "se", "p")

#Estimates
bw_cad_cause = mr_cause_run(Xgwas = data.frame("SNP" = bw_all$SNP, "beta" = bw_all$beta_bw, "se" = bw_all$se_bw,
                                               "EA" = bw_all$EA, "NEA" = bw_all$NEA,
                                               "chr" = bw_all$chr, "pos" = bw_all$pos, "p" = bw_all$p_bw),
                            Ygwas = data.frame("SNP" = cad_dat$SNP, "beta" = cad_dat$beta_cad, "se" = cad_dat$se_cad,
                                               "EA" = cad_dat$EA_cad, "NEA" = cad_dat$NEA_cad),
                            sig_snps = bw_sig)

fg_cad_cause = mr_cause_run(Xgwas = data.frame("SNP" = fg_all$SNP, "beta" = fg_all$beta_fg, "se" = fg_all$se_fg,
                                               "EA" = fg_all$EA, "NEA" = fg_all$NEA,
                                               "chr" = fg_all$chr, "pos" = fg_all$pos, "p" = fg_all$p_fg),
                            Ygwas = data.frame("SNP" = cad_dat$SNP, "beta" = cad_dat$beta_cad, "se" = cad_dat$se_cad,
                                               "EA" = cad_dat$EA_cad, "NEA" = cad_dat$NEA_cad),
                            sig_snps = fg_sig)

ldl_cad_cause = mr_cause_run(Xgwas = data.frame("SNP" = ldl_all$SNP, "beta" = ldl_all$beta_ldl, "se" = ldl_all$se_ldl,
                                                "EA" = ldl_all$EA, "NEA" = ldl_all$NEA,
                                                "chr" = ldl_all$chr, "pos" = ldl_all$pos, "p" = ldl_all$p_ldl),
                             Ygwas = data.frame("SNP" = cad_dat$SNP, "beta" = cad_dat$beta_cad, "se" = cad_dat$se_cad,
                                                "EA" = cad_dat$EA_cad, "NEA" = cad_dat$NEA_cad),
                             sig_snps = ldl_sig)

tc_cad_cause = mr_cause_run(Xgwas = data.frame("SNP" = tc_all$SNP, "beta" = tc_all$beta_tc, "se" = tc_all$se_tc,
                                               "EA" = tc_all$EA, "NEA" = tc_all$NEA,
                                               "chr" = tc_all$chr, "pos" = tc_all$pos, "p" = tc_all$p_tc),
                            Ygwas = data.frame("SNP" = cad_dat$SNP, "beta" = cad_dat$beta_cad, "se" = cad_dat$se_cad,
                                               "EA" = cad_dat$EA_cad, "NEA" = cad_dat$NEA_cad),
                            sig_snps = tc_sig)

bw_t2d_cause = mr_cause_run(Xgwas = data.frame("SNP" = bw_all$SNP, "beta" = bw_all$beta_bw, "se" = bw_all$se_bw,
                                               "EA" = bw_all$EA, "NEA" = bw_all$NEA,
                                               "chr" = bw_all$chr, "pos" = bw_all$pos, "p" = bw_all$p_bw),
                            Ygwas = data.frame("SNP" = t2d_dat$SNP, "beta" = t2d_dat$beta_t2d, "se" = t2d_dat$se_t2d,
                                               "EA" = t2d_dat$EA_t2d, "NEA" = t2d_dat$NEA_t2d),
                            sig_snps = bw_sig)

fg_t2d_cause = mr_cause_run(Xgwas = data.frame("SNP" = fg_all$SNP, "beta" = fg_all$beta_fg, "se" = fg_all$se_fg,
                                               "EA" = fg_all$EA, "NEA" = fg_all$NEA,
                                               "chr" = fg_all$chr, "pos" = fg_all$pos, "p" = fg_all$p_fg),
                            Ygwas = data.frame("SNP" = t2d_dat$SNP, "beta" = t2d_dat$beta_t2d, "se" = t2d_dat$se_t2d,
                                               "EA" = t2d_dat$EA_t2d, "NEA" = t2d_dat$NEA_t2d),
                            sig_snps = fg_sig)

set.seed(20230413)
ldl_t2d_cause = mr_cause_run(Xgwas = data.frame("SNP" = ldl_all$SNP, "beta" = ldl_all$beta_ldl, "se" = ldl_all$se_ldl,
                                                "EA" = ldl_all$EA, "NEA" = ldl_all$NEA,
                                                "chr" = ldl_all$chr, "pos" = ldl_all$pos, "p" = ldl_all$p_ldl),
                             Ygwas = data.frame("SNP" = t2d_dat$SNP, "beta" = t2d_dat$beta_t2d, "se" = t2d_dat$se_t2d,
                                                "EA" = t2d_dat$EA_t2d, "NEA" = t2d_dat$NEA_t2d),
                             sig_snps = ldl_sig)

tc_t2d_cause = mr_cause_run(Xgwas = data.frame("SNP" = tc_all$SNP, "beta" = tc_all$beta_tc, "se" = tc_all$se_tc,
                                               "EA" = tc_all$EA, "NEA" = tc_all$NEA,
                                               "chr" = tc_all$chr, "pos" = tc_all$pos, "p" = tc_all$p_tc),
                            Ygwas = data.frame("SNP" = t2d_dat$SNP, "beta" = t2d_dat$beta_t2d, "se" = t2d_dat$se_t2d,
                                               "EA" = t2d_dat$EA_t2d, "NEA" = t2d_dat$NEA_t2d),
                            sig_snps = tc_sig)

bw_ast_cause = mr_cause_run(Xgwas = data.frame("SNP" = bw_all$SNP, "beta" = bw_all$beta_bw, "se" = bw_all$se_bw,
                                               "EA" = bw_all$EA, "NEA" = bw_all$NEA,
                                               "chr" = bw_all$chr, "pos" = bw_all$pos, "p" = bw_all$p_bw),
                            Ygwas = data.frame("SNP" = ast_dat$SNP, "beta" = ast_dat$beta_ast, "se" = ast_dat$se_ast,
                                               "EA" = ast_dat$EA_ast, "NEA" = ast_dat$NEA_ast),
                            sig_snps = bw_sig)

fg_ast_cause = mr_cause_run(Xgwas = data.frame("SNP" = fg_all$SNP, "beta" = fg_all$beta_fg, "se" = fg_all$se_fg,
                                               "EA" = fg_all$EA, "NEA" = fg_all$NEA,
                                               "chr" = fg_all$chr, "pos" = fg_all$pos, "p" = fg_all$p_fg),
                            Ygwas = data.frame("SNP" = ast_dat$SNP, "beta" = ast_dat$beta_ast, "se" = ast_dat$se_ast,
                                               "EA" = ast_dat$EA_ast, "NEA" = ast_dat$NEA_ast),
                            sig_snps = fg_sig)

ldl_ast_cause = mr_cause_run(Xgwas = data.frame("SNP" = ldl_all$SNP, "beta" = ldl_all$beta_ldl, "se" = ldl_all$se_ldl,
                                               "EA" = ldl_all$EA, "NEA" = ldl_all$NEA,
                                               "chr" = ldl_all$chr, "pos" = ldl_all$pos, "p" = ldl_all$p_ldl),
                            Ygwas = data.frame("SNP" = ast_dat$SNP, "beta" = ast_dat$beta_ast, "se" = ast_dat$se_ast,
                                               "EA" = ast_dat$EA_ast, "NEA" = ast_dat$NEA_ast),
                            sig_snps = ldl_sig)

tc_ast_cause = mr_cause_run(Xgwas = data.frame("SNP" = tc_all$SNP, "beta" = tc_all$beta_tc, "se" = tc_all$se_tc,
                                               "EA" = tc_all$EA, "NEA" = tc_all$NEA,
                                               "chr" = tc_all$chr, "pos" = tc_all$pos, "p" = tc_all$p_tc),
                            Ygwas = data.frame("SNP" = ast_dat$SNP, "beta" = ast_dat$beta_ast, "se" = ast_dat$se_ast,
                                               "EA" = ast_dat$EA_ast, "NEA" = ast_dat$NEA_ast),
                            sig_snps = tc_sig)

########################################################################################################################
#Second applied example
########################################################################################################################
ad_dat = read.csv('beta_se_all.csv')
mvmrob = mr_mvinput(bx = cbind(ad_dat$hill_beta, ad_dat$okbay_beta, ad_dat$ukbb_beta),
                    bxse = cbind(ad_dat$hill_se, ad_dat$okbay_se, ad_dat$ukbb_se),
                    by = ad_dat$lambert_beta, byse = ad_dat$lambert_se, exposure = c("Int", "Edu", "Inc"))
set.seed(20230816)
ad_ivw = mr_mvivw(mvmrob, nx = c(248482, 293723, 361194))
set.seed(20230816)
ad_median = mr_mvmedian(mvmrob)
set.seed(20230816)
ad_cml = mr_mvcML(mvmrob, n = 54162, DP = FALSE)
set.seed(20230816)
ad_cml_dp = mr_mvcML(mvmrob, n = 54162, DP = TRUE)
set.seed(20230816)
ad_horse = mvmr_horse(mvmrob, no_ini = 3, variable.names = "theta")
