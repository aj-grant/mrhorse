# mr-horse
Provides the R code to implement the methods described in

Grant, AJ and Burgess, S (2023). [A Bayesian approach to Mendelian randomization using summary statistics in the univariable and multivariable settings with correlated pleiotropy](https://www.biorxiv.org/content/10.1101/2023.05.30.542988v1). *bioRxiv*, doi: https://doi.org/10.1101/2023.05.30.542988

The R code to reproduce the simulation results and applied examples are also provided.

The paper introduces a method for performing both univariable and multivariable Mendelian randomization to examine the
causal effect of one or more exposures on an outcome using genetic association data. The methods take as input estimates
of the associations between genetic instruments and the exposure(s) and outcome from GWAS summary statistics. Causal
effect estimation is performed in a Bayesian framework using a horseshoe shrinkage prior to account for both correlated
and uncorrelated pleiotropic effects.

The scripts in the [Simulation scripts](https://github.com/aj-grant/mrhorse/tree/main/Simulation%20scripts) and
[Applied example](https://github.com/aj-grant/mrhorse/tree/main/Applied%20example) folders are self-contained, however
the MR-Horse and MVMR-Horse methods may also be implemented using the code provided in mr_horse.R.

## Installation
In order to run the code, [JAGS](https://sourceforge.net/projects/mcmc-jags/) first needs to be installed.

The functions use the R2jags package.

```
install.packages("R2jags")
library(R2jags)
```

The file mr_horse.R can either be copied and sourced locally, or sourced via:
```
library(devtools)
source_URL("https://github.com/aj-grant/mr_horse.R")
```

## Implementation
The function `mr_horse()` implements the univariable MR-Horse method. The required input is a data frame with column headings:

- `betaY`: estimate of genetic association with the outcome
- `betaYse`: standard error of the estimate of genetic association with the outcome
- `betaX`: estimate of genetic association with the exposure
- `betaXse`: standard error of the estimate of genetic association with the exposure

Each row in the data frame represents a genetic variant.

Optional arguments are:

- `no_ini`: number of chains to run (default = 3)
- `variable.names`: vector of parameters to save in the MCMC output. The causal effect is "theta", and will always be saved.
Other relevant parameters include "alpha" (pleiotropic effects for each genetic variant) and "rho" (correlations between
pleiotropic effects and genetic variant-exposure effects for each variant)
- `n.iter`: number of iterations (in addition to burn-in, default = 10000)
- `n.burnin`: number of iterations for burn-in (default = 10000)

Output from the `mr_horse()` function include:

- `$MR_Estimate`: a data frame with the causal effect estimate (which is the posterior mean), standard deviation (i.e.,
the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value
- `$MR_Coda`: full MCMC samples for all parameters in `variable.names`

JAGS plotting tools can be used with the MCMC output, for example using `traceplot()` and `densplot()`.

The function `mvmr_horse` implements the multivariable MVMR-Horse method. The required inputs are as above, but
where `betaX` is replaced by the K columns `betaX1`, `betaX2`, ..., and `betaXse` is replaced by the K columns
`betaX1se`, `betaX2se`, ... .

## Examples
### Univariable MR
The csv file dat_ex.csv contains a dataframe containing simulated genetic association estimates between 100 genetic
instruments and an exposure and outcome, as well as their corresponding standard errors. This is taken from the first
replication of the first simulation study (that is, where 20% of variants are pleiotropic, and pleiotropy is balanced).
The dataset can be analysed as follows.

```
data_ex = read.csv(file = "data_ex.csv", header = TRUE)
set.seed(20230531)
MREx = mr_horse(data_ex)
MREx$MR_Estimate
  Estimate    SD 2.5% quantile 97.5% quantile Rhat
1    0.097 0.018         0.063          0.132    1
traceplot(MREx$MR_Coda[, "theta"])
densplot(MREx$MR_Coda[, "theta"])
```

### Multivariable MR
The csv file dat_mv_ex.csv contains a dataframe containing simulated genetic association estimates between 100 genetic
instruments and two exposures and an outcome, as well as their corresponding standard errors. This is taken from the first
replication of the multivariable simulation study (that is, where 20% of variants are pleiotropic, and pleiotropy is balanced).
The dataset can be analysed as follows.

```
data_mv_ex = read.csv(file = "data_mv_ex.csv", header = TRUE)
set.seed(20230531)
MVMREx = mvmr_horse(data_mv_ex)
MVMREx$MR_Estimate
  Parameter Estimate    SD 2.5% quantile 97.5% quantile  Rhat
1  theta[1]    0.100 0.019         0.065          0.137 1.000
2  theta[2]    0.105 0.018         0.069          0.139 1.001
traceplot(MVMREx$MR_Coda[, "theta[1]"])
densplot(MVMREx$MR_Coda[, "theta[1]"])
```
