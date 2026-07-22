library(ggplot2)
library(bayesplot)
library(patchwork)
library(cowplot)
library(simsurv)
library(dplyr)
library(tidyverse)
install.packages("Rtools44")
library(Rtools44)
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(cmdstanr)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_version() # verify the version 
cmdstanr::check_cmdstan_toolchain() # verify setup

library(nlme)
library(JMbayes2)
library(grid)


sample_survdata <- read.csv(paste0("sample_survdata.csv")) # data for the survival submodel 
sample_longdata <- read.csv(paste0("sample_longdata.csv")) # data for the longitudianl submodel 
source(paste0("dynamic prediction function.R")) # function for dynamic prediction 

# stan model 
path <- "/Users/f.x.y/Library/Mobile Documents/com~apple~CloudDocs/UofT/PhD Thesis/publication/"
JM_2biom_app_mod = cmdstan_model(paste0(path, "clustered_JM_2biom_correlated_application_lkj_weibull.stan"))


#### Model fitting #### 

# Data for the Stan model
X1 <- sample_survdata$sex # covariates
X2 <- sample_survdata$std_age
X3 <- sample_survdata$smoke
X4 <- sample_survdata$molar
X5 <- sample_survdata$lower
y1 <- log(sample_longdata$max_ppd)             # longitudinal outcomes
y2 <- sample_longdata$mobility          # longitudinal outcomes
ID <- as.numeric(factor(sample_longdata$id, levels = unique(sample_longdata$id)))  # patient IDs
IDtooth <- as.numeric(factor(sample_longdata$id_tooth, levels = unique(sample_longdata$id_tooth))) # tooth IDs

ID_s <- as.numeric(factor(sample_survdata$id, levels = unique(sample_survdata$id))) 
IDtooth_s <- unique(IDtooth)

nid <- length(unique(ID))        # number of patients
nteeth <- length(X1)                     # total number of teeth

status <- sample_survdata$status            # status (1 = dead, 0 = alive)
eventtime <- sample_survdata$eventtime              # times to event in year
longdata_obstime <- sample_longdata$followup_year         # visit times for repeated observations
longdata_row <- length(y1)                   # total number of longitudinal outcomes


sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
        -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
        0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
        0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
        0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
length_GK <- length(sk)


JM_datalist <- list(
  N=longdata_row, nteeth=nteeth, nid=nid, y1=y1, y2=y2, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5,  
  ID=ID, IDtooth=IDtooth, ID_s=ID_s,
  IDtooth_s=IDtooth_s, visits=longdata_obstime, times=times, status=status, 
  sk=sk, wk=wk, length_GK=length_GK
)


start_t <- Sys.time()
fit_JM_cluster <- JM_2biom_app_mod$sample(data = JM_datalist,
                                         seed = 123456,
                                         chains = 4,
                                         parallel_chains = 4,
                                         iter_warmup = 1000,
                                         iter_sampling = 1000
)

end_t <- Sys.time()
stan_t <- end_t - start_t


clustered_jm.fit <- fit_JM_cluster$summary(
  variables = c("gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]","alpha1","alpha2",
                "beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[7]","beta1[5]","beta1[6]","sigma_e",
                "beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[7]","beta2[5]","beta2[6]",
                "sigma_ri", "xi",
                "sigma_bi[1]","sigma_bi[2]", "R_bi[1,2]",
                "sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]",
                "R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]"),
  posterior::default_summary_measures()[1:4],
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975)),
  posterior::default_convergence_measures()
)


