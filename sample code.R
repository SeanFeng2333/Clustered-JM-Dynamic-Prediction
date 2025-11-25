library(ggplot2)
library(bayesplot)
library(patchwork)
library(cowplot)
library(simsurv)
library(dplyr)
library(tidyverse)
library(cmdstanr)
library(nlme)
library(JMbayes2)
library(grid)


path <- "/Users/f.x.y/Library/Mobile Documents/com~apple~CloudDocs/UofT/PhD Thesis/publication/"

sample_survdata <- read.csv(paste0(path, "sample_survdata.csv")) # data for the survival submodel 
sample_longdata <- read.csv(paste0(path, "sample_longdata.csv")) # data for the longitudianl submodel 
source(paste0(path, "dynamic prediction function.R")) # function for dynamic prediction 

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


# data for stan 
JM_datalist <- list(
  N=longdata_row, nteeth=nteeth, nid=nid, y1=y1, y2=y2, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, ID=ID, IDtooth=IDtooth, ID_s=ID_s,
  IDtooth_s=IDtooth_s, visits=longdata_obstime, times=eventtime, status=status
)

# stan model 
JM_2biom_app_mod = cmdstan_model("/Users/f.x.y/Library/Mobile Documents/com~apple~CloudDocs/UofT/PhD Thesis/publication/clustered_JM_2biom_correlated_application.stan")


start_t <- Sys.time()
fit_JM_cluster <- JM_2biom_app_mod$sample(data = JM_datalist,
                                         seed = 123456,
                                         chains = 3,
                                         parallel_chains = 3,
                                         iter_warmup = 1000,
                                         iter_sampling = 1000
)

end_t <- Sys.time()
stan_t <- end_t - start_t

# diagnostic of the HMC sampling
fit_JM_cluster$diagnostic_summary()


draws_df_cluster_all <- fit_JM_cluster$draws(format = "df")
draws_df_cluster <- draws_df_cluster_all[c("beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[5]","beta1[6]","beta1[7]",
                                           "beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[5]","beta2[6]","beta2[7]",
                                           "gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]","alpha1","alpha2",
                                           "sigma_e",
                                           "sigma_bi[1]","sigma_bi[2]",
                                           "sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]",
                                           "sigma_ri",
                                           "R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]",
                                           "R_bi[1,2]")]

# summary of parameters 
summary_2biom <- posterior::summarize_draws(
  draws_df_cluster,
  posterior::default_summary_measures()[1:4],
  extra_quantiles = ~posterior::quantile2(., probs = c(.025, .975)),
  posterior::default_convergence_measures())


print(summary_2biom, n=37)


mcmc_hist(draws_df_cluster)

# trace plots 
traceplot1 <- bayesplot::mcmc_trace(draws_df_cluster, pars=c("beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[5]","beta1[6]","beta1[7]",
                                                                      "beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[5]","beta2[6]","beta2[7]",
                                                                      "gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]"))


traceplot2 <- bayesplot::mcmc_trace(draws_df_cluster, pars=c("alpha1","alpha2", "sigma_e",
                                                                      "sigma_bi[1]","sigma_bi[2]",                                                 
                                                                      "sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]",
                                                                      "sigma_ri",
                                                                      "R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]",
                                                                      "R_bi[1,2]"))



#### dynamic prediction ####

# landmark time 
predtime <- 8


# prediction time window 
window <- seq(0, max(sample_survdata$eventtime)-predtime, 0.5)
fwindow <- as.vector(predtime + window)


# select a subject and a tooth to make prediction 
new_i <- 2
new_ij <- 213

# data of the sampled subject
sample_survdata <- sample_survdata %>% 
  mutate(i_index = match(id, sort(unique(id))),
         ij_index = match(id_tooth, sort(unique(id_tooth))))


new_ij_index <- sample_survdata[sample_survdata$id_tooth==new_ij,]["ij_index"][[1]]
new_i_index <- sample_survdata[sample_survdata$id_tooth==new_ij,]["i_index"][[1]]


object <- draws_df_cluster_all

newdata <- sample_longdata[sample_longdata$id_tooth == new_ij, ]


# sampling output of the selected subject 
new_ij_draws_df <- draws_df_cluster_all[c("beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[5]","beta1[6]","beta1[7]",
                                                   "beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[5]","beta2[6]","beta2[7]",
                                                   "gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]","alpha1","alpha2",
                                                   "sigma_e",
                                                   "sigma_bi[1]","sigma_bi[2]",
                                                   "sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]",
                                                   "sigma_ri",
                                                   "R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]",
                                                   "R_bi[1,2]", 
                                                   paste0("bi[",new_i_index,",1]"),paste0("bi[",new_i_index,",2]"),
                                                   paste0("bij[",new_ij_index,",1]"),paste0("bij[",new_ij_index,",2]"), 
                                                   paste0("bij[",new_ij_index,",3]"),paste0("bij[",new_ij_index,",4]"))]


summary_new_ij_draws <- posterior::summarize_draws(
  new_ij_draws_df,
  posterior::default_summary_measures()[1:4],
  extra_quantiles = ~posterior::quantile2(., probs = c(.025, .975)),
  posterior::default_convergence_measures())
print(summary_new_ij_draws, n=43)

new_ij_draws.est <- summary_new_ij_draws$mean
new_ij_draws.ci_lower <- summary_new_ij_draws$q2.5
new_ij_draws.ci_upper <- summary_new_ij_draws$q97.5



X_long <- c(newdata$sex, newdata$std_age, newdata$smoke, newdata$molar, newdata$lower)
previous.time <- newdata$followup_year  # observed times before prediction 


fitted_y1.est <- matrix(c(rep(1, length(previous.time)), 
                          X_long, previous.time), ncol = 7, byrow = F) %*% new_ij_draws.est[c(1:7)] + 
  matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
           previous.time), ncol = 3, byrow = F) %*% c(new_ij_draws.est[38], new_ij_draws.est[40:41])

fitted_y1.ci_lower <- matrix(c(rep(1, length(previous.time)), 
                               X_long, previous.time), ncol = 7, byrow = F) %*% new_ij_draws.ci_lower[c(1:7)] + 
  matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
           previous.time), ncol = 3, byrow = F) %*% c(new_ij_draws.ci_lower[38], new_ij_draws.ci_lower[40:41])

fitted_y1.ci_upper <- matrix(c(rep(1, length(previous.time)), 
                               X_long, previous.time), ncol = 7, byrow = F) %*% new_ij_draws.ci_upper[c(1:7)] + 
  matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
           previous.time), ncol = 3, byrow = F) %*% c(new_ij_draws.ci_upper[38], new_ij_draws.ci_upper[40:41])


fitted_y2.mu <- matrix(c(rep(1, length(previous.time)), 
                         X_long, previous.time), ncol = 7, byrow = F) %*% new_ij_draws.est[c(8:14)] + 
  matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
           previous.time), ncol = 3, byrow = F) %*% c(new_ij_draws.est[39], new_ij_draws.est[42:43])

fitted_y2.mu.ci_lower <- matrix(c(rep(1, length(previous.time)), 
                                  X_long, previous.time), ncol = 7, byrow = F) %*% new_ij_draws.ci_lower[c(8:14)] + 
  matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
           previous.time), ncol = 3, byrow = F) %*% c(new_ij_draws.ci_lower[39], new_ij_draws.ci_lower[42:43])

fitted_y2.mu.ci_upper <- matrix(c(rep(1, length(previous.time)), 
                                  X_long, previous.time), ncol = 7, byrow = F) %*% new_ij_draws.ci_upper[c(8:14)] + 
  matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
           previous.time), ncol = 3, byrow = F) %*% c(new_ij_draws.ci_upper[39], new_ij_draws.ci_upper[42:43])

fitted_y2 <- exp(fitted_y2.mu)/(1+exp(fitted_y2.mu))
fitted_y2.ci_lower <- exp(fitted_y2.mu.ci_lower)/(1+exp(fitted_y2.mu.ci_lower))
fitted_y2.ci_upper <- exp(fitted_y2.mu.ci_upper)/(1+exp(fitted_y2.mu.ci_upper))


newdata_fitted <- data.frame(id = newdata$id,
                             tooth = newdata$tooth,
                             idtooth = newdata$id_tooth,
                             PPD = newdata$max_ppd,
                             logPPD = log(newdata$max_ppd),
                             mobility = newdata$mobility,
                             visit = newdata$visit,
                             year = newdata$followup_year,
                             status = newdata$status, 
                             eventtime = newdata$eventtime,
                             fitted_y1.est, fitted_y1.ci_lower, fitted_y1.ci_upper,
                             fitted_y2, fitted_y2.ci_lower, fitted_y2.ci_upper,
                             sex = newdata$sex, std_age = newdata$std_age, smoke = newdata$smoke, 
                             molar = newdata$molar, lower = newdata$lower
)


clusteredJM_unit_DP <- unitlevel_DP_JM(object=draws_df_cluster_all, 
                                       object_summary=summary_2biom,
                                       surv_data=sample_survdata, 
                                       long_data=sample_longdata, 
                                       type = c("SurvProb", "Density"), 
                                       pred_ij = new_ij, 
                                       idVar = "id", 
                                       simulate = TRUE, 
                                       survTimes = fwindow, 
                                       pred.time = predtime, 
                                       LeftTrunc_var = NULL, 
                                       MCMCnum = 200L, 
                                       CI.levels = c(0.025, 0.975), 
                                       log = FALSE, 
                                       scale = 1, 
                                       init.b = NULL, 
                                       seed = 1L)

predsurv_df <- data.frame(clusteredJM_unit_DP)


#### dynamic prediction plot ####

x_limits <- c(0, tail(fwindow,1)); x_breaks <- seq(0, tail(fwindow,1), by = 2)

p_top_left <- ggplot() +
  geom_point(data = newdata_fitted, aes(year, logPPD), alpha = 0.7) +
  geom_line(data = newdata_fitted, aes(year, y = fitted_y1.est), linewidth = 1, , color = "blue") +
  geom_vline(xintercept = predtime, linetype = 2) +
  scale_x_continuous(limits = c(0, fwindow[1]), breaks = seq(0, fwindow[1], by = 2), expand = c(0, 0)) +
  labs(x = NULL, y = "log(PPD)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x  = element_blank(),   
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),                 
        axis.line.x  = element_line(color = "black"), 
        axis.line.y  = element_line(color = "black"),
        axis.title.y = element_text(size = 12, margin = margin(r = 6)),
        #        panel.spacing = unit(0, "pt"),
        plot.margin = margin(10,2.5,5,10))    # tighten white space around panel



p_bottom_left <- ggplot() +
  geom_point(data = newdata_fitted, aes(year, mobility),
             alpha = 0.6) +
  geom_smooth(data = newdata_fitted, aes(year, y = fitted_y2), linewidth = 1, 
              se = FALSE, color = "blue") +
  geom_vline(xintercept = predtime, linetype = 2) +
  scale_x_continuous(limits = c(0, fwindow[1]), breaks = seq(0, fwindow[1], by = 2), expand = c(0, 0)) +
  scale_y_continuous("Mobility", limits = c(0, 1.1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = NULL)  +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x  = element_text(),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        axis.line.x  = element_line(color = "black"),
        axis.line.y  = element_line(color = "black"),
        axis.title.y = element_text(size = 12, margin = margin(r = 6)),
        plot.margin = margin(10,2.5,5,10))  


p_right <- ggplot(predsurv_df, aes(times, Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2) +
  geom_line(linewidth = 1, color = "red") +
  scale_x_continuous(limits = c(fwindow[1], 16), breaks = seq(10, 16, by = 2), expand = c(0, 0)) +
  scale_y_continuous("Survival probability", limits = c(0, 1), breaks = seq(0, 1, by = 0.1),position = "right") +
  labs(x = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.x  = element_text(),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        axis.line.x  = element_line(color = "black"),
        axis.line.y  = element_line(color = "black"),
        panel.spacing = unit(0, "pt"),
        plot.margin = margin(10,10,5,0))



left_col  <- plot_grid(p_top_left, p_bottom_left, ncol = 1, align = "v")
fig_body <- plot_grid(left_col, p_right, ncol = 2, rel_widths = c(fwindow[1]+2.5, 16-fwindow[1]+2.5))

xlab_grob <- ggdraw() +
  draw_label("Year", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5, size = 12)

final_fig <- plot_grid(fig_body, xlab_grob, ncol = 1, rel_heights = c(1, 0.07))
final_fig





