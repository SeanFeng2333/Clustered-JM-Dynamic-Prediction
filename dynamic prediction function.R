

predictive_metrics_JM <- function(object=clustered_jm.est, 
                                  surv_data=data_surv_full_longi, 
                                  long_data=data_long_full_longi, 
                                  type = c("SurvProb", "Density"), 
                                  idVar = "id", 
                                  simulate = TRUE, 
                                  survTimes = fwindow, 
                                  pred.time = predtime, 
                                  LeftTrunc_var = NULL, 
                                  MCMCnum = 200L, 
                                  CI.levels = c(0.025, 0.975), 
                                  log = FALSE, 
                                  scale = 1.6, 
                                  init.b = NULL, 
                                  seed = 1L) {
  
  
  TT <- sim_data_long$obstime 
  if (is.null(survTimes) || !is.numeric(survTimes)) {
    survTimes <- seq(min(TT), quantile(TT, 0.90) + 0.01, length.out = 35L)
  }
  
  #### prepare data ####
  timeVar <- "obstime"
  hasScale <- T # biomarker is continuous, has measurement error
  gammas <- clustered_jm.est[c("gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]")]
  betas1 <- clustered_jm.est[c("beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[5]","beta1[6]","beta1[7]")]
  betas2 <- clustered_jm.est[c("beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[5]","beta2[6]","beta2[7]")]
  xi <- clustered_jm.est[c("xi")]
  alpha1 <- clustered_jm.est[c("alpha1")]
  alpha2 <- clustered_jm.est[c("alpha2")]
  sigma_e <- clustered_jm.est["sigma_e"]
  sigma_bi <-clustered_jm.est[c("sigma_bi[1]","sigma_bi[2]")]
  sigma_bij <- clustered_jm.est[c("sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]")]
  sigma_ri <- clustered_jm.est["sigma_ri"]
  R_bi <- clustered_jm.est["R_bi[1,2]"]
  R_bij <- clustered_jm.est[c("R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]")]
  L_bi <- clustered_jm.est[c("L_bi[1,1]", "L_bi[1,2]", "L_bi[2,1]","L_bi[2,2]")]
  L_bij <- clustered_jm.est[c("L_bij[1,1]", "L_bij[1,2]", "L_bij[1,3]","L_bij[1,4]",
                              "L_bij[2,1]", "L_bij[2,2]", "L_bij[2,3]","L_bij[2,4]",
                              "L_bij[3,1]", "L_bij[3,2]", "L_bij[3,3]","L_bij[3,4]",
                              "L_bij[4,1]", "L_bij[4,2]", "L_bij[4,3]","L_bij[4,4]")]
  
  var_ri = sigma_ri^2
  var_e = sigma_e^2
  
  length_bi <- length(sigma_bi)
  Cor_bi <- diag(1, length_bi)
  Cor_bi[lower.tri(Cor_bi)] <- R_bi
  Cor_bi <- Cor_bi + t(Cor_bi) - diag(diag(Cor_bi))  # Make symmetric
  var_bi <- diag(sigma_bi) %*% Cor_bi %*%  diag(sigma_bi)
  
  length_bij <- length(sigma_bij)
  Cor_bij <- diag(1, length_bij)
  Cor_bij[lower.tri(Cor_bij)] <- R_bij
  Cor_bij <- Cor_bij + t(Cor_bij) - diag(diag(Cor_bij))  # Make symmetric
  var_bij <- diag(sigma_bij) %*% Cor_bij %*%  diag(sigma_bij)
  L_bi_matrix <- t(chol(var_bi))
  L_bij_matrix <- t(chol(var_bij))
  
  max.time <- max(TT)
  
  # create an indicator for censoring before the landmark time
  # data have observation before the landmark time
  long_data <- long_data %>% filter(followup_year <= pred.time) %>% 
    group_by(id_tooth) %>%
    mutate(total_visits = n()) 
  
  id_have_obs <- unique(long_data$id_tooth)
  
  surv_data <- surv_data %>% 
    filter(id_tooth %in% id_have_obs) %>%
    mutate(censor_ind = ifelse((eventtime_year <= pred.time & status == 0), 1, 0))
  
  ## filter ij's who are at risk and not being censored at pred.time
  pred_data_surv <- surv_data %>% filter(eventtime_year > pred.time)
  pred_data_long <- long_data %>% 
    filter((id_tooth %in% pred_data_surv$id_tooth) & (followup_year <= pred.time)) %>% group_by(id_tooth) 
  
  ## filter ij's who are censored between pred.time and the prediction window
  pred_data_surv_censored <- surv_data %>% filter(eventtime_year <= survTimes[length(survTimes)] & status == 0)
  pred_data_long_censored <- long_data %>% 
    filter((id_tooth %in% pred_data_surv_censored$id_tooth)) %>% group_by(id_tooth) 
  
  n <- length(TT) # number of total observation time
  n.tp <- length(pred.time)
  ncz <- 7
  

  
  
  #### useful functions ####
  
  rmvt <- function (n, mu, Sigma, df) {
    p <- length(mu)
    if (is.list(Sigma)) {
      ev <- Sigma$values
      evec <- Sigma$vectors
    } else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
    }
    X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                               matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
    if (n == 1L) drop(X) else t.default(X)
  }
  
  dmvt <- function (x, mu, Sigma = NULL, invSigma = NULL, df, log = FALSE, prop = TRUE) {
    if (!is.numeric(x)) 
      stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
      x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
      }
      #      if (!all(ev >= -1e-06 * abs(ev[1]))) 
      #        stop("'Sigma' is not positive definite")
      invSigma <- evec %*% (t(evec)/ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
      fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + 
                                                                log(df)) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
      if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else as.vector(((1 + quad)^(-(df + p)/2)))
    }
  }
  
  densLong_norm <- function (y, eta.y, scale, log = FALSE, data) {
    dnorm(x = y, mean = eta.y, sd = scale, log = log)
  }
  
  dmvnorm <- function (x, mu, Sigma = NULL, invSigma = NULL, log = FALSE, prop = TRUE) {
    if (!is.matrix(x))
      x <- rbind(x)
    if (is.matrix(mu)) {
      if (nrow(mu) != nrow(x))
        stop("incorrect dimensions for 'mu'.")
      p <- ncol(mu)
    } else {
      p <- length(mu)
      mu <- rep(mu, each = nrow(x))
    }
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors            
      }
      invSigma <- evec %*% (t(evec) / ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- - determinant(as.matrix(invSigma))$modulus
    }
    ss <- x - mu
    quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
    if (!prop)
      fact <- - 0.5 * (p * log(2 * pi) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - quad) else as.vector(- quad)
    } else {
      if (!prop) as.vector(exp(fact - quad)) else as.vector(exp(- quad))
    }
  }
  
  densLong_binomial <- function (y, eta.y, scale, log = FALSE, data) {
    y*eta.y - log(1+exp(eta.y))
  }
  
  
  
  S_u_t <- function(t1, t0, X_surv, betas1, betas2, gammas, 
                    xi, alpha1, alpha2,
                    b_Ai,b_Bi,b_Aij0,b_Aij1,b_Bij0,b_Bij1,b_ri){
    sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
            -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
            0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
    wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
            0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
            0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
    length_GK <- length(sk)
    
    log_S_t <- 0
    for (j in 1:length_GK) {
      sk_scaled_j = (sk[j]+1)*(t1-t0)/2+t0
      wk_scaled_j = wk[j]*(t1-t0)/2
      
      log_S_t <- log_S_t - wk_scaled_j*xi*sk_scaled_j^(xi-1)*exp((
        alpha1*betas1[7]+alpha1*b_Aij1+alpha2*betas2[7]+alpha2*b_Bij1)*sk_scaled_j)*
        exp(b_ri + c(1,X_surv) %*% gammas + alpha1*(c(1, X_surv) %*% betas1[1:6] + b_Ai+ b_Aij0)
            + alpha2*(c(1, X_surv) %*% betas2[1:6]+b_Bi+b_Bij0))
    }
    exp(log_S_t)
    
  }
  
    #### Dynamic prediction of all units ####
  
  pd_u_t_clustered <- data.frame()
  pd_u_t_clustered_censored <- data.frame()
  
  t0 <- Sys.time()
  for (i in 1:length(unique(surv_data$researchID))) {
    id_i <- unique(surv_data$researchID)[i]
    data_long_i <- long_data %>% filter((researchID==id_i) )
    data_surv_i <- surv_data %>% filter(researchID==id_i )

    id_tooth_i <- data_surv_i$id_tooth
    
    ind_censor_i <- data_surv_i[["censor_ind"]]
    yA_ij <- log(data_long_i$max_ppd)
    yB_ij <- data_long_i$mobility 
    X1_ij <- data_surv_i$sex
    X2_ij <- data_surv_i$std_age
    X3_ij <- data_surv_i$smoke
    X4_ij <- data_surv_i$molar
    X5_ij <- data_surv_i$lower
    
    IDtooth <- as.numeric(factor(data_long_i$id_tooth, levels = unique(data_long_i$id_tooth)))
    
    IDtooth_s <- unique(IDtooth)
    nteeth <- length(X1_ij)                     # total number of teeth
    
    previous.time_ij <- data_long_i$followup_year
    
    obstimes_ij <- data_surv_i$eventtime_year
    obstimes_ij <- ifelse(obstimes_ij <= pred.time, obstimes_ij, pred.time)
    obsstatus_ij <- data_surv_i$status
    obsstatus_ij <- ifelse(data_surv_i$eventtime_year <= pred.time, obsstatus_ij, 0)
    
    status <- obsstatus_ij            # vital status (1 = dead, 0 = alive)
    times <- obstimes_ij           # times to event
    longdata_obstime <- previous.time_ij       # visit times for repeated observations
    longdata_row <- length(yA_ij)                   # total number of longitudinal outcomes
    
    sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
            -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
            0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
    wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
            0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
            0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
    length_GK <- length(sk)
    
    DP_datalist <- list(
      N=longdata_row, nteeth=nteeth, y1=yA_ij, y2=yB_ij, X1=X1_ij, X2=X2_ij, X3=X3_ij, X4=X4_ij, X5=X5_ij, 
      visits=longdata_obstime, times=times, status=status, 
      IDtooth=IDtooth, 
      IDtooth_s=IDtooth_s, 
      sk=sk, wk=wk, length_GK=length_GK, 
      beta1=betas1, beta2=betas2, gamma=gammas, xi=xi, alpha1=alpha1, alpha2=alpha2, 
      sigma_e=sigma_e, sigma_ri=sigma_ri, 
      sigma_bi=sigma_bi, 
      sigma_bij=sigma_bij,
      L_bi=L_bi_matrix, 
      L_bij=L_bij_matrix)
    
    
    fit_DP_cluster <- analysis_DP_mod$sample(data = DP_datalist,
                                             seed = 1234,
                                             chains = 4,
                                             parallel_chains = 4,
                                             iter_warmup = 500,
                                             iter_sampling = 500)
    
    draws_DP_cluster0 <- fit_DP_cluster$draws(format = "df")
    
    
    b_Ai.new <- mean(c(draws_DP_cluster0["bi[1]"])[[1]])
    b_Bi.new <- mean(c(draws_DP_cluster0["bi[2]"])[[1]])
    b_Aij0.new <- colMeans(draws_DP_cluster0[,4:(3+nteeth)])
    b_Aij1.new <- colMeans(draws_DP_cluster0[,(3+nteeth+1):(3+2*nteeth)])
    b_Bij0.new <- colMeans(draws_DP_cluster0[,(3+2*nteeth+1):(3+3*nteeth)])
    b_Bij1.new <- colMeans(draws_DP_cluster0[,(3+3*nteeth+1):(3+4*nteeth)])
    b_ri.new <- mean(c(draws_DP_cluster0["ri"])[[1]])
    
    for (ij in 1:nteeth) {
      id_tooth_ij = id_tooth_i[ij]
      # for non-censored units before landmark time
      if (ind_censor_i[ij] == 0){
        SS_ij <- S_u_t(survTimes,pred.time,c(X1_ij[ij], X2_ij[ij], X3_ij[ij], X4_ij[ij], X5_ij[ij]), betas1, betas2, gammas, xi, 
                       alpha1, alpha2, 
                       b_Ai.new, b_Bi.new, 
                       b_Aij0.new[ij], b_Aij1.new[ij], 
                       b_Bij0.new[ij], b_Bij1.new[ij], b_ri.new)
        
        SS_ij <- cbind(id_tooth_ij, SS_ij)
        
        pd_u_t_clustered <- rbind(pd_u_t_clustered, SS_ij)
      }
      
      # for units that are censored between landmark time and prediction window
      eventtime_ij <- data_surv_i$eventtime_year[ij]
      status_ij <- data_surv_i$status[ij]
      
      SS_ij <-  (eventtime_ij>=pred.time) * (eventtime_ij<=survTimes) * (status_ij==0) * 
      S_u_t(survTimes,eventtime_ij,c(X1_ij[ij], X2_ij[ij], X3_ij[ij], X4_ij[ij], X5_ij[ij]), betas1, betas2, gammas, xi, 
              alpha1, alpha2, 
              b_Ai.new, b_Bi.new, 
              b_Aij0.new[ij], b_Aij1.new[ij], 
              b_Bij0.new[ij], b_Bij1.new[ij], b_ri.new)
      
      SS_ij <- cbind(id_tooth_ij, SS_ij)
      
      pd_u_t_clustered_censored <- rbind(pd_u_t_clustered_censored, SS_ij)
    }
  }
  t1 <- Sys.time()
  DP_runtime <- t1-t0
  
  
  pd_u_t_clustered_wide <- pd_u_t_clustered %>% 
    group_by(id_tooth_ij) %>%
    mutate(pred_id = row_number()) %>%   # create index 1–6 within each id
    pivot_wider(
      names_from = pred_id,
      values_from = SS_ij,
      names_prefix = "pred_"
    ) %>%
    ungroup()
  colnames(pd_u_t_clustered_wide) <- c("id_tooth", paste0("year", survTimes))
  
  
  pd_u_t_clustered_censored_wide <- pd_u_t_clustered_censored %>% 
    group_by(id_tooth_ij) %>%
    mutate(pred_id = row_number()) %>%   # create index 1–6 within each id
    pivot_wider(
      names_from = pred_id,
      values_from = SS_ij,
      names_prefix = "pred_"
    ) %>%
    ungroup()
  colnames(pd_u_t_clustered_censored_wide) <- c("id_tooth", paste0("year", survTimes))
  
  
  
  data_surv_noncensor <- surv_data %>% 
    dplyr::select(researchID, tooth, id_tooth, eventtime_year, status, censor_ind)
  data_surv_censored <- surv_data %>%
    dplyr::select(researchID, tooth, id_tooth, eventtime_year, status, censor_ind)
  
  
  pd_u_t_clustered_df <- merge(data_surv_noncensor, pd_u_t_clustered_wide, by="id_tooth")
  pd_u_t_clustered_df <- pd_u_t_clustered_df %>% filter(eventtime_year >= pred.time)
  
  pd_u_t_clustered_censored_df <- merge(data_surv_censored, pd_u_t_clustered_censored_wide, by="id_tooth")
  pd_u_t_clustered_censored_df <- pd_u_t_clustered_censored_df %>% filter(id_tooth %in% pd_u_t_clustered_df$id_tooth)
  
  
  #### Brier Score (pec, IPCW) ####
  
  BrierScore_clustered <- pec(as.matrix(pd_u_t_clustered_df[,c(7:(6+length(survTimes)))]),
                              formula = Surv(eventtime_year, status) ~ 1,
                              data = pred_data_surv,
                              cens.model = "marginal", exact = F,
                              times = survTimes[2:length(survTimes)], start = pred.time, reference = FALSE)
  
  pec_bs <- BrierScore_clustered$AppErr[[1]]
  

  
  #### Brier Score (model-based) ####
  pe <- c()
  num_risk_t <- nrow(pred_data_surv)
  for (i in 1:length(survTimes)){
    u_i <- survTimes[i]
    pred_data_surv$alive_ind_i <- ifelse(pred_data_surv$eventtime_year > u_i, 1, 0)
    pred_data_surv$event_ind_i <- 1 - pred_data_surv$alive_ind_i
    pred_data_surv$ind1_i <- pred_data_surv$alive_ind_i
    pred_data_surv$ind2_i <- ifelse((pred_data_surv$event_ind_i==1)&(pred_data_surv$status==1), 1, 0)
    pred_data_surv$ind3_i <- ifelse((pred_data_surv$event_ind_i==1)&(pred_data_surv$status==0), 1, 0)
    
    pd_u_t_clustered_i <- pd_u_t_clustered_df[,i+6]
    pd_u_t_clustered_censored_df_i <- pd_u_t_clustered_censored_df[,i+6]
    
    pe_i <- 1/num_risk_t * (pred_data_surv$ind1_i %*% (1 - pd_u_t_clustered_i)^2 + 
                              pred_data_surv$ind2_i %*% pd_u_t_clustered_i^2 + 
                              pred_data_surv$ind3_i %*% (pd_u_t_clustered_censored_df_i*(1 - pd_u_t_clustered_i)^2 + 
                                                           (1-pd_u_t_clustered_censored_df_i)*pd_u_t_clustered_i^2))
    pe <- c(pe, pe_i)
  }
  pe
  
  #### AUC (within group) ####
  AUC_w_matrix <- c()
  AUC_w <- c()
  num_within_cluster_comp_pairs <- c()
  for (u in 1:length(survTimes)){
    u_i <- survTimes[u]
    pred_data_surv$alive_ind_i <- ifelse(pred_data_surv$eventtime_year > u_i, 1, 0)
    pred_data_surv$event_ind_i <- 1 - pred_data_surv$alive_ind_i
    
    AUC_w_ui_vec <- c()
    
    num_within_cluster_comp_pairs_u <- 0
    
    for (id in unique(pred_data_surv$researchID)) {
      pred_data_surv_i <- pred_data_surv %>% filter(researchID==id)
      if (nrow(pred_data_surv_i) > 1 ) {
        pd_u_t_clustered_i <- pd_u_t_clustered_df %>% filter(researchID==id)
        pd_u_t_clustered_u_i <- pd_u_t_clustered_i[,u+6]
        pd_u_t_clustered_censored_df_i <- pd_u_t_clustered_censored_df%>%filter(researchID==id)
        pd_u_t_clustered_censored_df_u_i <- pd_u_t_clustered_censored_df_i[,u+6]
        
        
        pairs <- combn(pred_data_surv_i$id_tooth, 2)
        Tij1 <- pred_data_surv_i$eventtime_year[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        Tij2 <- pred_data_surv_i$eventtime_year[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        dij1 <- pred_data_surv_i$status[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        dij2 <- pred_data_surv_i$status[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_ij1 <- pd_u_t_clustered_u_i[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_ij2 <- pd_u_t_clustered_u_i[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_censored_ij1 <- pd_u_t_clustered_censored_df_u_i[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_censored_ij2 <- pd_u_t_clustered_censored_df_u_i[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        
        
        ind1 <- (Tij1 <= u_i & dij1 == 1) & Tij2 > u_i
        ind2 <- (Tij1 <= u_i & dij1 == 0) & Tij2 > u_i
        ind3 <- (Tij1 <= u_i & dij1 == 1) & (Tij2 > Tij1 & Tij2 <= u_i & dij2 == 0)
        ind4 <- (Tij1 <= u_i & dij1 == 0) & (Tij2 > Tij1 & Tij2 <= u_i & dij2 == 0)      
        weight1 <- 1
        weight2 <- 1 - pd_u_t_clustered_censored_ij1
        weight3 <- pd_u_t_clustered_censored_ij2
        weight4 <- (1 - pd_u_t_clustered_censored_ij1) * pd_u_t_clustered_censored_ij2
        
        ind_order_ij12 <- ifelse(pd_u_t_clustered_ij1 < pd_u_t_clustered_ij2, 1, 0)
        
        AUC_w_ui1 <- sum(ind_order_ij12 * ind1)/sum(ind1)
        AUC_w_ui1 <- ifelse(is.na(AUC_w_ui1), 0, AUC_w_ui1)
        AUC_w_ui2 <-  sum(ind_order_ij12 * ind2 * weight2)/sum(ind2 * weight2)
        AUC_w_ui2 <- ifelse(is.na(AUC_w_ui2), 0, AUC_w_ui2)
        AUC_w_ui3 <-  sum(ind_order_ij12 * ind3 * weight3)/sum(ind3 * weight3)
        AUC_w_ui3 <- ifelse(is.na(AUC_w_ui3), 0, AUC_w_ui3)
        AUC_w_ui4 <-  sum(ind_order_ij12 * ind4 * weight4)/sum(ind4 * weight4)
        AUC_w_ui4 <- ifelse(is.na(AUC_w_ui4), 0, AUC_w_ui4)
        
        AUC_w_ui <- AUC_w_ui1 + AUC_w_ui2 + AUC_w_ui3 + AUC_w_ui4
        AUC_w_ui_vec <- c(AUC_w_ui_vec, AUC_w_ui)
        
        num_within_cluster_comp_pairs_u <- num_within_cluster_comp_pairs_u + sum(ind1)
        
      }
    }
    AUC_w_matrix <- cbind(AUC_w_matrix, AUC_w_ui_vec)
    AUC_w_u <- mean(AUC_w_ui_vec)
    AUC_w <- c(AUC_w, AUC_w_u)
    
    num_within_cluster_comp_pairs <- c(num_within_cluster_comp_pairs, num_within_cluster_comp_pairs_u)
  }
  
  #### AUC (between group) #### 
  t0=Sys.time()
  AUC_b_matrix <- c()
  AUC_b <- c()
  num_between_cluster_comp_pairs <- c()
  for (u in 1:length(survTimes)){
    u_i <- survTimes[u]
    pred_data_surv$alive_ind_i <- ifelse(pred_data_surv$eventtime_year > u_i, 1, 0)
    pred_data_surv$event_ind_i <- 1 - pred_data_surv$alive_ind_i
    
    
    AUC_b1_ui_vec1 <- c()
    AUC_b1_ui_vec2 <- c()
    AUC_b2_ui_vec1 <- c()
    AUC_b2_ui_vec2 <- c()
    AUC_b3_ui_vec1 <- c()
    AUC_b3_ui_vec2 <- c()
    AUC_b4_ui_vec1 <- c()
    AUC_b4_ui_vec2 <- c()
    
    
    pairs <- combn(unique(pred_data_surv$researchID), 2)
    i1_all <- pairs[1,]
    i2_all <- pairs[2,]
    
    num_between_cluster_comp_pairs_u <- 0
    
    
    for (id in 1:length(i1_all)){
      i1 <- i1_all[id]
      i2 <- i2_all[id]
      pred_data_surv_i1 <- pred_data_surv %>% filter(researchID==i1)
      pred_data_surv_i2 <- pred_data_surv %>% filter(researchID==i2)
      
      ij1 <- pred_data_surv_i1$id_tooth
      ij2 <- pred_data_surv_i2$id_tooth
      ij12_comb <- expand.grid(ij1, ij2)
      names(ij12_comb) <- c("ij1", "ij2")
      
      pred_data_surv_ij1 <- pred_data_surv %>% filter(id_tooth %in% ij1)
      pred_data_surv_ij2 <- pred_data_surv %>% filter(id_tooth %in% ij2)
      
      pd_u_t_clustered_ij1 <- pd_u_t_clustered_df %>% filter(id_tooth %in% ij1)
      pd_u_t_clustered_u_ij1 <- pd_u_t_clustered_ij1[,c(1, u+6)]
      pd_u_t_clustered_censored_df_ij1 <- pd_u_t_clustered_censored_df %>% filter(id_tooth %in% ij1)
      pd_u_t_clustered_censored_df_u_ij1 <- pd_u_t_clustered_censored_df_ij1[,c(1, u+6)]
      
      pd_u_t_clustered_ij2 <- pd_u_t_clustered_df %>% filter(id_tooth %in% ij2)
      pd_u_t_clustered_u_ij2 <- pd_u_t_clustered_ij2[,c(1, u+6)]
      pd_u_t_clustered_censored_df_ij2 <- pd_u_t_clustered_censored_df %>% filter(id_tooth %in% ij2)
      pd_u_t_clustered_censored_df_u_ij2 <- pd_u_t_clustered_censored_df_ij2[,c(1, u+6)]
      
      
      Tij1 <- pred_data_surv_ij1$eventtime_year[match(ij12_comb[, 1], pred_data_surv_ij1$id_tooth)]
      Tij2 <- pred_data_surv_ij2$eventtime_year[match(ij12_comb[, 2], pred_data_surv_ij2$id_tooth)]
      dij1 <- pred_data_surv_ij1$status[match(ij12_comb[, 1], pred_data_surv_ij1$id_tooth)]
      dij2 <- pred_data_surv_ij2$status[match(ij12_comb[, 2], pred_data_surv_ij2$id_tooth)]
      pd_u_t_clustered_ij1 <- pd_u_t_clustered_u_ij1[,2][match(ij12_comb[, 1], pd_u_t_clustered_u_ij1$id_tooth)]
      pd_u_t_clustered_ij2 <- pd_u_t_clustered_u_ij2[,2][match(ij12_comb[, 2], pd_u_t_clustered_u_ij2$id_tooth)]
      pd_u_t_clustered_censored_ij1 <- pd_u_t_clustered_censored_df_u_ij1[,2][match(ij12_comb[, 1], pd_u_t_clustered_censored_df_u_ij1$id_tooth)]
      pd_u_t_clustered_censored_ij2 <- pd_u_t_clustered_censored_df_u_ij2[,2][match(ij12_comb[, 2], pd_u_t_clustered_censored_df_u_ij2$id_tooth)]
      
      
      ind1 <- (Tij1 <= u_i & dij1 == 1) & Tij2 > u_i
      ind2 <- (Tij1 <= u_i & dij1 == 0) & Tij2 > u_i
      ind3 <- (Tij1 <= u_i & dij1 == 1) & (Tij2 > Tij1 & Tij2 <= u_i & dij2 == 0)
      ind4 <- (Tij1 <= u_i & dij1 == 0) & (Tij2 > Tij1 & Tij2 <= u_i & dij2 == 0)      
      
      weight1 <- 1
      weight2 <- 1 - pd_u_t_clustered_censored_ij1
      weight3 <- pd_u_t_clustered_censored_ij2
      weight4 <- (1 - pd_u_t_clustered_censored_ij1) * pd_u_t_clustered_censored_ij2
      
      
      ind_order_ij12 <- ifelse(pd_u_t_clustered_ij1 < pd_u_t_clustered_ij2, 1, 0)
      
      AUC_b1_ui1 <- sum(ind_order_ij12 * ind1)
      AUC_b1_ui2 <- sum(ind1)
      
      AUC_b2_ui1 <-  sum(ind_order_ij12 * ind2 * weight2)
      AUC_b2_ui1 <- ifelse(is.na(AUC_b2_ui1), 0, AUC_b2_ui1)
      AUC_b2_ui2 <- sum(ind2 * weight2)
      
      AUC_b3_ui1 <-  sum(ind_order_ij12 * ind3 * weight3)
      AUC_b3_ui1 <- ifelse(is.na(AUC_b3_ui1), 0, AUC_b3_ui1)
      AUC_b3_ui2 <- sum(ind3 * weight3)
      
      AUC_b4_ui1 <-  sum(ind_order_ij12 * ind4 * weight4)
      AUC_b4_ui1 <- ifelse(is.na(AUC_b4_ui1), 0, AUC_b4_ui1)
      AUC_b4_ui2 <- sum(ind4 * weight4)
      
      AUC_b1_ui_vec1 <- c(AUC_b1_ui_vec1, AUC_b1_ui1)
      AUC_b1_ui_vec2 <- c(AUC_b1_ui_vec2, AUC_b1_ui2)
      AUC_b2_ui_vec1 <- c(AUC_b2_ui_vec1, AUC_b2_ui1)
      AUC_b2_ui_vec2 <- c(AUC_b2_ui_vec2, AUC_b2_ui2)
      AUC_b3_ui_vec1 <- c(AUC_b3_ui_vec1, AUC_b3_ui1)
      AUC_b3_ui_vec2 <- c(AUC_b3_ui_vec2, AUC_b3_ui2)
      AUC_b4_ui_vec1 <- c(AUC_b4_ui_vec1, AUC_b4_ui1)
      AUC_b4_ui_vec2 <- c(AUC_b4_ui_vec2, AUC_b4_ui2)
      
      num_between_cluster_comp_pairs_u <- num_between_cluster_comp_pairs_u + sum(ind1)
    }
    AUC_b1_u <- sum(AUC_b1_ui_vec1)/sum(AUC_b1_ui_vec2)
    AUC_b1_u <- ifelse(is.na(AUC_b1_u), 0, AUC_b1_u)
    
    AUC_b2_u <- sum(AUC_b2_ui_vec1)/sum(AUC_b2_ui_vec2)
    AUC_b2_u <- ifelse(is.na(AUC_b2_u), 0, AUC_b2_u)
    
    AUC_b3_u <- sum(AUC_b3_ui_vec1)/sum(AUC_b3_ui_vec2)
    AUC_b3_u <- ifelse(is.na(AUC_b3_u), 0, AUC_b3_u)
    
    AUC_b4_u <- sum(AUC_b4_ui_vec1)/sum(AUC_b4_ui_vec2)
    AUC_b4_u <- ifelse(is.na(AUC_b4_u), 0, AUC_b4_u)
    
    AUC_b_u <- AUC_b1_u + AUC_b2_u + AUC_b3_u + AUC_b4_u
    
    num_between_cluster_comp_pairs <- c(num_between_cluster_comp_pairs, num_between_cluster_comp_pairs_u)
    AUC_b <- c(AUC_b, AUC_b_u)
    AUC_b <- AUC_b/4
  }
  t1=Sys.time()
  auc_b_time = t1-t0
  
  
  AUC_o <- c()
  for (u in 1:length(survTimes)){
    num_within_cluster_comp_pairs_u <- num_within_cluster_comp_pairs[u]
    num_between_cluster_comp_pairs_u <- num_between_cluster_comp_pairs[u]
    num_overall_cluster_comp_pairs_u <- num_within_cluster_comp_pairs_u+num_between_cluster_comp_pairs_u
    AUC_w_u <- AUC_w[u]
    AUC_b_u <- AUC_b[u]
    
    AUC_o_u <- num_within_cluster_comp_pairs_u/num_overall_cluster_comp_pairs_u*AUC_w_u+
      num_between_cluster_comp_pairs_u/num_overall_cluster_comp_pairs_u*AUC_b_u
    AUC_o <- c(AUC_o, AUC_o_u)
  }
  
  
  
  #### AUC (IPCW, within group)  ####
  
  AUC_ipcw_w_matrix <- c()
  AUC_ipcw_w <- c()
  num_within_cluster_comp_pairs <- c()
  
  cens_data <- data.frame(Time = pred_data_surv$eventtime_year, cens_ind = 1 - pred_data_surv$status)
  censoring_dist <- survfit(Surv(Time, cens_ind) ~ 1, data = cens_data)

  for (u in 1:length(survTimes)){
    AUC_w_ui_vec <- c()
    
    num_within_cluster_comp_pairs_u <- 0
    
    
    Thoriz_u <- survTimes[u]
    
    for (id in unique(pred_data_surv$researchID)) {
      pred_data_surv_i <- pred_data_surv %>% filter(researchID==id)
      if (nrow(pred_data_surv_i) > 1 ) {
        pd_u_t_clustered_i <- pd_u_t_clustered_df %>% filter(researchID==id)
        pd_u_t_clustered_u_i <- pd_u_t_clustered_i[,u+6]
        pd_u_t_clustered_censored_df_i <- pd_u_t_clustered_censored_df%>%filter(researchID==id)
        pd_u_t_clustered_censored_df_u_i <- pd_u_t_clustered_censored_df_i[,u+6]
        
        
        pairs <- combn(pred_data_surv_i$id_tooth, 2)
        Tij1 <- pred_data_surv_i$eventtime_year[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        Tij2 <- pred_data_surv_i$eventtime_year[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        ij1 <- pred_data_surv_i$id_tooth[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        ij2 <- pred_data_surv_i$id_tooth[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        
        dij1 <- pred_data_surv_i$status[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        dij2 <- pred_data_surv_i$status[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_ij1 <- pd_u_t_clustered_u_i[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_ij2 <- pd_u_t_clustered_u_i[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_censored_ij1 <- pd_u_t_clustered_censored_df_u_i[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        pd_u_t_clustered_censored_ij2 <- pd_u_t_clustered_censored_df_u_i[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        
        comparable_pirs <- (Tij1 <= Thoriz_u & dij1 == 1) & Tij2 > Thoriz_u
        
        
        # subjects who had the event in the interval (Tstart, Thoriz)
        pred_data_surv_i$ind1 <- pred_data_surv_i$eventtime_year < Thoriz_u & pred_data_surv_i$status == 1
        ind1_1 <- pred_data_surv_i$ind1[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        ind1_2 <- pred_data_surv_i$ind1[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        
        
        death_dij1 <- pred_data_surv_i$ind1[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        death_dij2 <- pred_data_surv_i$ind1[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        
        # subjects who had the event after Thoriz_i
        pred_data_surv_i$ind2 <- pred_data_surv_i$eventtime_year > Thoriz_u
        ind2_1 <- pred_data_surv_i$ind2[match(pairs[1, ], pred_data_surv_i$id_tooth)]
        ind2_2 <- pred_data_surv_i$ind2[match(pairs[2, ], pred_data_surv_i$id_tooth)]
        
        # subjects who were censored in the interval (Tstart, Thoriz)
        pred_data_surv_i$ind3 <- pred_data_surv_i$eventtime_year < Thoriz_u & pred_data_surv_i$status == 0
        
        weights1 <- numeric(length(Tij1))
        weights2 <- numeric(length(Tij1))
        
        if (sum(ind1_1)>=1){
          ss1 <- summary(censoring_dist, times = Tij1[ind1_1])$surv/summary(censoring_dist, times = pred.time)$surv
          weights1[ind1_1] <- 1 / ss1
        }
        if (sum(ind1_2)>=1){
          ss2 <- summary(censoring_dist, times = Tij2[ind1_2])$surv/summary(censoring_dist, times = pred.time)$surv
          weights2[ind1_2] <- 1 / ss2
        }
        
        
        weights1[ind2_1] <- 1 / (summary(censoring_dist, times = Thoriz_u)$surv/summary(censoring_dist, times = pred.time)$surv ) 
        weights2[ind2_2] <- 1 / (summary(censoring_dist, times = Thoriz_u)$surv/summary(censoring_dist, times = pred.time)$surv ) 
        
        
        ind_order_ij12 <- ifelse(pd_u_t_clustered_ij1 < pd_u_t_clustered_ij2, 1, 0)
        
        # the denominator is 0 if no event occured
        
        AUC_w_ui <- sum(ind_order_ij12*death_dij1*(1-death_dij2)*weights1*weights2)/sum(death_dij1*(1-death_dij2)*weights1*weights2)
        AUC_w_ui_vec <- c(AUC_w_ui_vec, AUC_w_ui)
        
        num_within_cluster_comp_pairs_u <- num_within_cluster_comp_pairs_u + sum(comparable_pirs)
        
      }
    }
    AUC_ipcw_w_matrix <- cbind(AUC_ipcw_w_matrix, AUC_w_ui_vec)
    AUC_w_u <- mean(AUC_w_ui_vec, na.rm=T)
    AUC_ipcw_w <- c(AUC_ipcw_w, AUC_w_u)
    
    num_within_cluster_comp_pairs <- c(num_within_cluster_comp_pairs, num_within_cluster_comp_pairs_u)
  }
  
  
  #### AUC (IPCW, between group) #### 
  AUC_ipcw_b_matrix <- c()
  AUC_ipcw_b <- c()
  num_between_cluster_comp_pairs <- c()
  for (u in 1:length(survTimes)){
    u_i <- survTimes[u]
    pred_data_surv$alive_ind_i <- ifelse(pred_data_surv$eventtime_year > u_i, 1, 0)
    pred_data_surv$event_ind_i <- 1 - pred_data_surv$alive_ind_i
    
    pairs <- combn(unique(pred_data_surv$researchID), 2)
    i1_all <- pairs[1,]
    i2_all <- pairs[2,]
    
    num_between_cluster_comp_pairs_u <- 0
    
    Thoriz_u <- survTimes[u]
    
    AUC_ipcw_ui_vec1 <-c()
    AUC_ipcw_ui_vec2 <- c()
    
    
    for (id in 1:length(i1_all)){
      i1 <- i1_all[id]
      i2 <- i2_all[id]
      pred_data_surv_i1 <- pred_data_surv %>% filter(researchID==i1)
      pred_data_surv_i2 <- pred_data_surv %>% filter(researchID==i2)
      
      ij1 <- pred_data_surv_i1$id_tooth
      ij2 <- pred_data_surv_i2$id_tooth
      ij12_comb <- expand.grid(ij1, ij2)
      names(ij12_comb) <- c("ij1", "ij2")
      
      pred_data_surv_ij1 <- pred_data_surv %>% filter(id_tooth %in% ij1)
      pred_data_surv_ij2 <- pred_data_surv %>% filter(id_tooth %in% ij2)
      
      pd_u_t_clustered_ij1 <- pd_u_t_clustered_df %>% filter(id_tooth %in% ij1)
      pd_u_t_clustered_u_ij1 <- pd_u_t_clustered_ij1[,c(1, u+6)]
      pd_u_t_clustered_censored_df_ij1 <- pd_u_t_clustered_censored_df %>% filter(id_tooth %in% ij1)
      pd_u_t_clustered_censored_df_u_ij1 <- pd_u_t_clustered_censored_df_ij1[,c(1, u+6)]
      
      pd_u_t_clustered_ij2 <- pd_u_t_clustered_df %>% filter(id_tooth %in% ij2)
      pd_u_t_clustered_u_ij2 <- pd_u_t_clustered_ij2[,c(1, u+6)]
      pd_u_t_clustered_censored_df_ij2 <- pd_u_t_clustered_censored_df %>% filter(id_tooth %in% ij2)
      pd_u_t_clustered_censored_df_u_ij2 <- pd_u_t_clustered_censored_df_ij2[,c(1, u+6)]
      
      
      Tij1 <- pred_data_surv_ij1$eventtime_year[match(ij12_comb[, 1], pred_data_surv_ij1$id_tooth)]
      Tij2 <- pred_data_surv_ij2$eventtime_year[match(ij12_comb[, 2], pred_data_surv_ij2$id_tooth)]
      dij1 <- pred_data_surv_ij1$status[match(ij12_comb[, 1], pred_data_surv_ij1$id_tooth)]
      dij2 <- pred_data_surv_ij2$status[match(ij12_comb[, 2], pred_data_surv_ij2$id_tooth)]
      pd_u_t_clustered_ij1 <- pd_u_t_clustered_u_ij1[,2][match(ij12_comb[, 1], pd_u_t_clustered_u_ij1$id_tooth)]
      pd_u_t_clustered_ij2 <- pd_u_t_clustered_u_ij2[,2][match(ij12_comb[, 2], pd_u_t_clustered_u_ij2$id_tooth)]
      pd_u_t_clustered_censored_ij1 <- pd_u_t_clustered_censored_df_u_ij1[,2][match(ij12_comb[, 1], pd_u_t_clustered_censored_df_u_ij1$id_tooth)]
      pd_u_t_clustered_censored_ij2 <- pd_u_t_clustered_censored_df_u_ij2[,2][match(ij12_comb[, 2], pd_u_t_clustered_censored_df_u_ij2$id_tooth)]
      
      
      comparable_pirs <- (Tij1 <= Thoriz_u & dij1 == 1) & Tij2 > Thoriz_u
      
      
      # subjects who had the event in the interval (Tstart, Thoriz)
      pred_data_surv_ij1$ind1 <- pred_data_surv_ij1$eventtime_year < Thoriz_u & pred_data_surv_ij1$status == 1
      pred_data_surv_ij2$ind1 <- pred_data_surv_ij2$eventtime_year < Thoriz_u & pred_data_surv_ij2$status == 1
      ind1_1 <- pred_data_surv_ij1$ind1[match(ij12_comb[, 1], pred_data_surv_ij1$id_tooth)]
      ind1_2 <- pred_data_surv_ij2$ind1[match(ij12_comb[, 2], pred_data_surv_ij2$id_tooth)]
      
      
      death_dij1 <- pred_data_surv_ij1$ind1[match(ij12_comb[, 1],  pred_data_surv_ij1$id_tooth)]
      death_dij2 <- pred_data_surv_ij2$ind1[match(ij12_comb[, 2],  pred_data_surv_ij2$id_tooth)]
      
      # subjects who had the event after Thoriz_i
      pred_data_surv_ij1$ind2 <- pred_data_surv_ij1$eventtime_year > Thoriz_u
      pred_data_surv_ij2$ind2 <- pred_data_surv_ij2$eventtime_year > Thoriz_u
      
      ind2_1 <- pred_data_surv_ij1$ind2[match(ij12_comb[, 1], pred_data_surv_ij1$id_tooth)]
      ind2_2 <- pred_data_surv_ij2$ind2[match(ij12_comb[, 2], pred_data_surv_ij2$id_tooth)]
      
      # subjects who were censored in the interval (Tstart, Thoriz)
      # pred_data_surv_i$ind3 <- pred_data_surv_i$eventtime < Thoriz_u & pred_data_surv_i$status == 0
      
      weights1 <- numeric(length(Tij1))
      weights2 <- numeric(length(Tij1))
      
      if (sum(pred_data_surv_ij1$ind1)>=1){
        ss1 <- summary(censoring_dist, times = Tij1[ind1_1])$surv/summary(censoring_dist, times = pred.time)$surv
        weights1[ind1_1] <- 1 / ss1
      }
      
      if (sum(pred_data_surv_ij2$ind1)>=1){
        ss2 <- summary(censoring_dist, times = Tij2[ind1_2])$surv/summary(censoring_dist, times = pred.time)$surv
        weights2[ind1_2] <- 1 / ss2
      }
      
      
      weights1[ind2_1] <- 1 / (summary(censoring_dist, times = Thoriz_u)$surv/summary(censoring_dist, times = pred.time)$surv ) 
      weights2[ind2_2] <- 1 / (summary(censoring_dist, times = Thoriz_u)$surv/summary(censoring_dist, times = pred.time)$surv ) 
      
      
      ind_order_ij12 <- ifelse(pd_u_t_clustered_ij1 < pd_u_t_clustered_ij2, 1, 0)
      
      AUC_b1_ui1 <- sum(ind_order_ij12*death_dij1*(1-death_dij2)*weights1*weights2)
      AUC_b1_ui2 <- sum(death_dij1*(1-death_dij2)*weights1*weights2)
      
      AUC_ipcw_ui_vec1 <- c(AUC_ipcw_ui_vec1, AUC_b1_ui1)
      AUC_ipcw_ui_vec2 <- c(AUC_ipcw_ui_vec2, AUC_b1_ui2)
      
      num_between_cluster_comp_pairs_u <- num_between_cluster_comp_pairs_u + sum(comparable_pirs)
    }
    AUC_b1_u <- sum(AUC_ipcw_ui_vec1)/sum(AUC_ipcw_ui_vec2)
    AUC_b1_u <- ifelse(is.na(AUC_b1_u), 0, AUC_b1_u)
    
    AUC_b_u <- AUC_b1_u
    
    num_between_cluster_comp_pairs <- c(num_between_cluster_comp_pairs, num_between_cluster_comp_pairs_u)
    AUC_ipcw_b <- c(AUC_ipcw_b, AUC_b_u)
  }
  
  
  AUC_ipcw_o <- c()
  for (u in 1:length(survTimes)){
    num_within_cluster_comp_pairs_u <- num_within_cluster_comp_pairs[u]
    num_between_cluster_comp_pairs_u <- num_between_cluster_comp_pairs[u]
    num_overall_cluster_comp_pairs_u <- num_within_cluster_comp_pairs_u+num_between_cluster_comp_pairs_u
    AUC_w_u <- AUC_ipcw_w[u]
    AUC_b_u <- AUC_ipcw_b[u]
    
    AUC_o_u <- num_within_cluster_comp_pairs_u/num_overall_cluster_comp_pairs_u*AUC_w_u+
      num_between_cluster_comp_pairs_u/num_overall_cluster_comp_pairs_u*AUC_b_u
    AUC_ipcw_o <- c(AUC_ipcw_o, AUC_o_u)
  }
  
  
  #### AUC (IPCW, timeROC)####
 
  timeROC_AUC <- c()
  for (i in 1:(length(survTimes))-1) {
    u_i <- survTimes[i+1]
    window_i <- u_i - pred.time
    timeROC_i <- timeROC(T=pred_data_surv$eventtime_year-pred.time, 
                         delta = pred_data_surv$status, 
                         marker = 1-pd_u_t_clustered_df[,(7+i)],
                         times =fwindow[i+1]-pred.time,
                         cause = 1)
    timeROC_AUC_i <- timeROC_i$AUC[[2]]
    timeROC_AUC <- c(timeROC_AUC, timeROC_AUC_i)
  }
  
  
  #### return output ####
  AUC_PE_result <- cbind(
    times = survTimes,
    "PE_model_based" = pe,
    "AUC_model_based" = AUC_o,
    "PE_pec_IPCW" = pec_bs, 
    "AUC_frailtyROC" = frailtyROC_AUC, 
    "AUC_timeROC" = timeROC_AUC,
    "AUC_cIPCW" = AUC_ipcw_o
  )
  result_list <- list(AUC_PE_result,
                      pd_u_t_clustered_df,
                      pd_u_t_clustered_censored_df)
  return(result_list)
}





roc_plot_clusteredJM <- function(pd_u_t_clustered_df, 
                                 pred.time = predtime,
                                 Thoriz,
                                 optimal_cutoff = c("", "F1score", "Youden"), 
                                 legend = T){
  
  #### tvROC ####
  Thoriz_name <- paste0("year", Thoriz)
  u = which(colnames(pd_u_t_clustered_df) %in% Thoriz_name)
  
  p_u_t_ij <- 1-pd_u_t_clustered_df[,u]
  q_u_t_ij <- matrix(1-p_u_t_ij, ncol=1)
  rownames(q_u_t_ij) <- pd_u_t_clustered_df$id_tooth
  id <- pd_u_t_clustered_df[["id_tooth"]]
  Time <- pd_u_t_clustered_df[["eventtime_year"]]
  event <- pd_u_t_clustered_df[["status"]]
  f <- factor(id, levels = unique(id))
  Time <- tapply(Time, f, tail, 1L)
  event <- tapply(event, f, tail, 1L)
  names(Time) <- names(event) <- as.character(unique(id))
  thrs <- seq(0, 1, length = 501)
  Check <- lapply(seq_len(ncol(q_u_t_ij)), function (i) outer(q_u_t_ij[, i], thrs, "<"))
  Check_mean <- outer(matrix(1-p_u_t_ij, ncol=1), thrs, "<")    

  # subjects who died before Thoriz
  ind1 <- Time < Thoriz & event == 1
  # subjects who were censored in the interval (Tstart, Thoriz)
  ind2 <- Time < Thoriz & event == 0
  ind <- ind1 | ind2
  
  if (any(ind2)) {
    nams <- names(ind2[ind2])
    p_u_t_ij2 <- 1-pd_u_t_clustered_df[pd_u_t_clustered_df[["id_tooth"]]%in%nams, u]
    f <- factor(nams, levels = unique(nams))
    names(p_u_t_ij2) <- f
    p_u_t_ij2 <- tapply(p_u_t_ij2, f, tail, 1)
    nams2 <- names(ind2[ind2])
    ind[ind2] <- ind[ind2] * p_u_t_ij2[nams2]
  }
  
  # calculate sensitivity and specificity
  ntp <- lapply(Check, function (x) colSums(x * c(ind)))
  tp <- lapply(ntp, function (x) x / sum(ind))
  nTP <- colSums(Check_mean * c(ind))
  nFN <- sum(ind) - nTP
  TP <- nTP / sum(ind)
  nfp <- lapply(Check, function (x) colSums(x * c(1 - ind)))
  fp <- lapply(nfp, function (x) x / sum(1 - ind))
  nFP <- rowMeans(do.call("cbind", nfp))
  nFP <- colSums(Check_mean * c(1 - ind))
  nTN <- sum(1 - ind) - nFP
  FP <- nFP / sum(1 - ind)
  
  
  f1score <- 2 * nTP / (2 * nTP + nFN + nFP)
  F1score <- median(thrs[f1score == max(f1score)])
  youden <- TP - FP
  Youden <- median(thrs[youden == max(youden)])
  tvROCout <- list(TP = TP, FP = FP, nTP = nTP, nFN = nFN, nFP = nFP, nTN = nTN,
                   tp = do.call("cbind", tp), fp = do.call("cbind", fp),
                   thrs = thrs, F1score = F1score, Youden = Youden,
                   Tstart = pred.time, Thoriz = Thoriz, nr = length(unique(id)))
  class(tvROCout) <- "tvROC"
  tvROCout
  
  plot(FP, TP, type = "l")
  abline(a = 0, b = 1, lty = 3)
  optimal_cutoff <- match.arg(optimal_cutoff)
  if (optimal_cutoff == "F1score")
    abline(v = F1score, lty = 3, lwd = 2, col = 2)
  if (optimal_cutoff == "Youden")
    abline(v = Youden, lty = 3, lwd = 2, col = 2)
  if (legend) {
    legend("bottomright", c(paste("At time:", round(Thoriz, 1), "\n"),
                            paste("Using information up to time:",
                                  round(pred.time, 1))),
           bty = "n")
  }
  
  
  auc <- sum(0.5 * diff(as.vector(FP)) * (TP[-1L] + TP[-length(TP)]), na.rm = TRUE)
  
  
  
  
  
  roc_dat_t <- subset(roc_dat, time > t0 | status == 1)
  
  case <- roc_dat_t$status == 1 & roc_dat_t$time <= u
  control <- roc_dat_t$time > u
  
  thresholds <- sort(unique(roc_dat_t$risk))
  
  roc_curve <- data.frame(
    threshold = thresholds,
    sensitivity = sapply(thresholds, function(c) {
      mean(roc_dat_t$risk[case] > c)
    }),
    specificity = sapply(thresholds, function(c) {
      mean(roc_dat_t$risk[control] <= c)
    })
  )
  
  plot(
    1 - roc_curve$specificity,
    roc_curve$sensitivity,
    type = "l",
    xlab = "1 - Specificity",
    ylab = "Sensitivity",
    main = paste0("Time-dependent ROC: (", t0, ", ", u, "]")
  )
  abline(0, 1, lty = 2)
}








survpred_JM <- function(object=clustered_jm.est, 
                        draws_df_cluster=draws_df_cluster_chain34,
                        surv_data, 
                        long_data, 
                        type = c("SurvProb", "Density"), 
                        pred_id, 
                        pred_j, 
                        idVar = "id", 
                        simulate = TRUE, 
                        survTimes = fwindow, 
                        pred.time = predtime, 
                        LeftTrunc_var = NULL, 
                        MCMCnum = 200L,
                        MCMCnum_b = 200L,
                        CI.levels = c(0.025, 0.975), 
                        log = FALSE, 
                        scale = 1.6, 
                        init.b = NULL, 
                        seed = 1L) {
  
  newdata <- long_data[long_data$i == pred_id & long_data$j == pred_j, ]
  weight = rep(1, nrow(newdata))
  
  #### prepare data ####
  timeVar <- "obstime"
  hasScale <- T # biomarker is continuous, has measurement error
  gammas <- clustered_jm.est[c("gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]")]
  betas1 <- clustered_jm.est[c("beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[5]","beta1[6]","beta1[7]")]
  betas2 <- clustered_jm.est[c("beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[5]","beta2[6]","beta2[7]")]
  xi <- clustered_jm.est[c("xi")]
  alpha1 <- clustered_jm.est[c("alpha1")]
  alpha2 <- clustered_jm.est[c("alpha2")]
  sigma_e <- clustered_jm.est["sigma_e"]
  sigma_bi <-clustered_jm.est[c("sigma_bi[1]","sigma_bi[2]")]
  sigma_bij <- clustered_jm.est[c("sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]")]
  sigma_ri <- clustered_jm.est["sigma_ri"]
  
  R_bi <- clustered_jm.est["R_bi[1,2]"]
  R_bij <- clustered_jm.est[c("R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]")]
  
  
  var_ri = sigma_ri^2
  var_e = sigma_e^2
  
  length_bi <- length(sigma_bi)
  Cor_bi <- diag(1, length_bi)
  Cor_bi[lower.tri(Cor_bi)] <- R_bi
  Cor_bi <- Cor_bi + t(Cor_bi) - diag(diag(Cor_bi))  # Make symmetric
  var_bi <- diag(sigma_bi) %*% Cor_bi %*%  diag(sigma_bi)
  
  length_bij <- length(sigma_bij)
  Cor_bij <- diag(1, length_bij)
  Cor_bij[lower.tri(Cor_bij)] <- R_bij
  Cor_bij <- Cor_bij + t(Cor_bij) - diag(diag(Cor_bij))  # Make symmetric
  var_bij <- diag(sigma_bij) %*% Cor_bij %*%  diag(sigma_bij)
  L_bi_matrix <- t(chol(var_bi))
  L_bij_matrix <- t(chol(var_bij))
  

  # create an indicator for censoring before the landmark time
  surv_data <- surv_data %>% 
    mutate(censor_ind = ifelse((eventtime <= pred.time & status == 0), 1, 0))
  
  ## filter ij's who are at risk and not being censored at pred.time
  pred_data_surv <- surv_data %>% filter(eventtime > pred.time)
  pred_data_long <- long_data %>% 
    filter((id_tooth %in% pred_data_surv$id_tooth) & (followup_year <= pred.time)) %>% group_by(id_tooth) 
  
  ## filter ij's who are censored between pred.time and the prediction window
  pred_data_surv_censored <- surv_data %>% filter(eventtime_year <= survTimes[length(survTimes)] & status == 0)
  pred_data_long_censored <- long_data %>% 
    filter((id_tooth %in% pred_data_surv_censored$id_tooth)) %>% group_by(id_tooth) 
  
  n.tp <- length(pred.time)
  ncz <- 7
  
  #  environment(log.posterior.b) <- environment()
  
  
  
  #### useful functions ####
  
  rmvt <- function (n, mu, Sigma, df) {
    p <- length(mu)
    if (is.list(Sigma)) {
      ev <- Sigma$values
      evec <- Sigma$vectors
    } else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
    }
    X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                               matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
    if (n == 1L) drop(X) else t.default(X)
  }
  
  dmvt <- function (x, mu, Sigma = NULL, invSigma = NULL, df, log = FALSE, prop = TRUE) {
    if (!is.numeric(x)) 
      stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
      x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
      }
      #      if (!all(ev >= -1e-06 * abs(ev[1]))) 
      #        stop("'Sigma' is not positive definite")
      invSigma <- evec %*% (t(evec)/ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
      fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + 
                                                                log(df)) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
      if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else as.vector(((1 + quad)^(-(df + p)/2)))
    }
  }
  
  densLong_norm <- function (y, eta.y, scale, log = FALSE, data) {
    dnorm(x = y, mean = eta.y, sd = scale, log = log)
  }
  
  dmvnorm <- function (x, mu, Sigma = NULL, invSigma = NULL, log = FALSE, prop = TRUE) {
    if (!is.matrix(x))
      x <- rbind(x)
    if (is.matrix(mu)) {
      if (nrow(mu) != nrow(x))
        stop("incorrect dimensions for 'mu'.")
      p <- ncol(mu)
    } else {
      p <- length(mu)
      mu <- rep(mu, each = nrow(x))
    }
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors            
      }
      invSigma <- evec %*% (t(evec) / ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- - determinant(as.matrix(invSigma))$modulus
    }
    ss <- x - mu
    quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
    if (!prop)
      fact <- - 0.5 * (p * log(2 * pi) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - quad) else as.vector(- quad)
    } else {
      if (!prop) as.vector(exp(fact - quad)) else as.vector(exp(- quad))
    }
  }
  
  densLong_binomial <- function (y, eta.y, scale, log = FALSE, data) {
    y*eta.y - log(1+exp(eta.y))
  }
  
  
  
 
  S_u_t <- function(t1, t0, X_surv, betas1, betas2, gammas, 
                    xi, alpha1, alpha2,
                    b_Ai,b_Bi,b_Aij0,b_Aij1,b_Bij0,b_Bij1,b_ri){
    sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
            -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
            0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
    wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
            0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
            0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
    length_GK <- length(sk)
    
    log_S_t <- 0
    for (j in 1:length_GK) {
      sk_scaled_j = (sk[j]+1)*(t1-t0)/2+t0
      wk_scaled_j = wk[j]*(t1-t0)/2
      
      log_S_t <- log_S_t - wk_scaled_j*xi*sk_scaled_j^(xi-1)*exp((
        alpha1*betas1[7]+alpha1*b_Aij1+alpha2*betas2[7]+alpha2*b_Bij1)*sk_scaled_j)*
        exp(b_ri + c(1,X_surv) %*% gammas + alpha1*(c(1, X_surv) %*% betas1[1:6] + b_Ai+ b_Aij0)
            + alpha2*(c(1, X_surv) %*% betas2[1:6]+b_Bi+b_Bij0))
    }
    exp(log_S_t)
    
  }
  
  
  
  
  #### dynamic prediction (CI) #### 
  
  ## the specific ij of interest to make dynamic prediction
  pred_data_long_ij_target <- surv_data %>% filter(researchID == pred_id & tooth == pred_j)
  ## the corresponding cluster i
  pred_data_long_i <- long_data %>% filter(researchID == pred_id & followup_year <= pred.time)
  pred_data_surv_i <- surv_data %>% filter(researchID == pred_id)
  target_ij_index <- seq_len(nrow(pred_data_surv_i))[pred_data_surv_i$tooth == pred_j]
  
  
  obstimes_ij <- pred_data_surv_i$eventtime_year
  obstimes_ij <- ifelse(obstimes_ij <= pred.time, obstimes_ij, pred.time)
  obsstatus_ij <- pred_data_surv_i$status
  obsstatus_ij <- ifelse(pred_data_surv_i$eventtime_year <= pred.time, obsstatus_ij, 0)
  
  X_long_ij <- cbind(pred_data_long_i$sex, pred_data_long_i$std_age, pred_data_long_i$smoke, 
                     pred_data_long_i$molar, pred_data_long_i$lower)
  X_surv_ij <- cbind(pred_data_surv_i$sex, pred_data_surv_i$std_age, pred_data_surv_i$smoke, 
                     pred_data_surv_i$molar, pred_data_surv_i$lower)
  yA_ij <- log(pred_data_long_i$max_ppd)
  yB_ij <- pred_data_long_i$mobility 
  previous.time_ij <- pred_data_long_i$followup_year 
  n_i <- length(pred_data_surv_i$tooth)
  
  id_target <- pred_data_long_ij_target$id_tooth
  X_long_target <- c(pred_data_long_ij_target$sex, pred_data_long_ij_target$std_age, pred_data_long_ij_target$smoke, 
                     pred_data_long_ij_target$molar, pred_data_long_ij_target$lower)
  X_surv_target <- c(pred_data_long_ij_target$sex[1], pred_data_long_ij_target$std_age[1], pred_data_long_ij_target$smoke[1], 
                     pred_data_long_ij_target$molar[1], pred_data_long_ij_target$lower[1])
  #  Z <- data.frame(3)
  yA_target <- log(pred_data_long_ij_target$max_ppd)
  yB_target <- pred_data_long_ij_target$mobility
  previous.time_target <- pred_data_long_ij_target$followup_year  # observed times before prediction 
  
  
  
  out <- matrix(data=NA, nrow = MCMCnum, ncol = length(survTimes))
 
  mcmc <- draws_df_cluster
  mcmc_colname <- names(mcmc)
  mcmc_colname_no_b <- mcmc_colname[!grepl("^bi", mcmc_colname)]
  mcmc_no_b <- data.frame(mcmc[mcmc_colname_no_b])
  samples <- sample(length(mcmc$`beta1[1]`), MCMCnum)
  mcmc <- mcmc_no_b[samples,]
  
  SS_list <- vector("list", MCMCnum)
  CIF_list <- vector("list", MCMCnum)
  
  t0_test = Sys.time()
  for (m in 1:MCMCnum) { 
    # Step 1: extract parameter values
    betas1.new <- c(mcmc$beta1.1.[m],mcmc$beta1.2.[m],mcmc$beta1.3.[m],mcmc$beta1.4.[m],mcmc$beta1.5.[m],mcmc$beta1.6.[m],mcmc$beta1.7.[m])
    betas2.new <- c(mcmc$beta2.1.[m],mcmc$beta2.2.[m],mcmc$beta2.3.[m],mcmc$beta2.4.[m],mcmc$beta2.5.[m],mcmc$beta2.6.[m],mcmc$beta2.7.[m])
    gammas.new <- c(mcmc$gamma.1.[m], mcmc$gamma.2.[m], mcmc$gamma.3.[m],mcmc$gamma.4.[m], mcmc$gamma.5.[m], mcmc$gamma.6.[m])
    xi.new <- mcmc$xi[m]
    alpha1.new <- mcmc$alpha1[m]
    alpha2.new <- mcmc$alpha2[m]
    
    sigma_e.new <- mcmc$sigma_e[m]
    sigma_bi.new <- c(mcmc$sigma_bi.1.[m], mcmc$sigma_bi.2.[m])
    sigma_bij.new <- c(mcmc$sigma_bij.1.[m], mcmc$sigma_bij.2.[m], mcmc$sigma_bij.3.[m], mcmc$sigma_bij.4.[m])
    sigma_ri.new <- mcmc$sigma_ri[m]  
    R_bi <- mcmc$R_bi.1.2.[m]
    R_bij <- c(mcmc$R_bij.1.2.[m], mcmc$R_bij.1.3.[m], mcmc$R_bij.1.4.[m], 
               mcmc$R_bij.2.3.[m], mcmc$R_bij.2.4.[m], mcmc$R_bij.3.4.[m])
    
    var_ri.new = sigma_ri^2
    var_e.new = sigma_e^2
    
    length_bi <- length(sigma_bi)
    Cor_bi <- diag(1, length_bi)
    Cor_bi[lower.tri(Cor_bi)] <- R_bi
    Cor_bi <- Cor_bi + t(Cor_bi) - diag(diag(Cor_bi))  # Make symmetric
    var_bi.new <- diag(sigma_bi) %*% Cor_bi %*%  diag(sigma_bi)
    
    length_bij <- length(sigma_bij)
    Cor_bij <- diag(1, length_bij)
    Cor_bij[lower.tri(Cor_bij)] <- R_bij
    Cor_bij <- Cor_bij + t(Cor_bij) - diag(diag(Cor_bij))  # Make symmetric
    var_bij.new <- diag(sigma_bij) %*% Cor_bij %*%  diag(sigma_bij)
    
    L_bi_matrix.new <- t(chol(var_bi.new))
    L_bij_matrix.new <- t(chol(var_bij.new))
    
    
    # Step 2: simulate new random effects values
    
    # Data for the Stan model
    X1 <- X_surv_ij[,1]
    X2 <- X_surv_ij[,2]
    X3 <- X_surv_ij[,3]
    X4 <- X_surv_ij[,4]
    X5 <- X_surv_ij[,5]
    
    y1 <- yA_ij        # longitudinal outcomes
    y2 <- yB_ij          # longitudinal outcomes
    IDtooth <- as.numeric(factor(pred_data_long_i$id_tooth, levels = unique(pred_data_long_i$id_tooth)))
    
    IDtooth_s <- unique(IDtooth)
    nteeth <- length(X1)                     # total number of teeth
    
    status <- obsstatus_ij            # vital status (1 = dead, 0 = alive)
    times <- obstimes_ij           # times to event
    longdata_obstime <- previous.time_ij       # visit times for repeated observations
    longdata_row <- length(y1)                   # total number of longitudinal outcomes
    
    
    
    
    sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
            -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
            0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
    wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
            0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
            0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
    length_GK <- length(sk)
    
    
    DP_datalist <- list(
      N=longdata_row, nteeth=nteeth, y1=y1, y2=y2, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, IDtooth=IDtooth, 
      IDtooth_s=IDtooth_s, visits=longdata_obstime, times=times, status=status, 
      sk=sk, wk=wk, length_GK=length_GK, 
      beta1=betas1.new, beta2=betas2.new, gamma=gammas.new, xi=xi.new, 
      alpha1=alpha1.new, alpha2=alpha2.new, 
      sigma_e=sigma_e.new, sigma_ri=sigma_ri.new, 
      sigma_bi=sigma_bi.new, 
      sigma_bij=sigma_bij.new,
      L_bi=L_bi_matrix.new, 
      L_bij=L_bij_matrix.new)
    
    
    
    
    
    t0 <- Sys.time()
    fit_DP_cluster <- analysis_DP_mod$sample(data = DP_datalist,
                                             seed = 1234,
                                             chains = 4,
                                             parallel_chains = 4,
                                             iter_warmup = 200,
                                             iter_sampling = 200)
    t1 <- Sys.time()
    DP_runtime <- t1-t0
    # 10 sec for one person
    
    draws_DP_cluster0 <- fit_DP_cluster$draws(format = "df")
    
    samples_b <- sample(nrow(draws_DP_cluster0), MCMCnum_b)
    b_Ai.new <- c(draws_DP_cluster0["bi[1]"][samples_b,])[[1]]
    b_Bi.new <- c(draws_DP_cluster0["bi[2]"][samples_b,])[[1]]
    b_Aij0.new <- draws_DP_cluster0[,4:(3+nteeth)]
    b_Aij1.new <- draws_DP_cluster0[,(3+nteeth+1):(3+2*nteeth)]
    b_Bij0.new <- draws_DP_cluster0[,(3+2*nteeth+1):(3+3*nteeth)]
    b_Bij1.new <- draws_DP_cluster0[,(3+3*nteeth+1):(3+4*nteeth)]
    b_ri.new <- c(draws_DP_cluster0["ri"][samples_b,])[[1]]
    
    
    b_Aij0.new_target_ij <- c(b_Aij0.new[,target_ij_index][samples_b,])[[1]]
    b_Aij1.new_target_ij <- c(b_Aij1.new[,target_ij_index][samples_b,])[[1]]
    b_Bij0.new_target_ij <- c(b_Bij0.new[,target_ij_index][samples_b,])[[1]]
    b_Bij1.new_target_ij <- c(b_Bij1.new[,target_ij_index][samples_b,])[[1]]
    
    SS <- c()
    CIF <- c()
    for (mb in 1:MCMCnum_b) {
      # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)
      SS_mb <- S_u_t(survTimes,pred.time,X_surv_target, betas1.new, betas2.new, gammas.new, xi.new, 
                     alpha1.new, alpha2.new, 
                     b_Ai.new[mb], b_Bi.new[mb], 
                     b_Aij0.new_target_ij[mb], b_Aij1.new_target_ij[mb], 
                     b_Bij0.new_target_ij[mb], b_Bij1.new_target_ij[mb], b_ri.new[mb])
      CIF_mb <- 1 - SS_mb
      
      
      SS <- rbind(SS, SS_mb)
      CIF <- rbind(CIF, CIF_mb)
    }
    
    SS_list[[m]] = SS
    CIF_list[[m]] <- CIF
  }
  Sys.time() - t0_test
  
  SS_mat <- do.call(rbind, SS_list)
  CIF_mat <- do.call(rbind, CIF_list)
  
  pred.cil_SS<- c()
  pred.ciu_SS <- c()
  pred.median_SS <- c()
  pred.mean_SS <- c()
  pred.cil_CIF <- c()
  pred.ciu_CIF <- c()
  pred.median_CIF <- c()
  
  for (i in 1:length(survTimes)) {
    pred.cil_SS[i] <- quantile(SS_mat[,i], probs = CI.levels[1])
    pred.ciu_SS[i] <- quantile(SS_mat[,i], probs = CI.levels[2])
    pred.median_SS[i] <- median(SS_mat[,i])
    pred.mean_SS[i] <- mean(SS_mat[,i])
    pred.cil_CIF[i] <- quantile(CIF_mat[,i], probs = CI.levels[1])
    pred.ciu_CIF[i] <- quantile(CIF_mat[,i], probs = CI.levels[2])
    pred.median_CIF[i] <- median(CIF_mat[,i])
  }
  
  result_200_200 <- cbind(
    pred.cil_SS,
    pred.ciu_SS,
    pred.median_SS,
    pred.mean_SS,
    pred.cil_CIF,
    pred.ciu_CIF,
    pred.median_CIF
  )
  
  #### return output ####
  result <- cbind(
    times = survTimes,
    "Mean" = pred.mean_SS,
    "Median" = pred.median_SS,
    "Lower" = pred.cil_SS,
    "Upper" = pred.ciu_SS
  )
  
  return(result)
}



